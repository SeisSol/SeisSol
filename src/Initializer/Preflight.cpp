// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Preflight.h"

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Filesystem.h"
#include "Common/Json.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/Parameters/ParameterOverrides.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Monitoring/RunSummary.h"

#include <algorithm>
#include <cctype>
#include <csignal>
#include <cstddef>
#include <exception>
#include <fcntl.h>
#include <fty/fty.hpp>
#include <optional>
#include <signal.h>
#include <string>
#include <string_view>
#include <system_error>
#include <unistd.h>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yaml-cpp/yaml.h>

#ifdef USE_HDF
#include <hdf5.h>
#endif

namespace seissol::initializer {

namespace {

constexpr std::string_view Schema = "seissol-check/1";

// What an error that stops the check leaves behind: the report as it stood, and where it goes.
// logError() raises SIGTRAP before it aborts; a signal handler can only write what is there
// already, so the report is rendered ahead, at every step of the check.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::string stoppedReport;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::string stoppedReportPath;

void stopCheck(int /*signal*/) {
  if (!stoppedReportPath.empty()) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
    const int file = ::open(stoppedReportPath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (file >= 0) {
      std::size_t written = 0;
      while (written < stoppedReport.size()) {
        const auto result =
            ::write(file, stoppedReport.data() + written, stoppedReport.size() - written);
        if (result <= 0) {
          break;
        }
        written += static_cast<std::size_t>(result);
      }
      ::close(file);
    }
  }
  ::_exit(1);
}

class Report {
  public:
  explicit Report(std::string path) : path_(std::move(path)) {
    root_["schema"] = std::string(Schema);
    root_["status"] = "running";
    root_["problems"] = json::Value::array();
    stoppedReportPath = path_;
  }

  json::Value& operator[](std::string_view key) { return root_[key]; }

  //! @brief Something that keeps the run from starting.
  void problem(const std::string& text) {
    logWarning() << text.c_str();
    root_["problems"].push(text);
  }

  //! @brief Starts the step @p name , which an error may stop.
  void step(std::string_view name) {
    json::Value stopped = root_;
    stopped["status"] = "failed";
    stopped["stoppedAt"] = std::string(name);
    stoppedReport = stopped.dump() + "\n";
  }

  //! @brief Writes the report, if asked to, and returns the exit code.
  int finish() {
    const bool passed = root_["problems"].size() == 0;
    root_["status"] = passed ? "passed" : "failed";
    if (!path_.empty() && !writeFileAtomically(path_, root_.dump() + "\n")) {
      logWarning() << "Could not write the report of the check to" << path_;
    }
    const auto problems = root_["problems"].size();
    if (passed) {
      logInfo() << "The check passed.";
    } else {
      logInfo() << "The check failed:" << problems << (problems == 1 ? "problem." : "problems.");
    }
    return passed ? 0 : 1;
  }

  private:
  std::string path_;
  json::Value root_;
};

std::string joined(const std::vector<std::string>& parts) {
  std::string text;
  for (const auto& part : parts) {
    text += (text.empty() ? "" : ".") + part;
  }
  return text;
}

/**
 * @brief Whether a run opens the file that the parameter @p field names: some files are read only
 * by a part of SeisSol that the parameters may switch off.
 */
bool runOpens(const std::vector<std::string>& field,
              const parameters::SeisSolParameters& parameters) {
  const auto name = joined(field);
  const auto& faultOutput = parameters.drParameters.outputPointType;
  if (name == "output.rfilename") {
    return parameters.output.receiverParameters.enabled;
  }
  if (name == "pickpoint.ppfilename") {
    return parameters.drParameters.isDynamicRuptureEnabled &&
           (faultOutput == parameters::OutputType::AtPickpoint ||
            faultOutput == parameters::OutputType::AtPickpointAndElementwise);
  }
  if (name.starts_with("equations.plasticityfilename")) {
    return parameters.model.plasticity;
  }
  if (name == "inicondition.filename") {
    return parameters.initialization.type == parameters::InitializationType::Easi;
  }
  return true;
}

bool isYamlFile(const std::string& path) {
  auto extension = std::filesystem::path(path).extension().string();
  std::ranges::transform(extension, extension.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return extension == ".yaml" || extension == ".yml";
}

json::Value meshSection(const parameters::MeshParameters& mesh, Report& report) {
  json::Value section;
  section["file"] = mesh.meshFileName;
  if (mesh.meshFormat == parameters::MeshFormat::CubeGenerator) {
    section["format"] = "cubeGenerator";
    return section;
  }
  section["format"] = "puml";
#ifdef USE_HDF
  std::error_code error;
  if (!std::filesystem::is_regular_file(mesh.meshFileName, error)) {
    // which the check of the files has reported already
    return section;
  }
  // HDF5 would print its error stack for every dataset that is not there
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  const hid_t file = H5Fopen(mesh.meshFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    report.problem("The mesh file " + mesh.meshFileName + " is not an HDF5 file.");
    return section;
  }
  // the number of rows of a dataset, if it is there
  const auto rows = [&](const char* name) -> std::optional<std::size_t> {
    if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
      return {};
    }
    const hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dataset < 0) {
      return {};
    }
    std::optional<std::size_t> count;
    const hid_t space = H5Dget_space(dataset);
    const int dimensions = H5Sget_simple_extent_ndims(space);
    if (dimensions >= 1) {
      std::vector<hsize_t> extent(static_cast<std::size_t>(dimensions));
      H5Sget_simple_extent_dims(space, extent.data(), nullptr);
      count = static_cast<std::size_t>(extent[0]);
    }
    H5Sclose(space);
    H5Dclose(dataset);
    return count;
  };
  for (const auto* name : {"/connect", "/geometry", "/group", "/boundary"}) {
    if (!rows(name).has_value()) {
      report.problem("The mesh file " + mesh.meshFileName + " has no dataset " + name + ".");
    }
  }
  const auto cells = rows("/connect");
  const auto vertices = rows("/geometry");
  section["cells"] = cells.has_value() ? json::Value(cells.value()) : json::Value();
  section["vertices"] = vertices.has_value() ? json::Value(vertices.value()) : json::Value();
  H5Fclose(file);
  if (cells.has_value() && vertices.has_value()) {
    logInfo() << "The mesh has" << cells.value() << "cells and" << vertices.value() << "vertices.";
  }
#else
  report.problem("This build of SeisSol cannot read the PUML mesh " + mesh.meshFileName +
                 ", for it was built without HDF5.");
#endif
  return section;
}

json::Value outputSection(const parameters::SeisSolParameters& parameters, Report& report) {
  json::Value section;
  const auto& output = parameters.output;
  section["prefix"] = output.prefix;

  std::error_code error;
  const auto directory = std::filesystem::absolute(output.prefix, error).parent_path();
  const bool exists = std::filesystem::is_directory(directory, error);
  section["directory"] = directory.string();
  section["exists"] = exists;
  // a run creates the directory, with its parents, in the closest one that exists
  auto existing = directory;
  while (!existing.empty() && existing != existing.root_path() &&
         !std::filesystem::exists(existing, error)) {
    existing = existing.parent_path();
  }
  const bool writable = std::filesystem::is_directory(existing, error) &&
                        ::access(existing.string().c_str(), W_OK | X_OK) == 0;
  section["writable"] = writable;
  if (!writable) {
    report.problem(exists ? "The output directory " + directory.string() + " is not writable."
                          : "The output directory " + directory.string() +
                                " cannot be created in " + existing.string() + ".");
  }

  auto& enabled = section["enabled"] = json::Value::array();
  if (output.waveFieldParameters.enabled) {
    enabled.push("wavefield");
  }
  if (output.freeSurfaceParameters.enabled) {
    enabled.push("freeSurface");
  }
  if (output.receiverParameters.enabled) {
    enabled.push("receivers");
  }
  const auto faultOutput = parameters.drParameters.outputPointType;
  if (parameters.drParameters.isDynamicRuptureEnabled) {
    if (faultOutput == parameters::OutputType::AtPickpoint ||
        faultOutput == parameters::OutputType::AtPickpointAndElementwise) {
      enabled.push("faultReceivers");
    }
    if (faultOutput == parameters::OutputType::Elementwise ||
        faultOutput == parameters::OutputType::AtPickpointAndElementwise) {
      enabled.push("faultElementwise");
    }
  }
  if (output.energyParameters.enabled) {
    enabled.push("energy");
  }
  if (output.checkpointParameters.enabled) {
    enabled.push("checkpoint");
  }
  return section;
}

} // namespace

std::optional<std::string> extractCheckOption(int& argc, char** argv) {
  constexpr std::string_view Option = "--check";
  std::optional<std::string> report;
  int kept = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string_view argument(argv[i]);
    if (argument == Option) {
      report = std::string();
    } else if (argument.substr(0, Option.size() + 1) == std::string(Option) + "=") {
      report = std::string(argument.substr(Option.size() + 1));
    } else {
      argv[kept] = argv[i];
      ++kept;
    }
  }
  argc = kept;
  argv[argc] = nullptr;
  return report;
}

int preflight(const std::string& parameterFile,
              const std::vector<std::string>& overrides,
              const std::optional<std::string>& checkpoint,
              const std::string& reportPath) {
  Report report(reportPath);
  report["build"] = monitoring::buildInformation();
  {
    auto& invocation = report["invocation"];
    std::error_code error;
    const auto absolute = std::filesystem::absolute(parameterFile, error);
    invocation["parameterFile"] = error ? parameterFile : absolute.string();
    const auto workingDirectory = std::filesystem::current_path(error);
    invocation["workingDirectory"] = error ? std::string() : workingDirectory.string();
    invocation["overrides"] = json::Value::object();
    invocation["checkpoint"] =
        checkpoint.has_value() ? json::Value(checkpoint.value()) : json::Value();
  }
  // an error ends the check through the SIGTRAP of the logger; SIGTRAP can be caught, so that
  // setting its handler cannot fail
  static_cast<void>(std::signal(SIGTRAP, stopCheck));

  // what a run reads, read as a run reads it; an error stops the check here
  report.step("parameters");
  YAML::Node node;
  try {
    node = fty::Loader<fty::AsLowercase>{}.load(parameterFile);
  } catch (const std::exception& error) {
    logError() << "Error while reading the parameter file:" << std::string(error.what());
  }
  std::vector<parameters::ParameterOverride> parsedOverrides;
  for (const auto& text : overrides) {
    const auto parsed = parameters::parseOverride(text);
    if (!parsed.has_value()) {
      logError() << "The override" << text << "is not of the form section.field=value.";
    }
    parsedOverrides.push_back(parsed.value_or(parameters::ParameterOverride{}));
    report["invocation"]["overrides"][parsedOverrides.back().name()] = parsedOverrides.back().value;
  }
  parameters::applyOverrides(node, parsedOverrides);
  parameters::ParameterReader reader(node, parameterFile, false);
  const auto seissolParameters = [&] {
    try {
      return parameters::readSeisSolParameters(&reader);
    } catch (const std::exception& error) {
      logError() << "Uncaught error while parsing the parameter file:" << std::string(error.what());
    }
    return parameters::SeisSolParameters();
  }();
  reader.warnUnknown();
  for (const auto& entry : parameters::unknownOverrides(reader, parsedOverrides)) {
    logError() << "The override" << entry.name()
               << "names a parameter that SeisSol does not read; check its spelling.";
  }
  auto& configurations = report["configurations"] = json::Value::array();
  for (const auto config : seissolParameters.model.configs()) {
    configurations.push(configName(configValue(config)));
  }
  report["endTime"] = seissolParameters.timeStepping.endTime;
  report["unknownFields"] = json::Value::arrayOf(reader.unknownFields());

  // the files the parameters name
  report.step("files");
  const auto paths = reader.paths();
  auto& files = report["files"] = json::Value::array();
  for (const auto& entry : paths) {
    if (entry.path.empty()) {
      continue;
    }
    const auto name = joined(entry.field);
    std::error_code error;
    const bool exists = std::filesystem::exists(entry.path, error);
    const bool opened = runOpens(entry.field, seissolParameters);
    auto& file = files.push(json::Value::object());
    file["parameter"] = name;
    file["path"] = entry.path;
    file["exists"] = exists;
    file["used"] = opened;
    if (!exists && opened) {
      report.problem("The file " + entry.path + " that " + name + " names does not exist.");
    }
  }

  // the YAML files among them, the easi models, parse; what they hold is checked by the run
  report.step("yaml");
  for (const auto& entry : paths) {
    std::error_code error;
    if (entry.path.empty() || !isYamlFile(entry.path) ||
        !runOpens(entry.field, seissolParameters) ||
        !std::filesystem::is_regular_file(entry.path, error)) {
      continue;
    }
    try {
      YAML::LoadFile(entry.path);
    } catch (const std::exception& exception) {
      report.problem("The YAML file " + entry.path + " that " + joined(entry.field) +
                     " names does not parse: " + exception.what());
    }
  }

  report.step("mesh");
  report["mesh"] = meshSection(seissolParameters.mesh, report);

  if (checkpoint.has_value()) {
    report.step("checkpoint");
    std::error_code error;
    if (!std::filesystem::is_regular_file(checkpoint.value(), error)) {
      report.problem("The checkpoint " + checkpoint.value() + " does not exist.");
    }
  }

  report.step("output");
  report["output"] = outputSection(seissolParameters, report);

  static_cast<void>(std::signal(SIGTRAP, SIG_DFL));
  return report.finish();
}

} // namespace seissol::initializer
