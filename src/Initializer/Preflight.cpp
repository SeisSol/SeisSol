// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Preflight.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Filesystem.h"
#include "Common/Json.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "Equations/Datastructures.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/FaceMap.h"
#include "Initializer/InitProcedure/Internal/FaceTypeCheck.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/Parameters/ParameterOverrides.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Parameters/SourceParameters.h"
#include "Initializer/PointMapper.h"
#include "Initializer/PreflightMesh.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "Model/CommonDatastructures.h"
#include "Monitoring/FailureReport.h"
#include "Monitoring/RunSummary.h"
#include "ResultWriter/ReceiverWriter.h"
#include "SourceTerm/FSRMReader.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <easi/Component.h>
#include <easi/YAMLParser.h>
#include <exception>
#include <filesystem>
#include <fty/fty.hpp>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yaml-cpp/yaml.h>

// whether a directory is writable: POSIX has access(); Windows tries to create a file there
#ifdef _WIN32
#include <fstream>
#include <random>
#else
#include <unistd.h>
#endif

#ifdef USE_HDF
#include <hdf5.h>
#endif

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
#include "SourceTerm/NRF.h"
#include "SourceTerm/NRFReader.h"
#endif

namespace seissol::initializer {

namespace {

constexpr std::string_view Schema = "seissol-check/1";

//! @brief Sets the abort handler of the logger for as long as it exists; then the one before again.
class AbortHandlerScope {
  public:
  explicit AbortHandlerScope(utils::Logger::AbortHandler handler)
      : previous_(utils::Logger::setAbortHandler(std::move(handler))) {}
  ~AbortHandlerScope() { utils::Logger::setAbortHandler(std::move(previous_)); }
  AbortHandlerScope(const AbortHandlerScope&) = delete;
  AbortHandlerScope(AbortHandlerScope&&) = delete;
  AbortHandlerScope& operator=(const AbortHandlerScope&) = delete;
  AbortHandlerScope& operator=(AbortHandlerScope&&) = delete;

  private:
  utils::Logger::AbortHandler previous_;
};

class Report {
  public:
  explicit Report(std::string path) : path_(std::move(path)) {
    root_["schema"] = std::string(Schema);
    root_["status"] = "running";
    root_["problems"] = json::Value::array();
    root_["warnings"] = json::Value::array();
  }

  json::Value& operator[](std::string_view key) { return root_[key]; }

  //! @brief Something that keeps the run from starting.
  void problem(const std::string& text) {
    logWarning() << text.c_str();
    root_["problems"].push(text);
    setFailure();
  }

  //! @brief Something that does not keep the run from starting, but may not be as intended.
  void warning(const std::string& text) {
    logWarning() << text.c_str();
    root_["warnings"].push(text);
    setFailure();
  }

  //! @brief Starts the step @p name , which an error or a crash may stop.
  void step(std::string_view name) {
    step_ = std::string(name);
    setFailure();
  }

  /**
   * @brief Writes the report of a check that the error of @p record stopped: the step it stopped,
   * and the message and the place of the error.
   */
  void stop(const utils::LogRecord& record) {
    // with what the step has found so far
    setFailure();
    monitoring::writeFailureReport(monitoring::failureOf(record));
  }

  //! @brief Writes the report, if asked to, and returns the exit code.
  int finish() {
    monitoring::clearFailureReports();
    const bool passed = root_["problems"].size() == 0;
    root_["status"] = passed ? "passed" : "failed";
    write();
    const auto problems = root_["problems"].size();
    const auto warnings = root_["warnings"].size();
    if (passed && warnings > 0) {
      logInfo() << "The check passed, with" << warnings
                << (warnings == 1 ? "warning." : "warnings.");
    } else if (passed) {
      logInfo() << "The check passed.";
    } else {
      logInfo() << "The check failed:" << problems << (problems == 1 ? "problem." : "problems.");
    }
    return passed ? 0 : 1;
  }

  private:
  void write() {
    if (!path_.empty() && !writeFileAtomically(path_, root_.dump() + "\n")) {
      logWarning() << "Could not write the report of the check to" << path_;
    }
  }

  //! @brief The report as it stands, with @p status and the placeholder in @p key after the step.
  [[nodiscard]] json::Value with(const char* status, const char* key) const {
    auto changed = json::Value::object();
    for (const auto& [name, value] : root_.asObject()) {
      changed[name] = name == "status" ? json::Value(status) : value;
      if (name == "status") {
        changed["stoppedAt"] = step_;
        changed[key] = monitoring::FailurePlaceholder;
      }
    }
    return changed;
  }

  /**
   * @brief Sets what a failure and a request to stop write: the report as it stands, failed with
   * the error, or terminated with the signal, each after the step it stopped.
   */
  void setFailure() {
    if (path_.empty()) {
      return;
    }
    monitoring::setFailureReport(path_, with("failed", "error"));
    monitoring::setTerminationReport(path_, with("terminated", "terminatedBy"));
  }

  std::string path_;
  json::Value root_;
  std::string step_;
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

//! @brief Whether files can be created in the directory @p path .
bool directoryWritable(const std::filesystem::path& path) {
#ifdef _WIN32
  // access does not know the permissions of Windows (access control lists), and its read-only
  // attribute does not keep anyone from creating files in a directory; so try it
  std::random_device seed;
  const auto probe = path / ("seissol-check-" + std::to_string(seed()));
  bool created = false;
  {
    const std::ofstream stream(probe);
    created = stream.good();
  }
  std::error_code error;
  std::filesystem::remove(probe, error);
  return created;
#else
  return ::access(path.string().c_str(), W_OK | X_OK) == 0;
#endif
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
  const bool writable =
      std::filesystem::is_directory(existing, error) && directoryWritable(existing);
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

//! @brief Whether the parameter @p name names an easi model, which the run reads with easi.
bool isEasiModel(const std::string& name) {
  return name == "equations.materialfilename" || name.starts_with("equations.plasticityfilename") ||
         name == "equations.boundaryfilename" || name.starts_with("dynamicrupture.modelfilename") ||
         name == "inicondition.filename";
}

/**
 * @brief Where easi looks for the file @p entry that the model file @p model names: beside the
 * model file, if it is there, or else as it is given (see YAMLParser::parse, getFileName).
 */
std::filesystem::path easiPath(const std::string& entry, const std::filesystem::path& model) {
  std::filesystem::path path(entry);
  if (path.is_relative()) {
    const auto beside = model.parent_path() / path;
    std::error_code error;
    if (std::filesystem::exists(beside, error)) {
      return beside;
    }
  }
  return path;
}

//! @brief What the files of an easi model refer to.
struct EasiReferences {
  //! whether a component reads a file with ASAGI
  bool asagi{false};
  //! the files of the model, the one named by the parameter and those it includes
  std::set<std::string> files;
};

/**
 * @brief Checks that the files the node @p node of the model file @p file refers to exist: those it
 * includes (!Include), whose nodes it checks as well, and those that its components read (!ASAGI,
 * !SCECFile). @p model is the file the parameter names.
 */
void checkEasiReferences(const YAML::Node& node,
                         const std::filesystem::path& file,
                         const std::string& model,
                         EasiReferences& references,
                         Report& report) {
  std::error_code error;
  const auto& tag = node.Tag();
  if (tag == "!Include" && node.IsScalar()) {
    const auto included = easiPath(node.as<std::string>(), file);
    if (!std::filesystem::is_regular_file(included, error)) {
      report.problem("The easi model " + model + " includes " + node.as<std::string>() + " (in " +
                     file.string() + "), which does not exist.");
      return;
    }
    const auto canonical = std::filesystem::weakly_canonical(included, error).string();
    if (references.files.insert(canonical).second) {
      try {
        checkEasiReferences(YAML::LoadFile(included.string()), included, model, references, report);
      } catch (const std::exception& exception) {
        report.problem("The file " + included.string() + " that the easi model " + model +
                       " includes does not parse: " + exception.what());
      }
    }
    return;
  }
  if ((tag == "!ASAGI" || tag == "!SCECFile") && node.IsMap() && node["file"] &&
      node["file"].IsScalar()) {
    references.asagi = references.asagi || tag == "!ASAGI";
    const auto entry = node["file"].as<std::string>();
    if (!std::filesystem::is_regular_file(easiPath(entry, file), error)) {
      report.problem("The easi model " + model + " reads " + entry + " (with " + tag + ", in " +
                     file.string() + "), which does not exist.");
    }
  }
  if (node.IsMap()) {
    for (const auto& entry : node) {
      checkEasiReferences(entry.second, file, model, references, report);
    }
  } else if (node.IsSequence()) {
    for (const auto& entry : node) {
      checkEasiReferences(entry, file, model, references, report);
    }
  }
}

/**
 * @brief Parses the easi model @p model as the run does, if it can be without MPI, and checks that
 * a material file supplies the parameters of the materials of the run.
 */
json::Value checkEasiModel(const std::string& parameter,
                           const std::string& modelFile,
                           const parameters::SeisSolParameters& parameters,
                           Report& report,
                           std::unique_ptr<easi::Component>* parsed) {
  json::Value section;
  section["parameter"] = parameter;
  section["file"] = modelFile;
  EasiReferences references;
  std::error_code error;
  references.files.insert(std::filesystem::weakly_canonical(modelFile, error).string());
  try {
    checkEasiReferences(YAML::LoadFile(modelFile), modelFile, modelFile, references, report);
  } catch (const std::exception& /*exception*/) {
    // which the check of the YAML files reported already
    section["parsed"] = false;
    return section;
  }
  section["asagi"] = references.asagi;
  if (references.asagi) {
    // ASAGI reads its files with MPI, which the check does not have
    section["parsed"] = false;
    return section;
  }

  const bool timeDependent =
      parameter == "inicondition.filename" && parameters.initialization.hasTime;
  std::unique_ptr<easi::Component> component;
  try {
#if defined(EASI_VERSION_MAJOR) && (EASI_VERSION_MINOR >= 5 || EASI_VERSION_MAJOR > 1)
    easi::YAMLParser parser(timeDependent ? std::set<std::string>{"t", "x", "y", "z"}
                                          : std::set<std::string>{"x", "y", "z"});
#else
    // (an easi too old for time-dependent conditions stops the run, see EasiLoader)
    easi::YAMLParser parser(3);
#endif
    component.reset(parser.parse(modelFile));
  } catch (const std::exception& exception) {
    report.problem("The easi model " + modelFile + " that " + parameter +
                   " names does not parse: " + exception.what());
    section["parsed"] = false;
    return section;
  }
  section["parsed"] = true;
  const auto supplied = component->suppliedParameters();
  section["supplied"] =
      json::Value::arrayOf(std::vector<std::string>(supplied.begin(), supplied.end()));

  // the parameters the materials of the configurations of the run need
  const auto needs = [&](const std::vector<std::string>& missing, const std::string& what) {
    if (missing.empty()) {
      return;
    }
    std::string list;
    for (const auto& name : missing) {
      list += (list.empty() ? "" : ", ") + name;
    }
    report.problem("The easi model " + modelFile + " that " + parameter +
                   " names does not supply " + list + ", which " + what + " needs.");
  };
  const bool material = parameter == "equations.materialfilename";
  const bool plasticity = parameter.starts_with("equations.plasticityfilename") ||
                          (material && parameters.model.plasticity &&
                           std::ranges::find(parameters.model.plasticityFileNames, modelFile) !=
                               parameters.model.plasticityFileNames.end());
  for (const auto config : parameters.model.configs()) {
    dispatchConfig(config, [&](auto cfg) {
      using MaterialT = model::MaterialOf<decltype(cfg)>;
      const auto name = configName(configValue(config));
      if (material) {
        needs(missingMaterialParameters<MaterialT>(supplied),
              "the material " + std::string(MaterialT::Text) + " of the configuration " + name);
      }
      if (plasticity && parameters.model.plasticity && MaterialT::SupportsPlasticity) {
        needs(missingMaterialParameters<model::Plasticity>(supplied),
              "the plasticity of the configuration " + name);
      }
    });
  }
  if (parsed != nullptr) {
    *parsed = std::move(component);
  }
  return section;
}

//! @brief What the evaluation of the material model for the cells found.
struct ModelEvaluation {
  std::size_t cells{0};
  std::size_t invalidCells{0};
  std::optional<std::size_t> firstInvalidCell;
  double minTimestep{std::numeric_limits<double>::infinity()};
  std::size_t minTimestepCell{0};
  double maxTimestep{0};
};

/**
 * @brief Evaluates the material model @p material for the cells of @p chunk , as the run does for
 * their time steps (see computeTimesteps), and takes the time steps and the cells without a wave
 * speed into @p evaluation .
 */
void evaluateMaterial(const meshscan::CellChunk& chunk,
                      easi::Component& material,
                      const parameters::SeisSolParameters& parameters,
                      ModelEvaluation& evaluation) {
  std::vector<std::size_t> valid;
  for (std::size_t i = 0; i < chunk.vertices.size(); ++i) {
    if (chunk.valid[i] != 0) {
      valid.push_back(i);
    }
  }
  const auto all = CellToVertexArray::fromVectors(chunk.vertices, chunk.groups);
  const auto& model = parameters.model;
  const auto configs = model.configs();
  for (const auto config : configs) {
    std::vector<std::size_t> cells;
    for (const auto cell : valid) {
      if (configs.size() == 1 || model.configOfGroup(chunk.groups[cell]) == config) {
        cells.push_back(cell);
      }
    }
    if (cells.empty()) {
      continue;
    }
    const auto cellsOfConfig = CellToVertexArray::subset(all, cells);
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using MaterialT = model::MaterialOf<Cfg>;
      const auto queryGen = getBestQueryGenerator<MaterialT>(
          model.useCellHomogenizedMaterial, cellsOfConfig, Cfg::ConvergenceOrder);
      std::vector<MaterialT> materials(cellsOfConfig.size);
      MaterialParameterDB<MaterialT> database;
      database.setMaterialVector(&materials);
      database.evaluate(material, *queryGen);

      for (std::size_t i = 0; i < cellsOfConfig.size; ++i) {
        const auto cell = chunk.first + cells[i];
        const double speed = materials[i].getMaxWaveSpeed();
        ++evaluation.cells;
        if (!(speed > 0) || !std::isfinite(speed)) {
          ++evaluation.invalidCells;
          if (!evaluation.firstInvalidCell.has_value()) {
            evaluation.firstInvalidCell = cell;
          }
          continue;
        }
        const auto cellMaxTimestep =
            std::min(materials[i].maximumTimestep(), parameters.timeStepping.maxTimestepWidth);
        const double timestep = computeCellTimestep(cellsOfConfig.elementCoordinates(i),
                                                    speed,
                                                    parameters.timeStepping.cfl,
                                                    cellMaxTimestep,
                                                    Cfg::ConvergenceOrder);
        if (timestep < evaluation.minTimestep) {
          evaluation.minTimestep = timestep;
          evaluation.minTimestepCell = cell;
        }
        evaluation.maxTimestep = std::max(evaluation.maxTimestep, timestep);
      }
    });
  }
}

//! @brief The number @p value as text, as short as it can be to read back the same.
std::string shortest(double value) {
  std::array<char, 32> buffer{};
  const auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value);
  return {buffer.data(), result.ptr};
}

//! @brief The point @p point as text, (x, y, z).
std::string formatted(const Eigen::Vector3d& point) {
  return "(" + shortest(point(0)) + ", " + shortest(point(1)) + ", " + shortest(point(2)) + ")";
}

json::Value pointValue(const Eigen::Vector3d& point) {
  auto value = json::Value::array();
  for (int d = 0; d < 3; ++d) {
    value.push(point(d));
  }
  return value;
}

//! @brief What the scan tells of why the point @p point lies in no cell; empty if nothing.
std::string missHint(const Eigen::Vector3d& point,
                     const meshscan::CellPoints::Found& found,
                     const meshscan::MeshScan& scan) {
  std::ostringstream text;
  if (found.surface.has_value() && point(2) > found.surface.value()) {
    text << ": it lies " << shortest(point(2) - found.surface.value())
         << " above the free surface, which is at z = " << shortest(found.surface.value())
         << " there";
    return text.str();
  }
  constexpr std::array<char, 3> Axes{'x', 'y', 'z'};
  bool outside = false;
  for (int d = 0; d < 3; ++d) {
    if (point(d) < scan.min(d) || point(d) > scan.max(d)) {
      text << (outside ? ", " : ": it lies outside the mesh, ") << Axes[d] << " = "
           << shortest(point(d)) << " not in [" << shortest(scan.min(d)) << ", "
           << shortest(scan.max(d)) << "]";
      outside = true;
    }
  }
  return text.str();
}

//! the points a problem names one by one at most; the report lists all of them
constexpr std::size_t NamedPoints = 20;

/**
 * @brief Reports the points of @p points that lie in no cell, as the run locates them (see
 * findUniqueMeshIds), which stops the run.
 */
json::Value cellPointsSection(const char* what,
                              const std::string& file,
                              const meshscan::CellPoints& points,
                              const meshscan::MeshScan& scan,
                              Report& report) {
  json::Value section;
  section["file"] = file;
  section["count"] = points.points.size();
  auto& missing = section["missing"] = json::Value::array();
  for (std::size_t i = 0; i < points.points.size(); ++i) {
    const auto& found = points.found[i];
    if (found.cell.has_value()) {
      continue;
    }
    const auto hint = missHint(points.points[i], found, scan);
    auto& entry = missing.push(json::Value::object());
    entry["index"] = i;
    entry["point"] = pointValue(points.points[i]);
    entry["surface"] =
        found.surface.has_value() ? json::Value(found.surface.value()) : json::Value();
    entry["hint"] = hint;
    if (missing.size() <= NamedPoints) {
      std::ostringstream text;
      text << "The " << what << " " << i << " at " << formatted(points.points[i]) << " of " << file
           << " lies in no cell of the mesh" << hint << ".";
      report.problem(text.str());
    }
  }
  if (missing.size() > NamedPoints) {
    report.problem("Another " + std::to_string(missing.size() - NamedPoints) + " " + what +
                   "s of " + file + " lie in no cell of the mesh (see the report).");
  }
  section["found"] = points.points.size() - missing.size();
  return section;
}

/**
 * @brief Reports the on-fault receivers that lie on no fault face, which stops the run, and those
 * that the run moves farther onto the fault than the face is large.
 */
json::Value faultPointsSection(const std::string& file,
                               const meshscan::FaultPoints& points,
                               const meshscan::MeshScan& scan,
                               Report& report) {
  json::Value section;
  section["file"] = file;
  section["count"] = points.points.size();
  // (both members first: adding one may move the other)
  section["missing"] = json::Value::array();
  section["moved"] = json::Value::array();
  auto& missing = section["missing"];
  auto& moved = section["moved"];
  for (std::size_t i = 0; i < points.points.size(); ++i) {
    const Eigen::Vector3d point(points.points[i][0], points.points[i][1], points.points[i][2]);
    const auto& found = points.found[i];
    if (!found.has_value()) {
      auto& entry = missing.push(json::Value::object());
      entry["index"] = i;
      entry["point"] = pointValue(point);
      if (missing.size() <= NamedPoints) {
        report.problem("The on-fault receiver " + std::to_string(i) + " at " + formatted(point) +
                       " of " + file + " lies on no fault face" +
                       (scan.faultFaces == 0
                            ? std::string(": the mesh has no fault faces")
                            : std::string(": no face of the fault lies across it, along "
                                          "its normal; it lies beyond the fault")) +
                       ".");
      }
    } else if (found->distance > found->faceSize) {
      auto& entry = moved.push(json::Value::object());
      entry["index"] = i;
      entry["point"] = pointValue(point);
      entry["distance"] = found->distance;
      if (moved.size() <= NamedPoints) {
        std::ostringstream text;
        text << "The on-fault receiver " << i << " at " << formatted(point) << " of " << file
             << " lies " << shortest(found->distance)
             << " away from the fault, farther than the fault face is large; the run moves it "
                "onto that face.";
        report.warning(text.str());
      }
    }
  }
  if (missing.size() > NamedPoints) {
    report.problem("Another " + std::to_string(missing.size() - NamedPoints) +
                   " on-fault receivers of " + file + " lie on no fault face (see the report).");
  }
  if (moved.size() > NamedPoints) {
    report.warning("Another " + std::to_string(moved.size() - NamedPoints) +
                   " on-fault receivers of " + file +
                   " lie farther from the fault than its faces are large (see the report).");
  }
  section["found"] = points.points.size() - missing.size();
  return section;
}

/**
 * @brief Reports what the scan of the whole mesh found about its faces and cells, and adds it to
 * the @p section of the mesh.
 */
void checkScannedMesh(const parameters::SeisSolParameters& parameters,
                      const meshscan::MeshScan& scan,
                      json::Value& section,
                      Report& report) {
  auto& scanned = section["scan"];
  scanned["seconds"] = scan.seconds;
  scanned["bytes"] = scan.bytes;
  scanned["cells"] = scan.cells;
  scanned["vertices"] = scan.vertices;
  scanned["min"] = pointValue(scan.min);
  scanned["max"] = pointValue(scan.max);
  scanned["smallestInsphere"] = scan.smallestInsphere;
  scanned["smallestCell"] = scan.smallestCell;
  scanned["largestInsphere"] = scan.largestInsphere;
  scanned["faultFaces"] = scan.faultFaces;
  scanned["oneSidedFaultFaces"] = scan.oneSidedFaultFaces;
  scanned["faultFacesWithoutPlusSide"] = scan.faultFacesWithoutPlusSide;
  auto& tags = scanned["tags"] = json::Value::object();
  for (const auto& [tag, faces] : scan.tags) {
    tags[std::to_string(tag)] = faces;
  }
  auto& faceTypes = scanned["faceTypes"] = json::Value::object();
  for (const auto& [type, faces] : scan.faceTypes) {
    faceTypes[faceTypeName(type)] = faces;
  }
  auto& groups = scanned["groups"] = json::Value::object();
  for (const auto& [group, cells] : scan.groups) {
    groups[std::to_string(group)] = cells;
  }
  {
    std::ostringstream text;
    text << "The mesh spans " << formatted(scan.min) << " to " << formatted(scan.max)
         << "; the radii of the inspheres of its cells range from " << scan.smallestInsphere
         << " (the cell " << scan.smallestCell << " of the file) to " << scan.largestInsphere
         << ".";
    logInfo() << text.str().c_str();
  }

  if (scan.invalidCells > 0) {
    report.problem("The mesh has " + std::to_string(scan.invalidCells) +
                   " cells with vertices that its file does not have.");
  }
  if (scan.misorientedCells > 0) {
    report.problem("The mesh has " + std::to_string(scan.misorientedCells) +
                   " cells that the run rejects for their orientation, e.g. the cell " +
                   std::to_string(scan.firstMisorientedCell.value_or(0)) +
                   " of the file: cells without volume, or all cells if the scaling of the mesh "
                   "(meshnml.scalingmatrixx, y, z) mirrors it.");
  }

  for (const auto& [tag, faces] : scan.unknownTags) {
    report.problem(
        "The mesh has " + std::to_string(faces) + " faces with the boundary condition tag " +
        std::to_string(tag) + ", which the face map does not know" +
        (parameters.mesh.faceMapFile.has_value() ? " (" + parameters.mesh.faceMapFile.value() + ")"
                                                 : std::string(" (the default one)")) +
        "; the run stops at them.");
  }

  // what the cells of each configuration have to support
  std::vector<internal::FaceTypeCensus> census(builtConfigCount());
  for (const auto& [group, types] : scan.groupFaceTypes) {
    census[parameters.model.configOfGroup(group)].present |= types;
  }
  for (const auto& problem : internal::faceTypeProblems(census, parameters.initialization.type)) {
    report.problem(
        "The boundary conditions of the mesh need what this build cannot do: " + problem + ".");
  }

  for (const auto& [group, config] : parameters.model.groupConfigs) {
    if (!scan.groups.contains(group)) {
      report.warning("The configuration map (equations.configmap) gives the group " +
                     std::to_string(group) + " a configuration, but no cell of the mesh has it.");
    }
  }
  if (parameters.drParameters.isDynamicRuptureEnabled && scan.faultFaces == 0) {
    report.warning("Dynamic rupture is on, but the mesh has no fault faces.");
  }
  if (scan.oneSidedFaultFaces > 0) {
    report.problem("The mesh has " + std::to_string(scan.oneSidedFaultFaces) +
                   " fault faces that only one of their cells has as such: on the boundary of the "
                   "mesh, where the run stops at them as at broken faces, or where the other cell "
                   "gives the face another boundary condition tag.");
  }
  if (scan.faultFacesWithoutPlusSide > 0) {
    report.warning("For " + std::to_string(scan.faultFacesWithoutPlusSide) +
                   " fault faces, the reference point (dynamicrupture.xref, yref, zref) tells "
                   "neither side to be the plus side: it lies in their plane, or is normal to "
                   "them with refpointmethod = 1. Which side is which is arbitrary there.");
  }
  const auto freeSurfaces = [&] {
    std::size_t faces = 0;
    for (const auto type : {FaceType::FreeSurface, FaceType::FreeSurfaceGravity}) {
      const auto found = scan.faceTypes.find(type);
      faces += found == scan.faceTypes.end() ? 0 : found->second;
    }
    return faces;
  }();
  if (parameters.output.freeSurfaceParameters.enabled && freeSurfaces == 0) {
    report.warning("The free-surface output is on, but the mesh has no faces of the free surface.");
  }
}

/**
 * @brief Reads the whole mesh, and checks its faces and cells and where the points of the run lie
 * in it (see meshscan::scanMesh).
 */
void checkMeshContent(const parameters::SeisSolParameters& parameters,
                      CheckDepth depth,
                      easi::Component* material,
                      Report& report) {
  const auto& output = parameters.output;
  const auto& rupture = parameters.drParameters;
  std::error_code error;

  // the points the run locates in the mesh, read as the run reads them
  report.step("points");
  meshscan::CellPoints receivers;
  receivers.tolerance = ReceiverTolerance;
  meshscan::CellPoints sources;
  sources.tolerance = SourceTolerance;
  meshscan::FaultPoints faultReceivers;
  const auto& receiverFile = output.receiverParameters.fileName;
  const bool withReceivers =
      output.receiverParameters.enabled && std::filesystem::is_regular_file(receiverFile, error);
  if (withReceivers) {
    try {
      receivers.points = writer::parseReceiverFile(receiverFile);
    } catch (const std::exception& exception) {
      report.problem("The receiver file " + receiverFile + " does not parse: " + exception.what());
    }
  }
  const auto pickpointFile = output.pickpointParameters.pickpointFileName.value_or("");
  const bool withFaultReceivers =
      rupture.isDynamicRuptureEnabled &&
      (rupture.outputPointType == parameters::OutputType::AtPickpoint ||
       rupture.outputPointType == parameters::OutputType::AtPickpointAndElementwise) &&
      std::filesystem::is_regular_file(pickpointFile, error);
  if (withFaultReceivers) {
    try {
      faultReceivers.points = dr::readPickpoints(pickpointFile);
    } catch (const std::exception& exception) {
      report.problem("The on-fault receiver file " + pickpointFile +
                     " does not parse: " + exception.what());
    }
  }
  const auto& sourceFile = parameters.source.fileName;
  const bool withSources = parameters.source.type != parameters::PointSourceType::None &&
                           std::filesystem::is_regular_file(sourceFile, error);
  if (withSources && parameters.source.type == parameters::PointSourceType::FsrmSource) {
    try {
      sourceterm::FSRMSource fsrm;
      fsrm.read(sourceFile);
      sources.points = fsrm.centers;
    } catch (const std::exception& exception) {
      report.problem("The point source file " + sourceFile +
                     " does not parse: " + exception.what());
    }
  }
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
  if (withSources && parameters.source.type == parameters::PointSourceType::NrfSource) {
    sourceterm::NRF nrf;
    sourceterm::readNRF(sourceFile.c_str(), nrf);
    sources.points = nrf.centres;
  }
#endif

  report.step("meshScan");
  auto& mesh = report["mesh"];
  const auto* cells = mesh.find("cells");
  const auto* vertices = mesh.find("vertices");
  if (parameters.mesh.meshFormat != parameters::MeshFormat::PUML || cells == nullptr ||
      cells->isNull() || vertices == nullptr || vertices->isNull()) {
    // where the mesh cannot be read, the check of its head reported why already
    return;
  }
  meshscan::MeshScanSettings settings;
  settings.file = parameters.mesh.meshFileName;
  settings.boundaryFormat = parameters.mesh.pumlBoundaryFormat;
  settings.displacement = parameters.mesh.displacement;
  settings.scaling = parameters.mesh.scaling;
  settings.faceMap = parameters.mesh.faceMapFile.has_value()
                         ? readFaceMap(parameters.mesh.faceMapFile.value())
                         : defaultFaceMap();
  settings.dynamicRupture = rupture.isDynamicRuptureEnabled;
  settings.referencePoint = rupture.referencePoint;
  settings.refPointMethod = rupture.refPointMethod;
  ModelEvaluation evaluation;
  const bool evaluate = depth == CheckDepth::Model && material != nullptr;
  if (evaluate) {
    // the material at the points of the cells, which the run homogenizes the material over, takes
    // about 3 KB per cell
    settings.chunk = std::size_t{1} << 16;
    settings.visitCells = [&](const meshscan::CellChunk& chunk) {
      evaluateMaterial(chunk, *material, parameters, evaluation);
    };
  }
  logInfo() << "Reading the whole mesh, with" << receivers.points.size() << "receivers,"
            << faultReceivers.points.size() << "on-fault receivers and" << sources.points.size()
            << "point sources" << (evaluate ? "; evaluating the material model." : ".");
  const auto scan = meshscan::scanMesh(settings, receivers, sources, faultReceivers);
  logInfo() << "Read the whole mesh in" << scan.seconds << "seconds.";

  checkScannedMesh(parameters, scan, mesh, report);
  auto& points = report["points"] = json::Value::object();
  if (withReceivers) {
    points["receivers"] = cellPointsSection("receiver", receiverFile, receivers, scan, report);
  }
  if (withFaultReceivers) {
    points["faultReceivers"] = faultPointsSection(pickpointFile, faultReceivers, scan, report);
  }
  if (withSources) {
    points["sources"] = cellPointsSection("point source", sourceFile, sources, scan, report);
  }

  if (depth == CheckDepth::Model) {
    auto& section = report["model"];
    section["file"] = parameters.model.materialFileName;
    section["evaluated"] = evaluate;
    if (!evaluate) {
      report.warning("The check does not evaluate the material model " +
                     parameters.model.materialFileName +
                     ": it reads files with ASAGI, which needs MPI, or it does not parse.");
      return;
    }
    section["cells"] = evaluation.cells;
    section["invalidCells"] = evaluation.invalidCells;
    if (evaluation.invalidCells > 0) {
      report.problem("The material model " + parameters.model.materialFileName + " gives " +
                     std::to_string(evaluation.invalidCells) +
                     " cells no positive and finite wave speed, e.g. the cell " +
                     std::to_string(evaluation.firstInvalidCell.value_or(0)) + " of the file.");
    }
    if (evaluation.minTimestep < std::numeric_limits<double>::infinity()) {
      section["minTimestep"] = evaluation.minTimestep;
      section["minTimestepCell"] = evaluation.minTimestepCell;
      section["maxTimestep"] = evaluation.maxTimestep;
      const double steps = parameters.timeStepping.endTime / evaluation.minTimestep;
      section["steps"] = steps;
      std::ostringstream text;
      text << "The time steps of the cells range from " << evaluation.minTimestep << " (the cell "
           << evaluation.minTimestepCell << " of the file) to " << evaluation.maxTimestep
           << "; up to the end time, the cells with the smallest one take " << std::ceil(steps)
           << " steps.";
      logInfo() << text.str().c_str();
    }
  }
}

std::string depthName(CheckDepth depth) {
  switch (depth) {
  case CheckDepth::Files:
    return "files";
  case CheckDepth::Mesh:
    return "mesh";
  case CheckDepth::Model:
    return "model";
  }
  return "unknown";
}

} // namespace

std::optional<CheckOptions> extractCheckOptions(int& argc, char** argv) {
  constexpr std::string_view Option = "--check";
  constexpr std::string_view DepthOption = "--check-depth=";
  std::optional<CheckOptions> options;
  int kept = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string_view argument(argv[i]);
    if (argument == Option) {
      options = options.value_or(CheckOptions{});
    } else if (argument.substr(0, Option.size() + 1) == std::string(Option) + "=") {
      options = options.value_or(CheckOptions{});
      options->reportPath = std::string(argument.substr(Option.size() + 1));
    } else if (argument.starts_with(DepthOption)) {
      options = options.value_or(CheckOptions{});
      const auto depth = argument.substr(DepthOption.size());
      if (depth == "files") {
        options->depth = CheckDepth::Files;
      } else if (depth == "mesh") {
        options->depth = CheckDepth::Mesh;
      } else if (depth == "model") {
        options->depth = CheckDepth::Model;
      } else {
        logError() << "--check-depth is files, mesh or model, not" << std::string(depth);
      }
    } else {
      argv[kept] = argv[i];
      ++kept;
    }
  }
  argc = kept;
  argv[argc] = nullptr;
  return options;
}

int preflight(const std::string& parameterFile,
              const std::vector<std::string>& overrides,
              const std::optional<std::string>& checkpoint,
              const CheckOptions& options) {
  Report report(options.reportPath);
  report["depth"] = depthName(options.depth);
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
  // an error stops the check as it would stop a run; the report gives the step and the error
  const AbortHandlerScope stopOnError([&report](const utils::LogRecord& record) {
    report.stop(record);
    std::_Exit(1);
  });

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

  // the easi models: the files they refer to, whether they parse, and what they supply
  report.step("easi");
  auto& easiModels = report["easi"] = json::Value::array();
  std::unique_ptr<easi::Component> material;
  for (const auto& entry : paths) {
    std::error_code error;
    const auto name = joined(entry.field);
    if (entry.path.empty() || !isEasiModel(name) || !runOpens(entry.field, seissolParameters) ||
        !std::filesystem::is_regular_file(entry.path, error)) {
      continue;
    }
    easiModels.push(checkEasiModel(name,
                                   entry.path,
                                   seissolParameters,
                                   report,
                                   name == "equations.materialfilename" ? &material : nullptr));
  }

  report.step("mesh");
  report["mesh"] = meshSection(seissolParameters.mesh, report);

  if (checkpoint.has_value()) {
    report.step("checkpoint");
    std::error_code error;
    if (!std::filesystem::is_regular_file(checkpoint.value(), error)) {
      report.problem("The checkpoint " + checkpoint.value() + " does not exist.");
    }
    // (see setupCheckpointing)
    if (seissolParameters.model.configs().size() > 1) {
      report.problem("A run whose cells compute in several configurations cannot restart from a "
                     "checkpoint yet.");
    }
  }

  report.step("output");
  report["output"] = outputSection(seissolParameters, report);

  if (options.depth != CheckDepth::Files) {
    checkMeshContent(seissolParameters, options.depth, material.get(), report);
  }

  return report.finish();
}

} // namespace seissol::initializer
