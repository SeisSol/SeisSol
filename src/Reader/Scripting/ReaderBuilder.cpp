// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ReaderBuilder.h"

#include "Expr/SderivFrontend.h"
#include "Reader/Scripting/CompiledReader.h"
#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/EasiReader.h"
#include "Reader/Scripting/LuaReader.h"
#include "Reader/Scripting/LuaTracer.h"

#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::reader::scripting {

namespace {

/// Everything after the first colon. Not parts[1]: a path may contain further
/// colons, and splitting on all of them then taking one piece is how
/// "lua:/net:2/model.lua" quietly becomes "/net".
std::string stripPrefix(const std::string& path) {
  const auto colon = path.find(':');
  return colon == std::string::npos ? path : path.substr(colon + 1);
}

std::string readFile(const std::string& path) {
  std::ifstream file(path);
  if (!file) {
    // an unchecked stream would give an empty script, and a reader that
    // silently returns nothing
    logError() << "Could not open the script" << path << ".";
  }
  std::stringstream code;
  code << file.rdbuf();
  return code.str();
}

/// Try to trace, and say clearly what happened either way. The interpreted
/// reader is handed in as both oracle and fallback, so a program that traces
/// but disagrees with it at run time still lands on the path that works.
std::unique_ptr<DataReader> buildLua(const std::string& path) {
  const std::string code = readFile(stripPrefix(path));

  TraceFailure failure;
  const TraceOptions options;
  auto program = traceLuaModule(code, options, failure);
  if (!program.has_value()) {
    // TraceFailure::reason is documented as already formatted for a log line
    // and as carrying the position, so it is not decorated further here.
    logWarning() << "The Lua model" << path
                 << "could not be traced, and is evaluated through the interpreter instead:"
                 << failure.reason;
    return std::make_unique<LuaReader>(code);
  }

  logInfo() << "Traced the Lua model" << path << "into a compiled program.";
  return std::make_unique<CompiledReader>(std::move(*program), std::make_unique<LuaReader>(code));
}

/// A .sderiv module names its own outputs, so nothing is passed in and nothing
/// can drift out of step with the file.
///
/// No reference reader is handed to CompiledReader, and that is not an
/// omission: there is no interpreter for this language, so there is neither an
/// oracle to compare against nor a path to fall back TO. A failure is a
/// configuration error, which is what makes the frontend's parse diagnostics
/// load-bearing rather than a convenience.
std::unique_ptr<DataReader> buildSderiv(const std::string& path) {
  const std::string source = readFile(stripPrefix(path));
  auto program = expr::compileSderivModule(source);
  logInfo() << "Compiled the sderiv module" << path << "with" << program.outputs().size()
            << "outputs.";
  return std::make_unique<CompiledReader>(std::move(program), nullptr);
}

bool endsWith(const std::string& text, const std::string& suffix) {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/// The kind of a script path: its prefix, else its extension; easi without either.
std::string kindOf(const std::string& path) {
  const auto parts = utils::StringUtils::split(path, ':');
  if (parts.size() > 1 && (parts[0] == "easi" || parts[0] == "lua" || parts[0] == "sderiv")) {
    return parts[0];
  }
  if (endsWith(path, ".lua")) {
    return "lua";
  }
  if (endsWith(path, ".sderiv")) {
    return "sderiv";
  }
  return "easi";
}

/// The file of a script path, without its prefix.
std::string fileOf(const std::string& path) {
  const auto kind = kindOf(path);
  return path.compare(0, kind.size() + 1, kind + ":") == 0 ? stripPrefix(path) : path;
}

} // namespace

std::unique_ptr<DataReader> buildInterpretedReader(const std::string& path,
                                                   const std::vector<std::string>& defaultInArgs) {
  const auto kind = kindOf(path);
  if (kind == "lua") {
    return std::make_unique<LuaReader>(readFile(fileOf(path)));
  }
  if (kind == "easi") {
    return std::make_unique<EasiReader>(fileOf(path), defaultInArgs);
  }
  return nullptr;
}

std::optional<expr::Program> buildProgram(const std::string& path, std::string* reason) {
  const auto kind = kindOf(path);
  if (kind == "sderiv") {
    const std::string source = readFile(fileOf(path));
    try {
      return expr::compileSderivModule(source);
    } catch (const expr::SderivError& error) {
      logError() << "The sderiv module" << path << "does not compile at byte" << error.position()
                 << ":" << error.what();
    }
  }
  if (kind == "lua") {
    TraceFailure failure;
    auto program = traceLuaModule(readFile(fileOf(path)), {}, failure);
    if (!program.has_value() && reason != nullptr) {
      *reason = failure.reason;
    }
    return program;
  }
  if (reason != nullptr) {
    *reason = "an easi file does not compile to a program";
  }
  return std::nullopt;
}

std::unique_ptr<DataReader> buildReader(const std::string& path,
                                        const std::vector<std::string>& defaultInArgs) {

  logInfo() << "Reading script" << path << "with default args" << defaultInArgs;

  const auto parts = utils::StringUtils::split(path, ':');
  if (parts.size() == 1) {
    return std::make_unique<EasiReader>(path, defaultInArgs);
  }
  if (parts[0] == "easi") {
    // The easi walker is not built yet, so this stays interpreted. When it
    // lands it takes the same shape as buildLua: walk, and on refusal fall back
    // to EasiReader with a warning naming why.
    return std::make_unique<EasiReader>(stripPrefix(path), defaultInArgs);
  }
  if (parts[0] == "lua") {
    return buildLua(path);
  }
  if (parts[0] == "sderiv") {
    return buildSderiv(path);
  }

  logError() << "The script" << path << "does not have a built-in reader.";

  return nullptr;
}

} // namespace seissol::reader::scripting
