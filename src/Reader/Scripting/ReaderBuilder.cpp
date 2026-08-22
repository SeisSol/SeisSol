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
    // FIXED (Package 4): this used to open `path` WITH the "lua:" prefix still
    // on it and never check the stream, so a mistyped path produced an empty
    // script and a reader that silently returned nothing.
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

} // namespace

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
