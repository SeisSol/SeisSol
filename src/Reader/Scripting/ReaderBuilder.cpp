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
#include "Reader/Scripting/EasiWalker.h"
#include "Reader/Scripting/LuaFrontend.h"
#include "Reader/Scripting/LuaReader.h"

#include <cstdlib>
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

/// A Lua model: translated into a program by the Lua frontend (LuaFrontend.h), with the
/// interpreted reader as the oracle and the fallback of the compiled one; where the frontend does
/// not translate it, evaluated by Lua, with a warning that says why.
std::unique_ptr<DataReader> buildLua(const std::string& path) {
  const std::string code = readFile(path);
  const lua::LuaOrigin origin{path, 1, ':'};
  try {
    auto program = lua::translateLuaModule(code, origin);
    logInfo() << "The Lua model" << path << "is translated into a program (the Lua frontend).";
#ifdef USE_LUA
    return std::make_unique<CompiledReader>(std::move(program), std::make_unique<LuaReader>(code));
#else
    return std::make_unique<CompiledReader>(std::move(program), nullptr);
#endif
  } catch (const lua::LuaError& error) {
    const std::string where = error.line() > 0 ? origin.at(error.line()) + ": " : std::string();
#ifdef USE_LUA
    logWarning() << "The Lua model" << path << "is not translated, and is evaluated by Lua instead:"
                 << where + error.message();
    return std::make_unique<LuaReader>(code);
#else
    logError() << "The Lua model" << path
               << "is not translated, and SeisSol is built without Lua:" << where + error.message();
    return nullptr;
#endif
  }
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
  const std::string source = readFile(path);
  expr::SderivOptions options;
  options.file = path;
  auto program = expr::compileSderivModule(source, options);
  logInfo() << "Compiled the sderiv module" << path << "with" << program.outputs().size()
            << "outputs.";
  return std::make_unique<CompiledReader>(std::move(program), nullptr);
}

std::unique_ptr<DataReader> makeEasiReader(const std::string& path,
                                           const std::vector<std::string>& dimensions) {
  return std::make_unique<EasiReader>(path, dimensions);
}

/// An easi file: translated into a program by the easi walker where it translates the file, else
/// evaluated by easi, with a note that names the component it does not translate.
/// SEISSOL_EASI=easi evaluates every file by easi; SEISSOL_EASI=check compares the program with
/// easi at its first call, and falls back to easi where they differ.
std::unique_ptr<DataReader> buildEasi(const std::string& path,
                                      const std::vector<std::string>& dimensions) {
  const char* given = std::getenv("SEISSOL_EASI");
  const std::string mode = given == nullptr ? std::string("walker") : std::string(given);
  if (mode == "easi") {
    return makeEasiReader(path, dimensions);
  }
  if (mode != "walker" && mode != "check") {
    logError() << "SEISSOL_EASI is easi, walker or check, not" << mode;
  }
  try {
    // the whole model once, so that what the walker does not translate shows now, not at a call
    EasiQuery query;
    query.dimensions = dimensions;
    query.parameters = easiParameters(path, dimensions);
    static_cast<void>(easiProgram(path, query));
    auto reader = std::make_unique<EasiProgramReader>(
        path, dimensions, mode == "check" ? &makeEasiReader : nullptr);
    logInfo() << "The easi file" << path << "is translated into a program (the easi walker).";
    return reader;
  } catch (const EasiUnsupported& error) {
    logInfo() << "easi evaluates the file" << path << ":" << std::string(error.what());
  } catch (const EasiModelError& error) {
    logWarning() << "The easi walker does not read the file" << path << ":"
                 << std::string(error.what()) << "; easi evaluates it.";
  }
  return makeEasiReader(path, dimensions);
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
#ifdef USE_LUA
    return std::make_unique<LuaReader>(readFile(fileOf(path)));
#else
    return nullptr;
#endif
  }
  if (kind == "easi") {
    return std::make_unique<EasiReader>(fileOf(path), defaultInArgs);
  }
  return nullptr;
}

std::optional<expr::Program> buildProgram(const std::string& path,
                                          std::string* reason,
                                          const expr::SderivOptions& sderivOptions) {
  const auto kind = kindOf(path);
  if (kind == "sderiv") {
    const std::string source = readFile(fileOf(path));
    expr::SderivOptions options = sderivOptions;
    options.file = fileOf(path);
    try {
      return expr::compileSderivModule(source, options);
    } catch (const expr::SderivError& error) {
      // an error of the core format names its line and column itself
      if (error.position() >= 0) {
        logError() << "The sderiv module" << path << "does not compile at byte" << error.position()
                   << ":" << error.what();
      } else {
        logError() << "The sderiv module" << path << "does not compile:" << error.what();
      }
    }
  }
  if (kind == "lua") {
    const std::string file = fileOf(path);
    const lua::LuaOrigin origin{file, 1, ':'};
    try {
      return lua::translateLuaModule(readFile(file), origin);
    } catch (const lua::LuaError& error) {
      if (reason != nullptr) {
        *reason =
            (error.line() > 0 ? origin.at(error.line()) + ": " : std::string()) + error.message();
      }
      return std::nullopt;
    }
  }
  if (reason != nullptr) {
    *reason = "an easi file does not compile to a program";
  }
  return std::nullopt;
}

std::unique_ptr<DataReader> buildReader(const std::string& path,
                                        const std::vector<std::string>& defaultInArgs) {

  logInfo() << "Reading script" << path << "with default args" << defaultInArgs;

  // by its prefix, else by its extension, as buildProgram() reads it
  const auto kind = kindOf(path);
  if (kind == "lua") {
    return buildLua(fileOf(path));
  }
  if (kind == "sderiv") {
    return buildSderiv(fileOf(path));
  }
  return buildEasi(fileOf(path), defaultInArgs);
}

} // namespace seissol::reader::scripting
