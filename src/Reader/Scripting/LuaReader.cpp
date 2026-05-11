// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Scripting/LuaReader.h"

#include "Parallel/OpenMP.h"

#include <cassert>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <utils/logger.h>

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include <lauxlib.h>
#include <lua.h>
#include <lualib.h>
}

#ifdef SEISSOL_LUA_LMATHX
extern "C" int luaopen_mathx(lua_State* luaState);

// Lua changes a lot between the different versions, hence we need quite some ifdefs here
namespace {

void loadLmathx(lua_State* luaState) {
#if LUA_VERSION_NUM < 502
  lua_pushcfunction(luaState, luaopen_mathx);
  lua_pushliteral(luaState, "mathx");
  lua_call(luaState, 1, 0);
#else
  luaL_requiref(luaState, "mathx", luaopen_mathx, 1);
  lua_pop(luaState, 1);
#endif
}

std::string loadCodeLmathX() {
  // cf. example code for mathx: for Lua 5.3 and upwards, we need to supplement the math table.
  // (for 5.2 and lower, that's done automatically)
#if LUA_VERSION_NUM < 503
  return "";
#else
  return "for fname,fval in pairs(mathx) do math[fname] = fval end;\n";
#endif
  return "";
}

} // namespace
#else
namespace {
void loadLmathx(lua_State* luaState) {}

std::string loadCodeLmathX() { return ""; }
} // namespace
#endif

namespace seissol::reader::scripting {

namespace {
lua_State* loadScript(const std::string& code) {
  auto* luaState = luaL_newstate();
  luaL_openlibs(luaState);

  loadLmathx(luaState);

  const auto modCode = loadCodeLmathX() + code;

  const int status = luaL_dostring(luaState, modCode.data());

  if (status != LUA_OK) {
    logError() << "Loading a script failed:" << modCode;
  }

  return luaState;
}
} // namespace

LuaReader::~LuaReader() = default;

LuaReader::LuaReader(const std::string& code) : code_(code) {

  auto* luaState = loadScript(code_);
  // TODO: parse input/output params
}

void LuaReader::call(const scripting::DataTable& table) {

  const auto& entries = table.dataEntries();

  std::vector<std::size_t> inVarMap(inVars_.size());
  std::vector<std::size_t> outVarMap(outVars_.size());

  for (std::size_t i = 0; i < entries.size(); ++i) {
    for (std::size_t j = 0; j < inVars_.size(); ++j) {
      if (inVars_[j] == entries[i].name) {
        inVarMap[j] = i;
      }
    }
    for (std::size_t j = 0; j < outVars_.size(); ++j) {
      if (outVars_[j] == entries[i].name) {
        outVarMap[j] = i;
      }
    }
  }

  const auto numStates = OpenMP::threadCount();
  if (luaStates_.size() < numStates) {
    luaStates_.resize(numStates);
  }

#pragma omp parallel
  {
    auto& luaState = luaStates_[OpenMP::threadId()];

    if (luaState == nullptr) {
      luaState = loadScript(code_);
    }

    // Save stack size
    const auto top = lua_gettop(luaState);

#pragma omp for schedule(static)
    for (std::size_t point = 0; point < table.numPoints(); ++point) {
      // Push function and arguments to stack
      lua_getfield(luaState, -1, "evaluate");

      for (const auto& inIdx : inVarMap) {
        switch (entries[inIdx].datatype) {
        case DataType::F32: {
          lua_pushnumber(luaState, entries[inIdx].getValue<float>(point));
          break;
        }
        case DataType::F64: {
          lua_pushnumber(luaState, entries[inIdx].getValue<double>(point));
          break;
        }
        case DataType::I32: {
          lua_pushinteger(luaState, entries[inIdx].getValue<std::int32_t>(point));
          break;
        }
        case DataType::I64: {
          lua_pushinteger(luaState, entries[inIdx].getValue<std::int64_t>(point));
          break;
        }
        default: {
          logError() << "Error with Lua: unknown input datatype.";
        }
        }
      }

      if (lua_pcall(luaState, inVars_.size(), outVars_.size(), 0) != 0) {
        logError() << "Error running Lua function:" << lua_tostring(luaState, -1);
      }

      for (const auto& outIdx : outVarMap) {
        int res = 0;
        switch (entries[outIdx].datatype) {
        case DataType::F32: {
          const auto result = lua_tonumberx(luaState, -1, &res);
          if (res == 0) {
            logError() << "Incorrect Lua data return type.";
          }
          entries[outIdx].setValue<float>(point, result);
          break;
        }
        case DataType::F64: {
          const auto result = lua_tonumberx(luaState, -1, &res);
          if (res == 0) {
            logError() << "Incorrect Lua data return type.";
          }
          entries[outIdx].setValue<double>(point, result);
          break;
        }
        case DataType::I32: {
          const auto result = lua_tointegerx(luaState, -1, &res);
          if (res == 0) {
            logError() << "Incorrect Lua data return type.";
          }
          entries[outIdx].setValue<std::int32_t>(point, result);
          break;
        }
        case DataType::I64: {
          const auto result = lua_tointegerx(luaState, -1, &res);
          if (res == 0) {
            logError() << "Incorrect Lua data return type.";
          }
          entries[outIdx].setValue<std::int64_t>(point, result);
          break;
        }
        default: {
          logError() << "Error with Lua: unknown input datatype.";
        }
        }
        lua_pop(luaState, -1);
      }

      // Reset stack size to value before function call
      // This should avoid stack overflows
      lua_settop(luaState, top);
    }
  }
}

} // namespace seissol::reader::scripting
