// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LuaOracle.h"

#ifdef USE_LUA

#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Parallel/OpenMP.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

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

namespace seissol::unit_test {

using namespace seissol::reader::scripting;

namespace {

// The vocabulary of the earlier tracer, `ssol`, computing: lt, le, gt, ge, eq and the logic on
// booleans and numbers (a number other than zero is true), select, and the builtins of sderiv.
bool truthy(lua_State* luaState, int index) {
  if (lua_isboolean(luaState, index)) {
    return lua_toboolean(luaState, index) != 0;
  }
  return lua_tonumber(luaState, index) != 0.0;
}
int ssolLt(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(lua_tonumber(l, 1) < lua_tonumber(l, 2)));
  return 1;
}
int ssolLe(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(lua_tonumber(l, 1) <= lua_tonumber(l, 2)));
  return 1;
}
int ssolGt(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(lua_tonumber(l, 1) > lua_tonumber(l, 2)));
  return 1;
}
int ssolGe(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(lua_tonumber(l, 1) >= lua_tonumber(l, 2)));
  return 1;
}
int ssolEq(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(lua_tonumber(l, 1) == lua_tonumber(l, 2)));
  return 1;
}
int ssolAnd(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(truthy(l, 1) && truthy(l, 2)));
  return 1;
}
int ssolOr(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(truthy(l, 1) || truthy(l, 2)));
  return 1;
}
int ssolNot(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(!truthy(l, 1)));
  return 1;
}
int ssolSelect(lua_State* l) {
  lua_pushvalue(l, truthy(l, 1) ? 2 : 3);
  return 1;
}
template <expr::Fn F>
int ssolArithmetic(lua_State* l) {
  std::array<double, 2> args{};
  for (int i = 0; i < expr::arity(F); ++i) {
    args[static_cast<std::size_t>(i)] = luaL_checknumber(l, i + 1);
  }
  lua_pushnumber(l, expr::applyPw<double>(F, args.data()));
  return 1;
}

void registerSsol(lua_State* luaState) {
  lua_newtable(luaState);
  const std::array<std::pair<const char*, lua_CFunction>, 13> functions{
      {{"lt", ssolLt},
       {"le", ssolLe},
       {"gt", ssolGt},
       {"ge", ssolGe},
       {"eq", ssolEq},
       {"land", ssolAnd},
       {"lor", ssolOr},
       {"lnot", ssolNot},
       {"select", ssolSelect},
       {"sign", ssolArithmetic<expr::Fn::Sign>},
       {"mod", ssolArithmetic<expr::Fn::Mod>},
       {"atan2", ssolArithmetic<expr::Fn::Atan2>},
       {"pow", ssolArithmetic<expr::Fn::Pow>}}};
  for (const auto& [name, function] : functions) {
    lua_pushcfunction(luaState, function);
    lua_setfield(luaState, -2, name);
  }
  lua_pushnumber(luaState, M_PI);
  lua_setfield(luaState, -2, "pi");
  lua_pushnumber(luaState, 9.80665);
  lua_setfield(luaState, -2, "g");
  lua_setglobal(luaState, "ssol");
}

[[noreturn]] void throwLuaError(lua_State* luaState, const std::string& context) {
  std::string msg = context + ": ";
  if (lua_type(luaState, -1) == LUA_TSTRING) {
    msg += lua_tostring(luaState, -1);
  } else {
    msg += "(no message)";
  }
  lua_pop(luaState, 1);
  logError() << "Lua error:" << msg.c_str();
  throw;
}

LuaStateState loadScript(const std::string& code) {
  auto* luaState = luaL_newstate();

  if (luaState == nullptr) {
    throwLuaError(luaState, "creating state");
  }

  luaL_openlibs(luaState);

  loadLmathx(luaState);

  registerSsol(luaState);

  const auto modCode = loadCodeLmathX() + code;

  LuaStateState state{};

  if (luaL_loadbufferx(luaState, modCode.c_str(), modCode.size(), "script", "t") != LUA_OK) {
    throwLuaError(luaState, "loading module");
  }
  if (lua_pcall(luaState, 0, 1, 0) != LUA_OK) {
    throwLuaError(luaState, "executing module");
  }
  // Stack: [M]
  if (!lua_istable(luaState, -1)) {
    lua_pop(luaState, 1);
    logError() << "Lua: module did not return a table";
  }
  // Save M in the registry
  lua_pushvalue(luaState, -1);                             // [M, M]
  state.refModule = luaL_ref(luaState, LUA_REGISTRYINDEX); // [M]
  // Save M.evaluate in the registry too (fast path)
  lua_getfield(luaState, -1, "evaluate"); // [M, evaluate]
  if (!lua_isfunction(luaState, -1)) {
    lua_pop(luaState, 2);
    logError() << "Lua: module has no `evaluate` function";
  }
  state.refEvaluate = luaL_ref(luaState, LUA_REGISTRYINDEX); // [M]

  // An empty fields table: the oracle has no grids.
  lua_newtable(luaState);                                  // [M, fields]
  state.refFields = luaL_ref(luaState, LUA_REGISTRYINDEX); // [M]

  lua_pop(luaState, 1); // []

  state.luaState = luaState;

  return state;
}

// ---------------------------------------------------------------------------
// Helpers for inspecting the Lua module table
// ---------------------------------------------------------------------------

/// Pop the table at top of stack and turn it into vector<string>.
/// Expects 1-indexed string entries (Lua array convention).
std::vector<std::string> readStringArray(lua_State* luaState) {
  std::vector<std::string> out;
  if (!lua_istable(luaState, -1)) {
    lua_pop(luaState, 1);
    return out;
  }
  const lua_Integer n = luaL_len(luaState, -1);
  out.reserve(static_cast<std::size_t>(n));
  for (lua_Integer i = 1; i <= n; ++i) {
    lua_geti(luaState, -1, i);
    if (lua_type(luaState, -1) == LUA_TSTRING) {
      out.emplace_back(lua_tostring(luaState, -1));
    }
    lua_pop(luaState, 1);
  }
  lua_pop(luaState, 1);
  return out;
}

} // namespace

LuaOracle::~LuaOracle() {
  for (auto& state : luaStates_) {
    lua_close(state.luaState);
  }
}

LuaOracle::LuaOracle(const std::string& code) : code_(code) {

  const auto luaState = loadScript(code_);

  readMetadata(luaState);

  lua_close(luaState.luaState);
}

namespace {

/// The parameter names of the function at the top of the stack, which it leaves in place.
std::vector<std::string> parameterNames(lua_State* luaState) {
  lua_Debug debug{};
  lua_pushvalue(luaState, -1);
  lua_getinfo(luaState, ">u", &debug);
  std::vector<std::string> names;
  for (int i = 1; i <= debug.nparams; ++i) {
    const char* name = lua_getlocal(luaState, nullptr, i);
    if (name == nullptr) {
      break;
    }
    names.emplace_back(name);
  }
  return names;
}

/// The string keys of the table at `index`, sorted.
std::vector<std::string> sortedKeys(lua_State* luaState, int index) {
  std::vector<std::string> keys;
  lua_pushnil(luaState);
  while (lua_next(luaState, index) != 0) {
    if (lua_type(luaState, -2) == LUA_TSTRING) {
      keys.emplace_back(lua_tostring(luaState, -2));
    }
    lua_pop(luaState, 1);
  }
  std::sort(keys.begin(), keys.end());
  return keys;
}

} // namespace

void LuaOracle::readMetadata(const LuaStateState& state) {
  lua_rawgeti(state.luaState, LUA_REGISTRYINDEX, state.refModule); // [M]

  lua_getfield(state.luaState, -1, "version");
  if (lua_type(state.luaState, -1) == LUA_TSTRING) {
    version_ = lua_tostring(state.luaState, -1);
  }
  lua_pop(state.luaState, 1);

  lua_getfield(state.luaState, -1, "source_file");
  if (lua_type(state.luaState, -1) == LUA_TSTRING) {
    sourceFile_ = lua_tostring(state.luaState, -1);
  }
  lua_pop(state.luaState, 1);

  lua_getfield(state.luaState, -1, "input_parameters");
  const bool declaresInputs = !lua_isnil(state.luaState, -1);
  inputs_ = readStringArray(state.luaState);

  lua_getfield(state.luaState, -1, "output_parameters");
  const bool declaresOutputs = !lua_isnil(state.luaState, -1);
  outputs_ = readStringArray(state.luaState);

  lua_getfield(state.luaState, -1, "state"); // [M, state]
  if (lua_istable(state.luaState, -1)) {
    stateNames_ = sortedKeys(state.luaState, lua_gettop(state.luaState));
    for (const auto& name : stateNames_) {
      lua_getfield(state.luaState, -1, name.c_str());
      stateInitial_.push_back(lua_tonumber(state.luaState, -1));
      lua_pop(state.luaState, 1);
    }
  }
  lua_pop(state.luaState, 1); // [M]

  // The arguments of evaluate: the declared inputs in their order, as scripts have always had
  // them; or, as the tracer reads them, the parameter names after `fields`, of which those naming
  // a state receive its value of the previous call.
  lua_getfield(state.luaState, -1, "evaluate"); // [M, evaluate]
  const auto parameters = parameterNames(state.luaState);
  if (declaresInputs) {
    if (!stateNames_.empty()) {
      logError() << "Lua: a module with state names its inputs by the parameters of `evaluate`; "
                    "remove `input_parameters`.";
    }
    for (std::size_t i = 0; i < inputs_.size(); ++i) {
      arguments_.push_back(Argument{false, i});
    }
  } else {
    for (std::size_t i = 1; i < parameters.size(); ++i) {
      const auto found = std::find(stateNames_.begin(), stateNames_.end(), parameters[i]);
      if (found != stateNames_.end()) {
        arguments_.push_back(Argument{true, static_cast<std::size_t>(found - stateNames_.begin())});
      } else {
        arguments_.push_back(Argument{false, inputs_.size()});
        inputs_.push_back(parameters[i]);
      }
    }
  }

  // Outputs named by a returned table, and not declared, are found by one call at zero.
  if (!declaresOutputs) {
    lua_rawgeti(state.luaState, LUA_REGISTRYINDEX, state.refFields);
    for (const auto& argument : arguments_) {
      lua_pushnumber(state.luaState, argument.state ? stateInitial_[argument.index] : 0.0);
    }
    if (lua_pcall(state.luaState, static_cast<int>(arguments_.size()) + 1, 1, 0) != LUA_OK ||
        !lua_istable(state.luaState, -1)) {
      logError() << "Lua: name the outputs of `evaluate` with `output_parameters`, or return them "
                    "in a table.";
    }
    outputs_ = sortedKeys(state.luaState, lua_gettop(state.luaState));
  }
  lua_pop(state.luaState, 1); // [M]

  lua_pop(state.luaState, 1); // pop M
}

void LuaOracle::call(const reader::scripting::DataTable& table) {

  const auto& entries = table.dataEntries();

  std::vector<std::size_t> inVarMap(inputs_.size());
  std::vector<std::size_t> outVarMap(outputs_.size());

  for (std::size_t i = 0; i < entries.size(); ++i) {
    for (std::size_t j = 0; j < inputs_.size(); ++j) {
      if (inputs_[j] == entries[i].name) {
        inVarMap[j] = i;
      }
    }
    for (std::size_t j = 0; j < outputs_.size(); ++j) {
      if (outputs_[j] == entries[i].name) {
        outVarMap[j] = i;
      }
    }
  }

  if (!stateNames_.empty() && statePoints_ != table.numPoints()) {
    stateValues_.clear();
    for (const double initial : stateInitial_) {
      stateValues_.emplace_back(table.numPoints(), initial);
    }
    statePoints_ = table.numPoints();
  }

  // With positional results, a state takes the output of its name.
  std::vector<std::size_t> stateOutput(stateNames_.size(), outputs_.size());
  for (std::size_t s = 0; s < stateNames_.size(); ++s) {
    const auto found = std::find(outputs_.begin(), outputs_.end(), stateNames_[s]);
    stateOutput[s] = static_cast<std::size_t>(found - outputs_.begin());
  }

  const auto numStates = OpenMP::threadCount();
  if (luaStates_.size() < numStates) {
    luaStates_.resize(numStates);
  }

#pragma omp parallel
  {
    if (luaStates_[OpenMP::threadId()].luaState == nullptr) {
      luaStates_[OpenMP::threadId()] = loadScript(code_);
    }

    const auto& luaStateState = luaStates_[OpenMP::threadId()];
    const auto& luaState = luaStateState.luaState;

    // Save stack size
    const auto top = lua_gettop(luaState);
    std::vector<double> nextState(stateNames_.size());

#pragma omp for schedule(static)
    for (std::size_t point = 0; point < table.numPoints(); ++point) {
      // Push function and arguments to stack
      lua_rawgeti(luaState, LUA_REGISTRYINDEX, luaStateState.refEvaluate);
      lua_rawgeti(luaState, LUA_REGISTRYINDEX, luaStateState.refFields);

      for (const auto& argument : arguments_) {
        if (argument.state) {
          lua_pushnumber(luaState, stateValues_[argument.index][point]);
          continue;
        }
        const auto inIdx = inVarMap[argument.index];
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

      // +1 due to the `fields` field
      if (lua_pcall(luaState, static_cast<int>(arguments_.size()) + 1, LUA_MULTRET, 0) != LUA_OK) {
        logError() << "Error running Lua function:" << lua_tostring(luaState, -1);
      }

      // The results: positional, or one table that names them. Either way, push the outputs in
      // order, followed by the next values of the states.
      const int results = lua_gettop(luaState) - top;
      if (results == 1 && lua_istable(luaState, top + 1)) {
        const int named = top + 1;
        for (const auto& name : outputs_) {
          lua_getfield(luaState, named, name.c_str());
        }
        for (const auto& name : stateNames_) {
          lua_getfield(luaState, named, name.c_str());
        }
      } else {
        if (results < static_cast<int>(outputs_.size())) {
          logError() << "Lua: `evaluate` returned" << results << "values instead of"
                     << outputs_.size() << ".";
        }
        for (std::size_t i = 0; i < outputs_.size(); ++i) {
          lua_pushvalue(luaState, top + 1 + static_cast<int>(i));
        }
        for (const auto output : stateOutput) {
          if (output == outputs_.size()) {
            logError() << "Lua: a state gets no next value; return an output of its name.";
          }
          lua_pushvalue(luaState, top + 1 + static_cast<int>(output));
        }
      }
      const int first = lua_gettop(luaState) - static_cast<int>(outputs_.size()) -
                        static_cast<int>(stateNames_.size()) + 1;

      for (std::size_t i = 0; i < outVarMap.size(); ++i) {
        const auto& outIdx = outVarMap[i];

        const int argPos = first + static_cast<int>(i);

        int res = 0;
        switch (entries[outIdx].datatype) {
        case DataType::F32: {
          const auto result = lua_tonumberx(luaState, argPos, &res);
          if (res == 0) {
            logError() << "Incorrect Lua data return type.";
          }
          entries[outIdx].setValue<float>(point, result);
          break;
        }
        case DataType::F64: {
          const auto result = lua_tonumberx(luaState, argPos, &res);
          if (res == 0) {
            logError() << "Incorrect Lua data return type.";
          }
          entries[outIdx].setValue<double>(point, result);
          break;
        }
        case DataType::I32: {
          const auto result = lua_tointegerx(luaState, argPos, &res);
          if (res == 0) {
            logError() << "Incorrect Lua data return type.";
          }
          entries[outIdx].setValue<std::int32_t>(point, result);
          break;
        }
        case DataType::I64: {
          const auto result = lua_tointegerx(luaState, argPos, &res);
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
      }

      // All states update at once, after every one of them has been read.
      for (std::size_t s = 0; s < stateNames_.size(); ++s) {
        int res = 0;
        nextState[s] = lua_tonumberx(luaState, first + static_cast<int>(outputs_.size() + s), &res);
        if (res == 0) {
          logError() << "Lua: the state" << stateNames_[s] << "gets no number as its next value.";
        }
      }
      for (std::size_t s = 0; s < stateNames_.size(); ++s) {
        stateValues_[s][point] = nextState[s];
      }

      // Reset stack size to value before function call
      // This should avoid stack overflows
      lua_settop(luaState, top);
    }
  }
}

} // namespace seissol::unit_test

#endif
