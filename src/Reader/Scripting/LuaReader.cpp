// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Scripting/LuaReader.h"

#include "Parallel/OpenMP.h"
#include "Reader/Scripting/DataTable.h"

#include <cassert>
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

namespace seissol::reader::scripting {

namespace {

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

  const auto modCode = loadCodeLmathX() + code;

  const int status = luaL_dostring(luaState, modCode.data());

  if (status != LUA_OK) {
    logError() << "Loading a script failed:" << modCode;
  }

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

  // Build an empty fields table; bindField populates it.
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

/// Read M.field_specs (table at top of stack) into a vector<FieldSpec>.
std::vector<FieldSpec> readFieldSpecs(lua_State* luaState) {
  std::vector<FieldSpec> out;
  if (!lua_istable(luaState, -1)) {
    lua_pop(luaState, 1);
    return out;
  }
  const lua_Integer n = luaL_len(luaState, -1);
  out.reserve(static_cast<std::size_t>(n));
  for (lua_Integer i = 1; i <= n; ++i) {
    lua_geti(luaState, -1, i); // field_specs[i]
    if (!lua_istable(luaState, -1)) {
      lua_pop(luaState, 1);
      continue;
    }
    FieldSpec fs;
    lua_getfield(luaState, -1, "name");
    if (lua_type(luaState, -1) == LUA_TSTRING) {
      fs.name = lua_tostring(luaState, -1);
    }
    lua_pop(luaState, 1);
    lua_getfield(luaState, -1, "kind");
    if (lua_type(luaState, -1) == LUA_TSTRING) {
      fs.kind = lua_tostring(luaState, -1);
    }
    lua_pop(luaState, 1);
    lua_getfield(luaState, -1, "file");
    if (lua_type(luaState, -1) == LUA_TSTRING) {
      fs.file = lua_tostring(luaState, -1);
    }
    lua_pop(luaState, 1);
    lua_getfield(luaState, -1, "interpolation");
    if (lua_type(luaState, -1) == LUA_TSTRING) {
      fs.interpolation = lua_tostring(luaState, -1);
    }
    lua_pop(luaState, 1);
    lua_getfield(luaState, -1, "parameters");
    if (lua_istable(luaState, -1)) {
      fs.parameters = readStringArray(luaState); // pops the table
    } else {
      lua_pop(luaState, 1);
    }
    out.push_back(std::move(fs));
    lua_pop(luaState, 1); // pop field_specs[i]
  }
  lua_pop(luaState, 1); // pop the field_specs table
  return out;
}

extern "C" {

// launcher for field samplers (e.g. ASAGI); is called from Lua code
int fieldSampleTrampoline(lua_State* luaState) {
  // Stack: [self_table, coord_1, coord_2, ...]
  auto* closure =
      static_cast<LuaReader::FieldClosure*>(lua_touserdata(luaState, lua_upvalueindex(1)));
  if (closure != nullptr || !closure->sampler) {
    return luaL_error(luaState, "field sampler not bound");
  }
  const int nArgs = lua_gettop(luaState) - 1; // everything after self
  if (nArgs > 16) {
    return luaL_error(luaState, "field sampler: too many coords (%d > 16)", nArgs);
  }
  double coords[16];
  for (int i = 0; i < nArgs; ++i) {
    coords[i] = lua_tonumber(luaState, 2 + i);
  }
  if (closure->nOutputs > 16) {
    return luaL_error(luaState, "field sampler: too many outputs (%zu > 16)", closure->nOutputs);
  }
  double results[16] = {0};
  try {
    closure->sampler(coords, results);
  } catch (const std::exception& e) {
    return luaL_error(luaState, "field sampler threw: %s", e.what());
  } catch (...) {
    return luaL_error(luaState, "field sampler threw: unknown exception");
  }
  for (std::size_t i = 0; i < closure->nOutputs; ++i) {
    lua_pushnumber(luaState, results[i]);
  }
  return static_cast<int>(closure->nOutputs);
}

} // extern "C"

} // namespace

LuaReader::~LuaReader() = default;

LuaReader::LuaReader(const std::string& code) : code_(code) {

  const auto luaState = loadScript(code_);

  readMetadata(luaState);

  lua_close(luaState.luaState);
}

void LuaReader::readMetadata(const LuaStateState& state) {
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
  inputs_ = readStringArray(state.luaState);

  lua_getfield(state.luaState, -1, "output_parameters");
  outputs_ = readStringArray(state.luaState);

  lua_getfield(state.luaState, -1, "field_specs");
  fieldSpecs_ = readFieldSpecs(state.luaState);

  lua_pop(state.luaState, 1); // pop M
}

void LuaReader::bindField(const std::string& name, FieldSampler sampler) {
  // Find the matching field spec to know how many outputs to expect.
  const FieldSpec& spec = [&]() {
    for (const auto& s : fieldSpecs_) {
      if (s.name == name) {
        return s;
      }
    }
    // nothing found; fail
    logError() << "bindField: no field '" + name + "' declared in module";
    throw;
  }();
  auto closure = std::make_unique<FieldClosure>();
  closure->owner = this;
  closure->sampler = std::move(sampler);
  closure->nOutputs = spec.parameters.size();

  fieldClosures_[name] = std::move(closure);
}

void LuaReader::registerFields(const LuaStateState& state) {
  for (const auto& [name, cls] : fieldClosures_) {
    // Build a Lua field object: a table with a `:sample(...)` method, where
    // `sample` is a C-closure capturing the closure pointer as its upvalue.
    lua_rawgeti(state.luaState, LUA_REGISTRYINDEX, state.refFields); // [fields]
    lua_newtable(state.luaState);                                    // [fields, fobj]
    lua_pushlightuserdata(state.luaState, cls.get());                // [fields, fobj, ud]
    lua_pushcclosure(state.luaState, &fieldSampleTrampoline, 1);     // [fields, fobj, csample]
    lua_setfield(state.luaState, -2, "sample");                      // [fields, fobj]
    lua_setfield(state.luaState, -2, name.c_str());                  // [fields]
    lua_pop(state.luaState, 1);                                      // []
  }
}

void LuaReader::call(const scripting::DataTable& table) {

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

  const auto numStates = OpenMP::threadCount();
  if (luaStates_.size() < numStates) {
    luaStates_.resize(numStates);
  }

#pragma omp parallel
  {
    if (luaStates_[OpenMP::threadId()].luaState == nullptr) {
      luaStates_[OpenMP::threadId()] = loadScript(code_);
      registerFields(luaStates_[OpenMP::threadId()]);
    }

    const auto& luaStateState = luaStates_[OpenMP::threadId()];
    const auto& luaState = luaStateState.luaState;

    // Save stack size
    const auto top = lua_gettop(luaState);

#pragma omp for schedule(static)
    for (std::size_t point = 0; point < table.numPoints(); ++point) {
      // Push function and arguments to stack
      lua_getfield(luaState, -1, "evaluate");
      lua_rawgeti(luaState, LUA_REGISTRYINDEX, luaStateState.refEvaluate);
      lua_rawgeti(luaState, LUA_REGISTRYINDEX, luaStateState.refFields);

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

      // +1 due to the `fields` field
      if (lua_pcall(luaState, inputs_.size() + 1, outputs_.size(), 0) != LUA_OK) {
        logError() << "Error running Lua function:" << lua_tostring(luaState, -1);
      }

      const std::size_t nout = outputs_.size();
      for (std::size_t i = 0; i < outVarMap.size(); ++i) {
        const auto& outIdx = outVarMap[i];

        const int argPos = static_cast<int>(i) - static_cast<int>(nout);

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

      // Reset stack size to value before function call
      // This should avoid stack overflows
      lua_settop(luaState, top);
    }
  }
}

} // namespace seissol::reader::scripting
