// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Scripting/LuaTracer.h"

#ifdef USE_LUA

#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include <lauxlib.h>
#include <lua.h>
#include <lualib.h>
}

#if LUA_VERSION_NUM < 504
#error "LuaTracer requires Lua 5.4: it relies on __le having no not-(b<a) fallback"
#endif

namespace seissol::reader::scripting {

namespace {

using expr::Fn;
using expr::NodeId;
using Cause = TraceFailure::Cause;

constexpr const char* SymMetaName = "seissol.expr.sym";
constexpr const char* TracerRegistryKey = "seissol.expr.tracer";

// A traced value. Full userdata rather than a table: a table lets a script add
// fields to it silently, and luaL_testudata gives a cheap type tag.
struct Sym {
  NodeId id{expr::NoNode};
  // Produced by lt/le/eq/and/or/not, hence exactly 0.0 or 1.0. Tracked here and
  // not in the IR because Program.h deliberately has no per-node type; this is
  // a tracer-local invariant that costs the IR nothing.
  bool isBool{false};
};

// ---------------------------------------------------------------------------
// Tracing state, reachable from every C closure through the registry.
// ---------------------------------------------------------------------------
struct TraceState {
  expr::Program* program{nullptr};
  const TraceOptions* options{nullptr};

  // Grid declarations, in M.field_specs order; index doubles as expr::GridId.
  std::vector<FieldSpec> specs;
  std::vector<expr::GridId> gridIds;

  // Conditions built vs. conditions actually fed to select/land/lor/lnot.
  std::vector<NodeId> boolProduced;
  std::vector<NodeId> boolConsumed;

  // Set when a C closure raises through luaL_error, so the pcall handler can
  // report a cause instead of guessing from the message text.
  Cause cause{Cause::UntracedOperator};

  // the script line the raise came from. fail() used to set
  // only `cause`, so every TraceFailure from a net carried line == -1 and the
  // diagnostic pointed at the file but not into it. Captured in fail() because
  // that is the only place the Lua stack is still standing -- after lua_error
  // unwinds there is nothing left to ask.
  int line{-1};
};

TraceState* state(lua_State* luaState) {
  lua_getfield(luaState, LUA_REGISTRYINDEX, TracerRegistryKey);
  auto* s = static_cast<TraceState*>(lua_touserdata(luaState, -1));
  lua_pop(luaState, 1);
  return s;
}

expr::Arena& arena(lua_State* luaState) { return state(luaState)->program->arena(); }

Sym* testSym(lua_State* luaState, int index) {
  return static_cast<Sym*>(luaL_testudata(luaState, index, SymMetaName));
}

int pushSym(lua_State* luaState, NodeId id, bool isBool = false) {
  auto* sym = static_cast<Sym*>(lua_newuserdatauv(luaState, sizeof(Sym), 0));
  sym->id = id;
  sym->isBool = isBool;
  luaL_setmetatable(luaState, SymMetaName);
  return 1;
}

[[noreturn]] void fail(lua_State* luaState, Cause cause, const char* fmt, ...) {
  TraceState* tracer = state(luaState);
  tracer->cause = cause;

  // Level 1 is the Lua frame that called us; a metamethod raise from `z > 0`
  // sits directly under the C closure, which is the line the user wrote.
  lua_Debug debug;
  tracer->line = (lua_getstack(luaState, 1, &debug) != 0 && lua_getinfo(luaState, "l", &debug) != 0)
                     ? debug.currentline
                     : -1;
  va_list args;
  va_start(args, fmt);
  lua_pushvfstring(luaState, fmt, args);
  va_end(args);
  lua_error(luaState);
  std::abort(); // unreachable; lua_error does not return
}

// A usable operand: a traced value, or a plain number lifted to Const. Anything
// else is a hard error, which is where a nil from a mistyped access surfaces.
NodeId operand(lua_State* luaState, int index) {
  if (Sym* sym = testSym(luaState, index)) {
    return sym->id;
  }
  if (lua_type(luaState, index) == LUA_TNUMBER) {
    return arena(luaState).konst(lua_tonumber(luaState, index));
  }
  fail(luaState,
       Cause::UntracedOperator,
       "expression expected, got %s",
       luaL_typename(luaState, index));
}

void budgetCheck(lua_State* luaState) {
  TraceState* s = state(luaState);
  if (s->program->arena().size() > s->options->nodeBudget) {
    fail(luaState,
         Cause::BudgetExceeded,
         "trace exceeded the node budget (%d); an unrolled loop is probably larger than intended",
         static_cast<int>(s->options->nodeBudget));
  }
}

// ---------------------------------------------------------------------------
// Metamethods
// ---------------------------------------------------------------------------
template <Fn F>
int mmUnary(lua_State* luaState) {
  const NodeId a = operand(luaState, 1);
  budgetCheck(luaState);
  return pushSym(luaState, arena(luaState).pw(F, a));
}

template <Fn F>
int mmBinary(lua_State* luaState) {
  const NodeId a = operand(luaState, 1);
  const NodeId b = operand(luaState, 2);
  budgetCheck(luaState);
  return pushSym(luaState, arena(luaState).pw(F, a, b));
}

// `%` is Lua's FLOOR modulo, which is fmod pulled onto the sign of the divisor.
// Fn::Mod must carry exactly that: over 121 sign/magnitude pairs the obvious
// a - floor(a/b)*b disagrees with Lua on 14 of them, and the interpreter would
// then drift from the kernel with no visible cause.
int mmMod(lua_State* luaState) { return mmBinary<Fn::Mod>(luaState); }

int mmRejectCompare(lua_State* luaState) {
  const char* op = lua_tostring(luaState, lua_upvalueindex(1));
  const char* replacement = lua_tostring(luaState, lua_upvalueindex(2));
  fail(luaState,
       Cause::UntracedOperator,
       "`%s` on a traced value: the VM casts the result to a boolean, so this would silently "
       "collapse to one branch. Build the condition with ssol.%s(a, b) and consume it with "
       "ssol.select(c, a, b).",
       op,
       replacement);
}

int mmRejectNoMeaning(lua_State* luaState) {
  const char* op = lua_tostring(luaState, lua_upvalueindex(1));
  fail(luaState, Cause::UntracedOperator, "`%s` has no meaning on a traced value", op);
}

int mmTostring(lua_State* luaState) {
  const Sym* sym = testSym(luaState, 1);
  lua_pushfstring(luaState, "<traced #%d>", static_cast<int>(sym->id));
  return 1;
}

void rejectCompare(lua_State* luaState, const char* event, const char* op, const char* fn) {
  lua_pushstring(luaState, op);
  lua_pushstring(luaState, fn);
  lua_pushcclosure(luaState, mmRejectCompare, 2);
  lua_setfield(luaState, -2, event);
}

void rejectMeta(lua_State* luaState, const char* event, const char* op) {
  lua_pushstring(luaState, op);
  lua_pushcclosure(luaState, mmRejectNoMeaning, 1);
  lua_setfield(luaState, -2, event);
}

void registerSymMetatable(lua_State* luaState) {
  luaL_newmetatable(luaState, SymMetaName);
  const struct {
    const char* event;
    lua_CFunction fn;
  } arith[] = {{"__add", mmBinary<Fn::Add>},
               {"__sub", mmBinary<Fn::Sub>},
               {"__mul", mmBinary<Fn::Mul>},
               {"__div", mmBinary<Fn::Div>},
               {"__pow", mmBinary<Fn::Pow>},
               {"__mod", mmMod},
               {"__unm", mmUnary<Fn::Neg>},
               {"__tostring", mmTostring}};
  for (const auto& entry : arith) {
    lua_pushcfunction(luaState, entry.fn);
    lua_setfield(luaState, -2, entry.event);
  }

  // Net 1. __eq only fires when BOTH operands are full userdata and they are
  // not primitively equal, so `x == 0` against a number never reaches here --
  // that case belongs to net 3.
  rejectCompare(luaState, "__lt", "<", "lt");
  rejectCompare(luaState, "__le", "<=", "le");
  rejectCompare(luaState, "__eq", "==", "eq");

  for (const auto* op : {"__idiv",
                         "__concat",
                         "__len",
                         "__band",
                         "__bor",
                         "__bxor",
                         "__bnot",
                         "__shl",
                         "__shr",
                         "__index",
                         "__newindex",
                         "__call"}) {
    rejectMeta(luaState, op, op + 2); // strip the "__"
  }
  lua_pop(luaState, 1);
}

// ---------------------------------------------------------------------------
// ssol library, symbolic flavour
// ---------------------------------------------------------------------------
int compare(lua_State* luaState, Fn fn, bool swap) {
  const NodeId a = operand(luaState, swap ? 2 : 1);
  const NodeId b = operand(luaState, swap ? 1 : 2);
  const NodeId id = arena(luaState).pw(fn, a, b);
  state(luaState)->boolProduced.push_back(id);
  return pushSym(luaState, id, true);
}

int ssolLt(lua_State* l) { return compare(l, Fn::Lt, false); }
int ssolLe(lua_State* l) { return compare(l, Fn::Le, false); }
int ssolGt(lua_State* l) { return compare(l, Fn::Lt, true); }
int ssolGe(lua_State* l) { return compare(l, Fn::Le, true); }
int ssolEq(lua_State* l) { return compare(l, Fn::Eq, false); }

NodeId condition(lua_State* luaState, int index, const char* who) {
  if (Sym* sym = testSym(luaState, index)) {
    if (!sym->isBool) {
      fail(luaState,
           Cause::UntracedOperator,
           "%s: argument %d is a value, not a condition; build one with ssol.lt/le/gt/ge/eq",
           who,
           index);
    }
    state(luaState)->boolConsumed.push_back(sym->id);
    return sym->id;
  }
  if (lua_isboolean(luaState, index)) {
    return arena(luaState).konst(lua_toboolean(luaState, index) != 0 ? 1.0 : 0.0);
  }
  if (lua_type(luaState, index) == LUA_TNUMBER) {
    return operand(luaState, index);
  }
  fail(luaState,
       Cause::UntracedOperator,
       "%s: condition expected, got %s",
       who,
       luaL_typename(luaState, index));
}

// Conditions are 0.0/1.0, so the logical connectives are branchless: And/Or are
// Min/Max and negation is 1-c. That identity is what lets Select stay the only
// control-flow-shaped node in the IR.
int ssolAnd(lua_State* luaState) {
  const NodeId a = condition(luaState, 1, "ssol.land");
  const NodeId b = condition(luaState, 2, "ssol.land");
  return pushSym(luaState, arena(luaState).pw(Fn::And, a, b), true);
}
int ssolOr(lua_State* luaState) {
  const NodeId a = condition(luaState, 1, "ssol.lor");
  const NodeId b = condition(luaState, 2, "ssol.lor");
  return pushSym(luaState, arena(luaState).pw(Fn::Or, a, b), true);
}
int ssolNot(lua_State* luaState) {
  const NodeId c = condition(luaState, 1, "ssol.lnot");
  const NodeId one = arena(luaState).konst(1.0);
  return pushSym(luaState, arena(luaState).pw(Fn::Sub, one, c), true);
}
int ssolSelect(lua_State* luaState) {
  const NodeId c = condition(luaState, 1, "ssol.select");
  const NodeId a = operand(luaState, 2);
  const NodeId b = operand(luaState, 3);
  budgetCheck(luaState);
  return pushSym(luaState, arena(luaState).pw(Fn::Select, c, a, b));
}

// ---------------------------------------------------------------------------
// Shadow math
// ---------------------------------------------------------------------------
int mathVariadic(lua_State* luaState, Fn fn) {
  const int n = lua_gettop(luaState);
  luaL_argcheck(luaState, n >= 1, 1, "at least one argument");
  NodeId acc = operand(luaState, 1);
  for (int i = 2; i <= n; ++i) {
    acc = arena(luaState).pw(fn, acc, operand(luaState, i));
  }
  budgetCheck(luaState);
  return pushSym(luaState, acc);
}
int mathMin(lua_State* l) { return mathVariadic(l, Fn::Min); }
int mathMax(lua_State* l) { return mathVariadic(l, Fn::Max); }

int mathAtan(lua_State* luaState) {
  if (lua_gettop(luaState) >= 2 && !lua_isnoneornil(luaState, 2)) {
    return mmBinary<Fn::Atan2>(luaState);
  }
  return mmUnary<Fn::Atan>(luaState);
}

// Deliberately NOT lowered. exp(x)-1 and log(1+x) lose exactly the precision
// those functions exist for, so a lowered version would disagree with the
// interpreted reader precisely where someone reached for them.
int mathUnsupported(lua_State* luaState) {
  fail(luaState,
       Cause::UntracedOperator,
       "math.%s is not traceable: there is no IR op for it, and lowering it onto the ops that "
       "exist would change the result",
       lua_tostring(luaState, lua_upvalueindex(1)));
}

void buildShadowMath(lua_State* luaState) {
  lua_newtable(luaState);
  const struct {
    const char* name;
    lua_CFunction fn;
  } entries[] = {{"sqrt", mmUnary<Fn::Sqrt>},
                 {"abs", mmUnary<Fn::Abs>},
                 {"exp", mmUnary<Fn::Exp>},
                 {"log", mmUnary<Fn::Log>},
                 {"log2", mmUnary<Fn::Log2>},
                 {"log10", mmUnary<Fn::Log10>},
                 {"floor", mmUnary<Fn::Floor>},
                 {"ceil", mmUnary<Fn::Ceil>},
                 {"round", mmUnary<Fn::Round>},
                 {"sin", mmUnary<Fn::Sin>},
                 {"cos", mmUnary<Fn::Cos>},
                 {"tan", mmUnary<Fn::Tan>},
                 {"asin", mmUnary<Fn::Asin>},
                 {"acos", mmUnary<Fn::Acos>},
                 {"sinh", mmUnary<Fn::Sinh>},
                 {"cosh", mmUnary<Fn::Cosh>},
                 {"tanh", mmUnary<Fn::Tanh>},
                 {"asinh", mmUnary<Fn::Asinh>},
                 {"acosh", mmUnary<Fn::Acosh>},
                 {"atanh", mmUnary<Fn::Atanh>},
                 {"erf", mmUnary<Fn::Erf>},
                 {"atan", mathAtan},
                 {"pow", mmBinary<Fn::Pow>},
                 {"min", mathMin},
                 {"max", mathMax},
                 {"fmin", mathMin},
                 {"fmax", mathMax}};
  for (const auto& entry : entries) {
    lua_pushcfunction(luaState, entry.fn);
    lua_setfield(luaState, -2, entry.name);
  }

  // math.fmod is C truncation; `%` is floor. Two different functions, and only
  // the second has an op. math.random/randomseed are absent for a different
  // reason: a traced program must not be able to bake in a per-run value.
  for (const auto* name :
       {"fmod",     "expm1",     "log1p",     "hypot",     "cbrt",   "gamma",     "lgamma", "erfc",
        "copysign", "trunc",     "nextafter", "remainder", "fma",    "exp2",      "frexp",  "ldexp",
        "scalbn",   "tointeger", "type",      "ult",       "random", "randomseed"}) {
    lua_pushstring(luaState, name);
    lua_pushcclosure(luaState, mathUnsupported, 1);
    lua_setfield(luaState, -2, name);
  }

  lua_pushnumber(luaState, M_PI);
  lua_setfield(luaState, -2, "pi");
  lua_pushnumber(luaState, HUGE_VAL);
  lua_setfield(luaState, -2, "huge");
}

// ---------------------------------------------------------------------------
// Grid sampling: fields.<name>:sample(coords...) -> one Lookup per component
// ---------------------------------------------------------------------------
int traceSample(lua_State* luaState) {
  TraceState* s = state(luaState);
  const auto slot = static_cast<std::size_t>(lua_tointeger(luaState, lua_upvalueindex(1)));
  const int nCoords = lua_gettop(luaState) - 1; // drop `self`
  if (nCoords < 1 || nCoords > 6) {
    fail(luaState,
         Cause::UntracedOperator,
         "field `%s`: sampled with %d coordinates, the IR allows 1..6",
         s->specs[slot].name.c_str(),
         nCoords);
  }
  std::vector<NodeId> coords;
  coords.reserve(static_cast<std::size_t>(nCoords));
  for (int i = 0; i < nCoords; ++i) {
    coords.push_back(operand(luaState, 2 + i));
  }

  const auto components = static_cast<std::int32_t>(s->specs[slot].parameters.size());
  for (std::int32_t component = 0; component < components; ++component) {
    pushSym(luaState, arena(luaState).lookup(s->gridIds[slot], component, coords));
  }
  budgetCheck(luaState);
  return static_cast<int>(components);
}

// ---------------------------------------------------------------------------
// Deterministic pairs
// ---------------------------------------------------------------------------
// Lua seeds its string hash per state from an address and the clock, so an
// unshadowed `pairs` walks a table in a different order in every process. A
// script that accumulates over a table then gets a different summation order --
// hence a different fingerprint AND different rounding -- on every rank, which
// under a domain decomposition is two ranks computing different values at the
// same point.
int sortedNext(lua_State* luaState) {
  const int n = static_cast<int>(lua_tointeger(luaState, lua_upvalueindex(2)));
  const int i = static_cast<int>(lua_tointeger(luaState, lua_upvalueindex(3)));
  if (i >= n) {
    lua_pushnil(luaState);
    return 1;
  }
  lua_pushinteger(luaState, i + 1);
  lua_replace(luaState, lua_upvalueindex(3));
  lua_rawgeti(luaState, lua_upvalueindex(1), i + 1); // key
  lua_pushvalue(luaState, -1);
  lua_gettable(luaState, 1); // value
  return 2;
}

bool keyLess(lua_State* luaState, int lhs, int rhs) {
  const int tl = lua_type(luaState, lhs);
  const int tr = lua_type(luaState, rhs);
  if (tl != tr) {
    return tl == LUA_TNUMBER;
  }
  if (tl == LUA_TNUMBER) {
    return lua_tonumber(luaState, lhs) < lua_tonumber(luaState, rhs);
  }
  if (tl == LUA_TSTRING) {
    return std::strcmp(lua_tostring(luaState, lhs), lua_tostring(luaState, rhs)) < 0;
  }
  return false;
}

int deterministicPairs(lua_State* luaState) {
  luaL_checktype(luaState, 1, LUA_TTABLE);
  lua_newtable(luaState); // key array
  int n = 0;
  lua_pushnil(luaState);
  while (lua_next(luaState, 1) != 0) {
    lua_pop(luaState, 1); // value
    lua_pushvalue(luaState, -1);
    lua_rawseti(luaState, -3, ++n);
  }
  // Insertion sort: key counts are tiny here and it keeps the ordering rule in
  // one readable place.
  for (int i = 2; i <= n; ++i) {
    lua_rawgeti(luaState, -1, i);
    int j = i - 1;
    while (j >= 1) {
      lua_rawgeti(luaState, -2, j);
      const bool less = keyLess(luaState, -2, -1);
      if (!less) {
        lua_pop(luaState, 1);
        break;
      }
      lua_rawseti(luaState, -3, j + 1);
      --j;
    }
    lua_rawseti(luaState, -2, j + 1);
  }
  lua_pushinteger(luaState, n);
  lua_pushinteger(luaState, 0);
  lua_pushcclosure(luaState, sortedNext, 3);
  lua_pushvalue(luaState, 1);
  lua_pushnil(luaState);
  return 3;
}

// ---------------------------------------------------------------------------
// Concrete ssol, for the probe runs
// ---------------------------------------------------------------------------
bool truthy(lua_State* luaState, int index) {
  if (lua_isboolean(luaState, index)) {
    return lua_toboolean(luaState, index) != 0;
  }
  return lua_tonumber(luaState, index) != 0.0;
}
int concreteLt(lua_State* l) {
  lua_pushboolean(l, lua_tonumber(l, 1) < lua_tonumber(l, 2));
  return 1;
}
int concreteLe(lua_State* l) {
  lua_pushboolean(l, lua_tonumber(l, 1) <= lua_tonumber(l, 2));
  return 1;
}
int concreteGt(lua_State* l) {
  lua_pushboolean(l, lua_tonumber(l, 1) > lua_tonumber(l, 2));
  return 1;
}
int concreteGe(lua_State* l) {
  lua_pushboolean(l, lua_tonumber(l, 1) >= lua_tonumber(l, 2));
  return 1;
}
int concreteEq(lua_State* l) {
  lua_pushboolean(l, lua_tonumber(l, 1) == lua_tonumber(l, 2));
  return 1;
}
int concreteAnd(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(truthy(l, 1) && truthy(l, 2)));
  return 1;
}
int concreteOr(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(truthy(l, 1) || truthy(l, 2)));
  return 1;
}
int concreteNot(lua_State* l) {
  lua_pushboolean(l, static_cast<int>(!truthy(l, 1)));
  return 1;
}
// Branches inside a C function, so it emits no line event and cannot itself
// produce a false positive in the probe diff.
int concreteSelect(lua_State* l) {
  lua_pushvalue(l, truthy(l, 1) ? 2 : 3);
  return 1;
}

// ---------------------------------------------------------------------------
// Environment
// ---------------------------------------------------------------------------
void pushSsolTableImpl(lua_State* luaState, bool symbolic);

void pushSsolTable(lua_State* luaState, bool symbolic) { pushSsolTableImpl(luaState, symbolic); }

void buildEnvironment(lua_State* luaState, bool symbolic) {
  lua_newtable(luaState);
  for (const auto* global : {"assert",
                             "error",
                             "ipairs",
                             "next",
                             "select",
                             "tonumber",
                             "tostring",
                             "type",
                             "rawequal",
                             "rawlen",
                             "table",
                             "string"}) {
    lua_getglobal(luaState, global);
    if (lua_isnil(luaState, -1)) {
      lua_pop(luaState, 1);
      continue;
    }
    lua_setfield(luaState, -2, global);
  }
  lua_pushcfunction(luaState, deterministicPairs);
  lua_setfield(luaState, -2, "pairs");

  if (symbolic) {
    buildShadowMath(luaState);
  } else {
    lua_getglobal(luaState, "math");
  }
  lua_setfield(luaState, -2, "math");

  pushSsolTable(luaState, symbolic);
  lua_setfield(luaState, -2, "ssol");

  lua_pushvalue(luaState, -1);
  lua_setfield(luaState, -2, "_G");
}

void pushSsolTableImpl(lua_State* luaState, bool symbolic) {
  lua_newtable(luaState);
  const struct {
    const char* name;
    lua_CFunction symbolicFn;
    lua_CFunction concreteFn;
  } ssol[] = {{"lt", ssolLt, concreteLt},
              {"le", ssolLe, concreteLe},
              {"gt", ssolGt, concreteGt},
              {"ge", ssolGe, concreteGe},
              {"eq", ssolEq, concreteEq},
              {"land", ssolAnd, concreteAnd},
              {"lor", ssolOr, concreteOr},
              {"lnot", ssolNot, concreteNot},
              {"select", ssolSelect, concreteSelect}};
  for (const auto& entry : ssol) {
    // The probe runs need the same names computing CONCRETELY and returning
    // real booleans. That is exactly what makes `if ssol.gt(z, 0) then` show up
    // as a divergent line trace -- the symbolic run cannot see it, because a
    // userdata condition is unconditionally truthy.
    //
    // The concrete column is also what the interpreted LuaReader installs, so
    // the fallback path evaluates the identical functions rather than a second
    // implementation that agrees with this one only by inspection.
    lua_pushcfunction(luaState, symbolic ? entry.symbolicFn : entry.concreteFn);
    lua_setfield(luaState, -2, entry.name);
  }
}

// os, io, debug, require, load, dofile and package are absent by construction:
// the environment is built up from a whitelist rather than stripped down from
// the globals, so a library added to luaL_openlibs later cannot leak in.
lua_State* loadModule(const std::string& code, bool symbolic, std::string& error) {
  lua_State* luaState = luaL_newstate();
  if (luaState == nullptr) {
    error = "could not create a Lua state";
    return nullptr;
  }
  luaL_openlibs(luaState);
  registerSymMetatable(luaState);

  if (luaL_loadbufferx(luaState, code.c_str(), code.size(), "@model.lua", "t") != LUA_OK) {
    error = lua_tostring(luaState, -1);
    lua_close(luaState);
    return nullptr;
  }
  buildEnvironment(luaState, symbolic);
  lua_setupvalue(luaState, -2, 1); // _ENV

  if (lua_pcall(luaState, 0, 1, 0) != LUA_OK) {
    error = lua_tostring(luaState, -1);
    lua_close(luaState);
    return nullptr;
  }
  if (!lua_istable(luaState, -1)) {
    error = "the module did not return a table";
    lua_close(luaState);
    return nullptr;
  }
  return luaState;
}

// ---------------------------------------------------------------------------
// Module inspection
// ---------------------------------------------------------------------------
std::vector<std::string> stringArray(lua_State* luaState, int index) {
  std::vector<std::string> out;
  if (!lua_istable(luaState, index)) {
    return out;
  }
  const lua_Integer n = luaL_len(luaState, index);
  out.reserve(static_cast<std::size_t>(n));
  for (lua_Integer i = 1; i <= n; ++i) {
    lua_geti(luaState, index, i);
    if (lua_type(luaState, -1) == LUA_TSTRING) {
      out.emplace_back(lua_tostring(luaState, -1));
    }
    lua_pop(luaState, 1);
  }
  return out;
}

std::string stringField(lua_State* luaState, int index, const char* name) {
  std::string out;
  lua_getfield(luaState, index, name);
  if (lua_type(luaState, -1) == LUA_TSTRING) {
    out = lua_tostring(luaState, -1);
  }
  lua_pop(luaState, 1);
  return out;
}

std::vector<FieldSpec> readFieldSpecs(lua_State* luaState, int moduleIndex) {
  std::vector<FieldSpec> out;
  lua_getfield(luaState, moduleIndex, "field_specs");
  if (lua_istable(luaState, -1)) {
    const lua_Integer n = luaL_len(luaState, -1);
    for (lua_Integer i = 1; i <= n; ++i) {
      lua_geti(luaState, -1, i);
      if (lua_istable(luaState, -1)) {
        FieldSpec spec;
        spec.name = stringField(luaState, -1, "name");
        spec.kind = stringField(luaState, -1, "kind");
        spec.file = stringField(luaState, -1, "file");
        spec.variable = stringField(luaState, -1, "variable");
        spec.interpolation = stringField(luaState, -1, "interpolation");
        lua_getfield(luaState, -1, "parameters");
        spec.parameters = stringArray(luaState, -1);
        lua_pop(luaState, 1);
        out.push_back(std::move(spec));
      }
      lua_pop(luaState, 1);
    }
  }
  lua_pop(luaState, 1);
  return out;
}

struct ParamInfo {
  std::vector<std::string> names;
  int declared{0};
  bool vararg{false};
};

// Verified behaviour on 5.4.6: locals are NOT leaked, only parameters; method
// syntax inserts `self` first; a vararg function reports zero parameters; a
// stripped chunk yields no names while lua_getinfo still reports the count.
ParamInfo parameterNames(lua_State* luaState, int functionIndex) {
  ParamInfo info;
  lua_Debug debug;
  lua_pushvalue(luaState, functionIndex);
  lua_getinfo(luaState, ">u", &debug); // pops the function
  info.declared = debug.nparams;
  info.vararg = debug.isvararg != 0;

  lua_pushvalue(luaState, functionIndex);
  for (int i = 1;; ++i) {
    const char* name = lua_getlocal(luaState, nullptr, i);
    if (name == nullptr) {
      break;
    }
    info.names.emplace_back(name);
  }
  lua_pop(luaState, 1);
  return info;
}

// A module-scope local is an upvalue, not a table field, so no __newindex on M
// would ever see it being written. Comparing the scalar upvalues around the
// trace does.
std::vector<double> upvalueSnapshot(lua_State* luaState, int functionIndex) {
  std::vector<double> out;
  for (int i = 1;; ++i) {
    const char* name = lua_getupvalue(luaState, functionIndex, i);
    if (name == nullptr) {
      break;
    }
    if (std::strcmp(name, "_ENV") != 0) {
      if (lua_type(luaState, -1) == LUA_TNUMBER) {
        out.push_back(lua_tonumber(luaState, -1));
      } else if (lua_isboolean(luaState, -1)) {
        out.push_back(lua_toboolean(luaState, -1) != 0 ? 1.0 : 0.0);
      }
    }
    lua_pop(luaState, 1);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Net 3: probe runs
// ---------------------------------------------------------------------------
struct ProbeRecord {
  std::vector<int> lines;
  bool overflow{false};
};

// One hook cannot take a context pointer, and the probe is single-threaded and
// short-lived, so a thread-local is the honest mechanism here.
thread_local ProbeRecord* activeProbe = nullptr;

void lineHook(lua_State* luaState, lua_Debug* debug) {
  if (activeProbe == nullptr) {
    return;
  }
  if (activeProbe->lines.size() > (1U << 20U)) {
    activeProbe->overflow = true;
    return;
  }
  lua_getinfo(luaState, "l", debug);
  activeProbe->lines.push_back(debug->currentline);
}

// The ladder, not the count, is what decides whether a branch is seen: a probe
// only observes a threshold it straddles. A first version used four
// "interesting" vectors and missed `if z > -1000` because every draw sat above
// it. Any threshold strictly between two neighbouring rungs is caught; one
// beyond the ends is not, which is why net 4 exists.
const std::vector<double>& probeLadder() {
  static const std::vector<double> rungs = {
      -1e9, -1e6, -1e3, -1.0, -1e-3, 0.0, 1e-3, 1.0, 1e3, 1e6, 1e9};
  return rungs;
}

int probeSample(lua_State* luaState) {
  const auto components = static_cast<int>(lua_tointeger(luaState, lua_upvalueindex(1)));
  double accumulated = 0.0;
  for (int i = 2; i <= lua_gettop(luaState); ++i) {
    accumulated += lua_tonumber(luaState, i);
  }
  for (int component = 0; component < components; ++component) {
    lua_pushnumber(luaState, accumulated + component);
  }
  return components;
}

void pushProbeFields(lua_State* luaState, const std::vector<FieldSpec>& specs) {
  lua_newtable(luaState);
  for (const auto& spec : specs) {
    lua_newtable(luaState);
    lua_pushinteger(luaState, static_cast<lua_Integer>(spec.parameters.size()));
    lua_pushcclosure(luaState, probeSample, 1);
    lua_setfield(luaState, -2, "sample");
    lua_setfield(luaState, -2, spec.name.c_str());
  }
}

bool probe(const std::string& code,
           std::size_t inputCount,
           const std::vector<FieldSpec>& specs,
           TraceFailure& failure) {
  std::vector<ProbeRecord> records;
  records.reserve(probeLadder().size());

  for (const double rung : probeLadder()) {
    std::string error;
    lua_State* luaState = loadModule(code, /*symbolic=*/false, error);
    if (luaState == nullptr) {
      failure.cause = Cause::LoadError;
      failure.reason = "probe run failed to load the module: " + error;
      return false;
    }
    lua_getfield(luaState, -1, "evaluate");
    pushProbeFields(luaState, specs);
    for (std::size_t i = 0; i < inputCount; ++i) {
      lua_pushnumber(luaState, rung);
    }

    ProbeRecord record;
    activeProbe = &record;
    lua_sethook(luaState, lineHook, LUA_MASKLINE, 0);
    const int status = lua_pcall(luaState, static_cast<int>(inputCount) + 1, LUA_MULTRET, 0);
    lua_sethook(luaState, nullptr, 0, 0);
    activeProbe = nullptr;
    // A probe that errors is not evidence of branching -- a domain error at one
    // rung says nothing about control flow -- so it contributes no trace.
    if (status != LUA_OK) {
      record.lines.clear();
    }
    lua_close(luaState);
    records.push_back(std::move(record));
  }

  for (std::size_t i = 1; i < records.size(); ++i) {
    if (records[i].lines == records[0].lines) {
      continue;
    }
    std::size_t k = 0;
    while (k < records[0].lines.size() && k < records[i].lines.size() &&
           records[0].lines[k] == records[i].lines[k]) {
      ++k;
    }
    const int line = k < records[0].lines.size()   ? records[0].lines[k]
                     : k < records[i].lines.size() ? records[i].lines[k]
                                                   : -1;
    failure.cause = Cause::DataDependentFlow;
    failure.line = line;
    failure.reason = "control flow depends on the input values (the probes at " +
                     std::to_string(probeLadder()[0]) + " and " + std::to_string(probeLadder()[i]) +
                     " diverge near line " + std::to_string(line) +
                     "); a traced program has one shape for every point";
    return false;
  }
  return true;
}

} // namespace

void registerConcreteSsol(lua_State* luaState) {
  pushSsolTable(luaState, /*symbolic=*/false);
  lua_setglobal(luaState, "ssol");
}

// ---------------------------------------------------------------------------
std::optional<expr::Program>
    traceLuaModule(const std::string& code, const TraceOptions& options, TraceFailure& failure) {
  expr::Program program;
  TraceState traceState;
  traceState.program = &program;
  traceState.options = &options;

  std::string error;
  lua_State* luaState = loadModule(code, /*symbolic=*/true, error);
  if (luaState == nullptr) {
    failure.cause = Cause::LoadError;
    failure.reason = error;
    return std::nullopt;
  }
  const int moduleIndex = lua_gettop(luaState);

  auto abort = [&](Cause cause, std::string reason) {
    failure.cause = cause;
    failure.reason = std::move(reason);
    failure.line = traceState.line;
    lua_close(luaState);
    return std::nullopt;
  };

  traceState.specs = readFieldSpecs(luaState, moduleIndex);
  for (const auto& spec : traceState.specs) {
    if (spec.parameters.empty()) {
      return abort(Cause::LoadError,
                   "field `" + spec.name +
                       "` declares no parameters; the tracer needs the component count before it "
                       "can emit the lookups, so `parameters` is required for every kind");
    }
    const auto kind = reader::datafield::parseGridKind(spec.kind);
    if (!kind.has_value()) {
      return abort(Cause::LoadError,
                   "field `" + spec.name + "` declares unknown kind \"" + spec.kind + "\"");
    }
    // An absent `variable` is the documented default, not a diagnostic; an
    // unknown interpolation is the reverse. Distinguishing the two is why the
    // parser returns optional instead of falling back to Linear.
    const std::string interpolationName =
        spec.interpolation.empty() ? std::string("linear") : spec.interpolation;
    const auto interpolation = reader::datafield::parseInterpolation(interpolationName);
    if (!interpolation.has_value()) {
      return abort(Cause::LoadError,
                   "field `" + spec.name + "` declares unknown interpolation \"" +
                       spec.interpolation + "\"");
    }

    reader::datafield::GridDesc desc;
    // Grid.h is the contract. `parameters` stays local: it gives the tracer the
    // component count it needs to emit the lookups, but Kind::Lookup carries a
    // resolved integer, and Grid.h declines to hold a second name table for the
    // same reason Program.h infers the signature instead of trusting it.
    //
    // NOT YET EXPRESSIBLE in field_specs: boundary and timeAxis, both left at
    // their defaults (Clamp, static), exactly as on the sderiv side.
    desc.kind = *kind;
    desc.path = spec.file;
    desc.variable = spec.variable.empty() ? reader::datafield::DefaultGridVariable : spec.variable;
    desc.interpolation = *interpolation;
    traceState.gridIds.push_back(program.internGrid(desc));
  }

  lua_pushlightuserdata(luaState, &traceState);
  lua_setfield(luaState, LUA_REGISTRYINDEX, TracerRegistryKey);

  lua_getfield(luaState, moduleIndex, "evaluate");
  if (lua_isfunction(luaState, -1) == 0) {
    return abort(Cause::NoEvaluate, "the module has no `evaluate` function");
  }
  const int functionIndex = lua_gettop(luaState);

  const ParamInfo params = parameterNames(luaState, functionIndex);
  if (params.vararg) {
    return abort(Cause::NoEvaluate,
                 "`evaluate` is variadic; the input signature is read from its parameter names, "
                 "so `...` would trace as zero inputs and yield a constant program");
  }
  if (static_cast<int>(params.names.size()) != params.declared) {
    return abort(Cause::NoEvaluate,
                 "`evaluate` declares " + std::to_string(params.declared) +
                     " parameters but carries no debug names; was the chunk stripped?");
  }
  if (params.names.empty() || params.names.front() != "fields") {
    return abort(Cause::NoEvaluate,
                 std::string("`evaluate`'s first parameter must be named `fields`, found `") +
                     (params.names.empty() ? "(none)" : params.names.front()) +
                     "`; note that `function M:evaluate(...)` silently inserts `self` first");
  }
  // Script parameters are folded to Const here rather than becoming channels.
  // The choice deliberately does NOT live at trace time for the caller: tracing
  // happens once per script either way, and a caller that wants one kernel
  // across parametrisations simply leaves options.parameters empty and lets the
  // names stay ordinary inputs.
  std::vector<std::string> inputNames;
  std::vector<NodeId> argumentNodes;
  for (auto it = params.names.begin() + 1; it != params.names.end(); ++it) {
    const auto folded = std::find_if(options.parameters.begin(),
                                     options.parameters.end(),
                                     [&](const auto& p) { return p.first == *it; });
    if (folded != options.parameters.end()) {
      argumentNodes.push_back(program.arena().konst(folded->second));
    } else {
      argumentNodes.push_back(program.arena().field(*it));
      inputNames.push_back(*it);
    }
  }
  for (const auto& [name, value] : options.parameters) {
    if (std::find(params.names.begin() + 1, params.names.end(), name) == params.names.end()) {
      return abort(Cause::SignatureMismatch,
                   "parameter `" + name + "` was supplied but `evaluate` has no such parameter");
    }
  }
  const std::vector<double> upvaluesBefore = upvalueSnapshot(luaState, functionIndex);

  // fields table: one object per grid, `sample` closing over the spec slot.
  lua_newtable(luaState);
  for (std::size_t slot = 0; slot < traceState.specs.size(); ++slot) {
    lua_newtable(luaState);
    lua_pushinteger(luaState, static_cast<lua_Integer>(slot));
    lua_pushcclosure(luaState, traceSample, 1);
    lua_setfield(luaState, -2, "sample");
    lua_setfield(luaState, -2, traceState.specs[slot].name.c_str());
  }
  for (const NodeId node : argumentNodes) {
    pushSym(luaState, node);
  }

  const int base = functionIndex - 1;
  if (lua_pcall(luaState, static_cast<int>(argumentNodes.size()) + 1, LUA_MULTRET, 0) != LUA_OK) {
    const std::string message = lua_tostring(luaState, -1);
    return abort(traceState.cause, message);
  }

  const int returned = lua_gettop(luaState) - base;
  std::vector<NodeId> roots;
  roots.reserve(static_cast<std::size_t>(returned));
  for (int i = 1; i <= returned; ++i) {
    const int index = base + i;
    if (const Sym* sym = testSym(luaState, index)) {
      roots.push_back(sym->id);
    } else if (lua_type(luaState, index) == LUA_TNUMBER) {
      // A constant output is legal: a model may pin one quantity.
      roots.push_back(program.arena().konst(lua_tonumber(luaState, index)));
    } else {
      return abort(Cause::SignatureMismatch,
                   "`evaluate` returned a " + std::string(luaL_typename(luaState, index)) +
                       " at position " + std::to_string(i) + "; only numbers are outputs");
    }
  }

  lua_settop(luaState, moduleIndex);
  lua_getfield(luaState, moduleIndex, "evaluate");
  if (upvalueSnapshot(luaState, lua_gettop(luaState)) != upvaluesBefore) {
    return abort(Cause::SideEffect,
                 "`evaluate` mutated module state during the trace; a traced program is evaluated "
                 "once per point in unspecified order, so state carried between calls cannot mean "
                 "what it looks like it means");
  }
  lua_pop(luaState, 1);

  lua_getfield(luaState, moduleIndex, "output_parameters");
  const std::vector<std::string> outputNames = stringArray(luaState, -1);
  lua_pop(luaState, 1);
  if (outputNames.size() != roots.size()) {
    return abort(Cause::SignatureMismatch,
                 "`output_parameters` declares " + std::to_string(outputNames.size()) +
                     " names but `evaluate` returned " + std::to_string(roots.size()) + " values");
  }

  // Net 2.
  for (const NodeId id : traceState.boolProduced) {
    const bool consumed =
        std::find(traceState.boolConsumed.begin(), traceState.boolConsumed.end(), id) !=
            traceState.boolConsumed.end() ||
        std::find(roots.begin(), roots.end(), id) != roots.end();
    if (!consumed) {
      return abort(Cause::RawCondition,
                   "a condition was built but never passed to ssol.select/land/lor/lnot; a "
                   "traced condition placed in a raw `if` is always true, and the other branch "
                   "would be silently traced away");
    }
  }

  lua_close(luaState);

  // Net 3. Runs on a fresh state with the real math library, so it must come
  // after the symbolic state is closed.
  if (options.probe && !probe(code, argumentNodes.size(), traceState.specs, failure)) {
    return std::nullopt;
  }

  for (const auto& name : inputNames) {
    program.addInput(name, DataType::F64);
  }
  for (std::size_t i = 0; i < outputNames.size(); ++i) {
    program.addOutput(outputNames[i], DataType::F64, roots[i]);
  }
  try {
    expr::validate(program);
  } catch (const std::invalid_argument& e) {
    // A structural failure here is a tracer bug, not a script bug, but the
    // caller still has to fall back rather than abort a running job.
    failure.cause = Cause::LoadError;
    failure.reason = std::string("the traced program failed validation: ") + e.what();
    return std::nullopt;
  }
  return program;
}

} // namespace seissol::reader::scripting

#endif // USE_LUA
