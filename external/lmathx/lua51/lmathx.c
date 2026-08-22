// SPDX-FileCopyrightText: none
// SPDX-License-Identifier: CC0-1.0
/*
* lmathx.c
* C99 math functions for Lua
* Luiz Henrique de Figueiredo <lhf@tecgraf.puc-rio.br>
* 04 Apr 2010 22:53:32
* This code is hereby placed in the public domain.
*/

#define _GNU_SOURCE 1
#include <math.h>
#include <stdio.h>

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#define A(i)	luaL_checknumber(L,i)

static int Lacosh(lua_State *L)
{
  lua_pushnumber(L,acosh(A(1)));
  return 1;
}

static int Lasinh(lua_State *L)
{
  lua_pushnumber(L,asinh(A(1)));
  return 1;
}

static int Latanh(lua_State *L)
{
  lua_pushnumber(L,atanh(A(1)));
  return 1;
}

static int Lcbrt(lua_State *L)
{
  lua_pushnumber(L,cbrt(A(1)));
  return 1;
}

static int Lcopysign(lua_State *L)
{
  lua_pushnumber(L,copysign(A(1),A(2)));
  return 1;
}

static int Lerf(lua_State *L)
{
  lua_pushnumber(L,erf(A(1)));
  return 1;
}

static int Lerfc(lua_State *L)
{
  lua_pushnumber(L,erfc(A(1)));
  return 1;
}

static int Lexp2(lua_State *L)
{
  lua_pushnumber(L,exp2(A(1)));
  return 1;
}

static int Lexpm1(lua_State *L)
{
  lua_pushnumber(L,expm1(A(1)));
  return 1;
}

static int Lfdim(lua_State *L)
{
  lua_pushnumber(L,fdim(A(1),A(2)));
  return 1;
}

static int Lfma(lua_State *L)
{
  lua_pushnumber(L,fma(A(1),A(2),A(3)));
  return 1;
}

static int Lfmax(lua_State *L)
{
 int i,n=lua_gettop(L);
 lua_Number m=A(1);
 for (i=2; i<=n; i++) m=fmax(m,A(i));
 lua_pushnumber(L,m);
 return 1;
}

static int Lfmin(lua_State *L)
{
 int i,n=lua_gettop(L);
 lua_Number m=A(1);
 for (i=2; i<=n; i++) m=fmin(m,A(i));
 lua_pushnumber(L,m);
 return 1;
}

static int Lfpclassify(lua_State *L)
{
 switch (fpclassify(A(1)))
 {
  case FP_INFINITE:	lua_pushliteral(L,"inf");	break;
  case FP_NAN:		lua_pushliteral(L,"nan");	break;
  case FP_NORMAL:	lua_pushliteral(L,"normal");	break;
  case FP_SUBNORMAL:	lua_pushliteral(L,"subnormal");	break;
  case FP_ZERO:		lua_pushliteral(L,"zero");	break;
 }
 return 1;
}

static int Lhypot(lua_State *L)
{
  lua_pushnumber(L,hypot(A(1),A(2)));
  return 1;
}

static int Lisfinite(lua_State *L)
{
  lua_pushboolean(L,isfinite(A(1)));
  return 1;
}

static int Lisinf(lua_State *L)
{
  lua_pushboolean(L,isinf(A(1)));
  return 1;
}

static int Lisnan(lua_State *L)
{
  lua_pushboolean(L,isnan(A(1)));
  return 1;
}

static int Lisnormal(lua_State *L)
{
  lua_pushboolean(L,isnormal(A(1)));
  return 1;
}

static int Llgamma(lua_State *L)
{
  lua_pushnumber(L,lgamma(A(1)));
  return 1;
}

static int Llog1p(lua_State *L)
{
  lua_pushnumber(L,log1p(A(1)));
  return 1;
}

static int Llog2(lua_State *L)
{
  lua_pushnumber(L,log2(A(1)));
  return 1;
}

static int Llogb(lua_State *L)
{
  lua_pushnumber(L,logb(A(1)));
  return 1;
}

static int Lnearbyint(lua_State *L)
{
  lua_pushnumber(L,nearbyint(A(1)));
  return 1;
}

static int Lnextafter(lua_State *L)
{
  lua_pushnumber(L,nextafter(A(1),A(2)));
  return 1;
}

static int Lremainder(lua_State *L)
{
  lua_pushnumber(L,remainder(A(1),A(2)));
  return 1;
}

static int Lrint(lua_State *L)
{
  lua_pushnumber(L,rint(A(1)));
  return 1;
}

static int Lround(lua_State *L)
{
  lua_pushnumber(L,round(A(1)));
  return 1;
}

static int Lscalbn(lua_State *L)
{
  lua_pushnumber(L,scalbn(A(1),(int)A(2)));
  return 1;
}

static int Lsignbit(lua_State *L)
{
  lua_pushboolean(L,signbit(A(1)));
  return 1;
}

static int Ltgamma(lua_State *L)
{
  lua_pushnumber(L,tgamma(A(1)));
  return 1;
}

static int Ltrunc(lua_State *L)
{
  lua_pushnumber(L,trunc(A(1)));
  return 1;
}

static const luaL_Reg R[] =
{
	{ "acosh",	Lacosh },
	{ "asinh",	Lasinh },
	{ "atanh",	Latanh },
	{ "cbrt",	Lcbrt },
	{ "copysign",	Lcopysign },
	{ "erfc",	Lerfc },
	{ "erf",	Lerf },
	{ "exp2",	Lexp2 },
	{ "expm1",	Lexpm1 },
	{ "fdim",	Lfdim },
	{ "fma",	Lfma },
	{ "fmax",	Lfmax },
	{ "fmin",	Lfmin },
	{ "fpclassify",	Lfpclassify },
	{ "gamma",	Ltgamma },
	{ "hypot",	Lhypot },
	{ "isfinite",	Lisfinite },
	{ "isinf",	Lisinf },
	{ "isnan",	Lisnan },
	{ "isnormal",	Lisnormal },
	{ "lgamma",	Llgamma },
	{ "log1p",	Llog1p },
	{ "log2",	Llog2 },
	{ "logb",	Llogb },
	{ "nearbyint",	Lnearbyint },
	{ "nextafter",	Lnextafter },
	{ "remainder",	Lremainder },
	{ "rint",	Lrint },
	{ "round",	Lround },
	{ "scalbn",	Lscalbn },
	{ "signbit",	Lsignbit },
	{ "trunc",	Ltrunc },
	{ NULL,		NULL }
};


LUALIB_API int luaopen_mathx(lua_State *L)
{
 luaL_register(L,LUA_MATHLIBNAME,R);
 lua_pushnumber(L,INFINITY);
 lua_setfield(L,-2,"infinity");
 lua_pushnumber(L,NAN);
 lua_setfield(L,-2,"nan");
 return 1;
}
