/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Loop unrolling with templates.
 **/

#ifndef INITIALIZER_TREE_FOREACH_HPP_
#define INITIALIZER_TREE_FOREACH_HPP_
 
template<bool B> struct BoolToType {};
typedef BoolToType<true> ForEachContinue;
typedef BoolToType<false> ForEachTerminate;

// One argument
template<template<unsigned> class F, unsigned INDEX, unsigned COUNT, typename Arg1>
void ForEach(Arg1 arg1, ForEachContinue)
{
  F<INDEX>()(arg1);
  ForEach<F, INDEX+1, COUNT>(arg1, BoolToType<INDEX+1 < COUNT>());
}

template<template<unsigned> class, unsigned, unsigned, typename Arg1>
void ForEach(Arg1, ForEachTerminate) {}

template<template<unsigned> class F, unsigned COUNT, typename Arg1>
void ForEachClassMethod(Arg1 arg1)
{
  ForEach<F, 0, COUNT>(arg1, BoolToType<0 < COUNT>());
}

// Two arguments
template<template<unsigned> class F, unsigned INDEX, unsigned COUNT, typename Arg1, typename Arg2>
void ForEach(Arg1 arg1, Arg2 arg2, ForEachContinue)
{
  F<INDEX>()(arg1, arg2);
  ForEach<F, INDEX+1, COUNT>(arg1, arg2, BoolToType<INDEX+1 < COUNT>());
}

template<template<unsigned> class, unsigned, unsigned, typename Arg1, typename Arg2>
void ForEach(Arg1, Arg2, ForEachTerminate) {}

template<template<unsigned> class F, unsigned COUNT, typename Arg1, typename Arg2>
void ForEachClassMethod(Arg1 arg1, Arg2 arg2)
{
  ForEach<F, 0, COUNT>(arg1, arg2, BoolToType<0 < COUNT>());
}

// Three arguments
template<template<unsigned> class F, unsigned INDEX, unsigned COUNT, typename Arg1, typename Arg2, typename Arg3>
void ForEach(Arg1 arg1, Arg2 arg2, Arg3 arg3, ForEachContinue)
{
  F<INDEX>()(arg1, arg2, arg3);
  ForEach<F, INDEX+1, COUNT>(arg1, arg2, arg3, BoolToType<INDEX+1 < COUNT>());
}

template<template<unsigned> class, unsigned, unsigned, typename Arg1, typename Arg2, typename Arg3>
void ForEach(Arg1, Arg2, Arg3, ForEachTerminate) {}

template<template<unsigned> class F, unsigned COUNT, typename Arg1, typename Arg2, typename Arg3>
void ForEachClassMethod(Arg1 arg1, Arg2 arg2, Arg3 arg3)
{
  ForEach<F, 0, COUNT>(arg1, arg2, arg3, BoolToType<0 < COUNT>());
}

#endif
