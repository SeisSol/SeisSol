/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2019, SeisSol Group
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
 **/

#ifndef INITIALIZER_INTERFACEHELPER_H_
#define INITIALIZER_INTERFACEHELPER_H_

#include <easi/util/Magic.h>

namespace seissol {
  template<typename X>
  struct extract_type
  {
    typedef X type;
  };

  template<template<typename> typename F, typename X>
  struct extract_type<F<X>>
  {
    typedef X type;
  };
}

#define _LTSTREE_MEMBER_PTR(N, HANDLE_STRUCT, X) seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type* X = nullptr;
#define _LTSTREE_MEMBER_REF(N, HANDLE_STRUCT, X) seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type& X;
#define _LTSTREE_MEMBER_REF_CS(HANDLE_STRUCT, X) seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type& X
#define _LTSTREE_LOAD(N, HANDLE_STRUCT, X) X = tree.var(handleStruct.X);
#define _LTSTREE_INIT_MEMBER(HANDLE_STRUCT, X) X(X)
#define _LTSTREE_ACCESS(HANDLE_STRUCT, X) X[index]
#define _LTSTREE_LOOKUP(HANDLE_STRUCT, X) lut.lookup(handleStruct.X,meshId)
#define LTSTREE_GENERATE_INTERFACE(NAME, HANDLE_STRUCT, ...)  struct NAME { \
                                                                MAGIC_FOR_EACH(_LTSTREE_MEMBER_REF, HANDLE_STRUCT, __VA_ARGS__) \
                                                                NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_MEMBER_REF_CS, HANDLE_STRUCT, __VA_ARGS__)) \
                                                                  : MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_INIT_MEMBER, HANDLE_STRUCT, __VA_ARGS__) {} \
                                                                  template<typename T> \
                                                                  static NAME lookup(HANDLE_STRUCT const& handleStruct, T const& lut, unsigned meshId) { \
                                                                    return NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_LOOKUP, HANDLE_STRUCT, __VA_ARGS__)); \
                                                                  } \
                                                                struct Loader { \
                                                                  MAGIC_FOR_EACH(_LTSTREE_MEMBER_PTR, HANDLE_STRUCT, __VA_ARGS__) \
                                                                  template<typename T> \
                                                                  void load(HANDLE_STRUCT const & handleStruct, T& tree) { \
                                                                    MAGIC_FOR_EACH(_LTSTREE_LOAD, HANDLE_STRUCT, __VA_ARGS__) \
                                                                  } \
                                                                  NAME entry(unsigned index) { \
                                                                    return NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_ACCESS, HANDLE_STRUCT, __VA_ARGS__)); \
                                                                  } \
                                                                }; \
                                                              };

#endif
