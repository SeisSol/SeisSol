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
 * Typelists. See Andrei Alexandrescu -- Modern C++ Design.
 **/

#ifndef INITIALIZER_TREE_TYPELIST_HPP_
#define INITIALIZER_TREE_TYPELIST_HPP_ 

struct null_type {};

template<typename T, typename U>
struct typelist {
  typedef T type;
  typedef U next_type;
};

template<typename T, int INDEX>
struct get_type {
  typedef typename get_type<typename T::next_type, INDEX-1>::type type;
};

template<typename T>
struct get_type<T, 0> {
  typedef typename T::type type;
};

template<
  typename T1 = null_type, typename T2 = null_type, typename T3 = null_type,
  typename T4 = null_type, typename T5 = null_type, typename T6 = null_type,
  typename T7 = null_type, typename T8 = null_type, typename T9 = null_type,
  typename T10 = null_type, typename T11 = null_type, typename T12 = null_type,
  typename T13 = null_type, typename T14 = null_type, typename T15 = null_type
>
class make_typelist {
private:
  typedef typename make_typelist<T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15>::result tail;
public:
  typedef typelist<T1, tail> result;
};
template<>
class make_typelist<> {
public:
  typedef null_type result;
};

#endif
