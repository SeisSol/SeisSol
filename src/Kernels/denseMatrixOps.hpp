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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Common kernel-level functions
 **/

#ifndef KERNElS_DENSEMATRIXOPS_HPP_
#define KERNElS_DENSEMATRIXOPS_HPP_

#include <Kernels/precision.hpp>
#if defined(__SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cassert>

namespace seissol {
  namespace kernels {
    /** Stores X in Y with non-temporal hint.
     * 
     * @param numberOfReals The size of X and Y.
     * @param X
     * @param Y
     */
    inline void streamstore( unsigned numberOfReals,
                             real const* X,
                             real* Y )
    {
      assert((numberOfReals * sizeof(real)) % ALIGNMENT == 0);

#if defined(__AVX512F__) 
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=8 ) {
        _mm512_stream_pd(&(Y[l_dof]), _mm512_load_pd(&(X[l_dof])));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=16 ) {
        _mm512_stream_ps(&(Y[l_dof]), _mm512_load_ps(&(X[l_dof])));
      }
#else
#error no precision was defined 
#endif
#elif defined(__AVX__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=4 ) {
        _mm256_stream_pd(&(Y[l_dof]), _mm256_load_pd(&(X[l_dof])));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=8 ) {
        _mm256_stream_ps(&(Y[l_dof]), _mm256_load_ps(&(X[l_dof])));
      }
#else
#error no precision was defined 
#endif
#elif defined(__SSE3__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=2 ) {
        _mm_stream_pd(&(Y[l_dof]), _mm_load_pd(&(X[l_dof])));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=4 ) {
        _mm_stream_ps(&(Y[l_dof]), _mm_load_ps(&(X[l_dof])));
      }
#else
#error no precision was defined 
#endif
#elif defined(__MIC__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=8 ) {
        _mm512_storenrngo_pd(&(Y[l_dof]), _mm512_load_pd(&(X[l_dof])));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof+=16 ) {
        _mm512_storenrngo_ps(&(Y[l_dof]), _mm512_load_ps(&(X[l_dof])));
      }
#else
#error no precision was defined 
#endif
#else
      for( unsigned int l_dof = 0; l_dof < numberOfReals; l_dof++ ) {
        Y[l_dof] = X[l_dof];
      }
#endif
    }
    /**
     * Scales a (sub-)matrix and adds it to another matrix. Computes Y += scalar * X.
     * The number of rows must be a multiple of the alignment.
     * 
     * @param scalar The value that X shall be scaled with.
     * @param columns Number of columns of X and Y
     * @param X
     * @param LDX Size of leading dimension of X.
     * @param Y
     * @param LDY Size of leading dimension of Y.
     * @param XStore If XStore != NULL a copy of X will be placed at XStore.
     **/
    inline void aligned_axpy( real scalar,
                              unsigned rows,
                              unsigned columns,
                              real const* X,
                              unsigned LDX,
                              real* Y,
                              unsigned LDY,
                              real* XStore = NULL )
    {
      assert((rows * sizeof(real)) % ALIGNMENT == 0);

#if defined(DOUBLE_PRECISION)
#if defined(__AVX512F__)
      __m512d l_intrin_scalar = _mm512_broadcastsd_pd(_mm_load_sd(&scalar));
#elif defined(__AVX__)
      __m256d l_intrin_scalar = _mm256_broadcast_sd(&scalar);
#elif defined(__SSE3__)
      __m128d l_intrin_scalar = _mm_loaddup_pd(&scalar);
#elif defined(__MIC__)
      __m512d l_intrin_scalar = _mm512_extload_pd(&scalar, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#endif
#elif defined(SINGLE_PRECISION)
#if defined(__AVX512F__)
      __m512 l_intrin_scalar = _mm512_broadcastss_ps(_mm_load_ss(&scalar));
#elif defined(__AVX__)
      __m256 l_intrin_scalar = _mm256_broadcast_ss(&scalar);
#elif defined(__SSE3__)
      __m128 l_intrin_scalar = _mm_load_ss(&scalar);
      l_intrin_scalar = _mm_shuffle_ps(l_intrin_scalar, l_intrin_scalar, 0x00);
#elif defined(__MIC__)
      __m512 l_intrin_scalar = _mm512_extload_ps(&scalar, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#endif
#else
#error no precision was defined 
#endif

#ifndef NDEBUG
      long long temp_flops = columns*rows*2.0;
#ifdef _OPENMP
#pragma omp atomic
#endif
      libxsmm_num_total_flops += temp_flops;
#endif

      if (XStore == NULL) {
        for (unsigned int column = 0; column < columns; ++column) {
          unsigned yOffset = column * LDY;
          unsigned xOffset = column * LDX;
#if defined(__AVX512F__) || defined(__MIC__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 8) {
            __m512d x = _mm512_load_pd(&X[xOffset + row]);
            _mm512_store_pd(&Y[yOffset + row], _mm512_fmadd_pd(l_intrin_scalar, x, _mm512_load_pd(&Y[yOffset + row])));
          }
#elif defined(SINGLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 16) {
            __m512 x = _mm512_load_ps(&X[xOffset + row]);
            _mm512_store_ps(&Y[yOffset + row], _mm512_fmadd_ps(l_intrin_scalar, x, _mm512_load_ps(&Y[yOffset + row])));
          }
#else
#error no precision was defined 
#endif
#elif defined(__AVX2__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 4) {
            __m256d x = _mm256_load_pd(&X[xOffset + row]);
            _mm256_store_pd(&Y[yOffset + row], _mm256_fmadd_pd(l_intrin_scalar, x, _mm256_load_pd(&Y[yOffset + row])));
          }
#elif defined(SINGLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 8) {
            __m256 x = _mm256_load_ps(&X[xOffset + row]);
            _mm256_store_ps(&Y[yOffset + row], _mm256_fmadd_ps(l_intrin_scalar, x, _mm256_load_ps(&Y[yOffset + row])));
          }
#else
#error no precision was defined 
#endif
#elif defined(__AVX__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 4) {
            __m256d x = _mm256_load_pd(&X[xOffset + row]);
            _mm256_store_pd(&Y[yOffset + row], _mm256_add_pd(_mm256_mul_pd(l_intrin_scalar, x), _mm256_load_pd(&Y[yOffset + row])));
          }
#elif defined(SINGLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 8) {
            __m256 x = _mm256_load_ps(&X[xOffset + row]);
            _mm256_store_ps(&Y[yOffset + row], _mm256_add_ps(_mm256_mul_ps(l_intrin_scalar, x), _mm256_load_ps(&Y[yOffset + row])));
          }
#else
#error no precision was defined 
#endif
#elif defined(__SSE3__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 2) {
            __m128d x = _mm_load_pd(&X[xOffset + row]);
            _mm_store_pd(&Y[yOffset + row], _mm_add_pd(_mm_mul_pd(l_intrin_scalar, x), _mm_load_pd(&Y[yOffset + row])));
          }
#elif defined(SINGLE_PRECISION)    
          for (unsigned row = 0; row < rows; row += 4) {
            __m128 x = _mm_load_ps(&X[xOffset + row]);
            _mm_store_ps(&Y[yOffset + row], _mm_add_ps(_mm_mul_pd(l_intrin_scalar, x), _mm_load_ps(&Y[yOffset + row])));
          }
#else
#error no precision was defined 
#endif
#else
          for (unsigned row = 0; row < rows; ++row) {
            Y[yOffset + row] += scalar * X[xOffset + row];
          }
#endif
        }
      } else {
        for (unsigned int column = 0; column < columns; ++column) {
          unsigned yOffset = column * LDY;
          unsigned xOffset = column * LDX;
#if defined(__AVX512F__) || defined(__MIC__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 8) {
            __m512d x = _mm512_load_pd(&X[xOffset + row]);
            _mm512_store_pd(&Y[yOffset + row], _mm512_fmadd_pd(l_intrin_scalar, x, _mm512_load_pd(&Y[yOffset + row])));
#ifdef __AVX512F__
            _mm512_stream_pd(&XStore[xOffset + row], x);
#endif
#ifdef __MIC__
            _mm512_storenrngo_pd(&XStore[xOffset + row], x);
#endif
          }
#elif defined(SINGLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 16) {
            __m512 x = _mm512_load_ps(&X[xOffset + row]);
            _mm512_store_ps(&Y[yOffset + row], _mm512_fmadd_ps(l_intrin_scalar, x, _mm512_load_ps(&Y[yOffset + row])));
#ifdef __AVX512F__
            _mm512_stream_ps(&XStore[xOffset + row], x);
#endif
#ifdef __MIC__
            _mm512_storenrngo_ps(&XStore[xOffset + row], x);
#endif
          }
#else
#error no precision was defined 
#endif
#elif defined(__AVX2__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 4) {
            __m256d x = _mm256_load_pd(&X[xOffset + row]);
            _mm256_store_pd(&Y[yOffset + row], _mm256_fmadd_pd(l_intrin_scalar, x, _mm256_load_pd(&Y[yOffset + row])));
            _mm256_stream_pd(&XStore[xOffset + row], x);
          }
#elif defined(SINGLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 8) {
            __m256 x = _mm256_load_ps(&X[xOffset + row]);
            _mm256_store_ps(&Y[yOffset + row], _mm256_fmadd_ps(l_intrin_scalar, x, _mm256_load_ps(&Y[yOffset + row])));
            _mm256_stream_ps(&XStore[xOffset + row], x);
          }
#else
#error no precision was defined 
#endif
#elif defined(__AVX__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 4) {
            __m256d x = _mm256_load_pd(&X[xOffset + row]);
            _mm256_store_pd(&Y[yOffset + row], _mm256_add_pd(_mm256_mul_pd(l_intrin_scalar, x), _mm256_load_pd(&Y[yOffset + row])));
            _mm256_stream_pd(&XStore[xOffset + row], x);
          }
#elif defined(SINGLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 8) {
            __m256 x = _mm256_load_ps(&X[xOffset + row]);
            _mm256_store_ps(&Y[yOffset + row], _mm256_add_ps(_mm256_mul_ps(l_intrin_scalar, x), _mm256_load_ps(&Y[yOffset + row])));
            _mm256_stream_ps(&XStore[xOffset + row], x);
          }
#else
#error no precision was defined 
#endif
#elif defined(__SSE3__)
#if defined(DOUBLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 2) {
            __m128d x = _mm_load_pd(&X[xOffset + row]);
            _mm_store_pd(&Y[yOffset + row], _mm_add_pd(_mm_mul_pd(l_intrin_scalar, x), _mm_load_pd(&Y[yOffset + row])));
            _mm_stream_pd(&XStore[xOffset + row], x);
          }
#elif defined(SINGLE_PRECISION)
          for (unsigned row = 0; row < rows; row += 4) {
            __m128 x = _mm_load_ps(&X[xOffset + row]);
            _mm_store_ps(&Y[yOffset + row], _mm_add_ps(_mm_mul_pd(l_intrin_scalar, x), _mm_load_ps(&Y[yOffset + row])));
            _mm_stream_ps(&XStore[xOffset + row], x);
          }
#else
#error no precision was defined 
#endif
#else
          for (unsigned row = 0; row < rows; ++row) {
            Y[yOffset + row] += scalar * X[xOffset + row];
            XStore[xOffset + row] = X[xOffset + row];
          }
#endif
        }
      }
    }
  }
}

#endif
