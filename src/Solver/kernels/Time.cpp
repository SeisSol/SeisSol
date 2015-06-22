/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Time kernel of SeisSol.
 **/

#include "Time.h"

#include <matrix_kernels/sparse.h>
#include <matrix_kernels/dense.h>

#ifndef NDEBUG
#pragma message "compiling time kernel with assertions"
extern long long libxsmm_num_total_flops;
#endif

#include <cstring>
#include <cassert>
#include <stdint.h>
#if defined(__SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

seissol::kernels::Time::Time() {
  // compute the aligned number of basis functions and offsets of the derivatives
  m_derivativesOffsets[0] = 0;
  for( int l_order = 0; l_order < CONVERGENCE_ORDER; l_order++ ) {
    m_numberOfAlignedBasisFunctions[l_order] = getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order, ALIGNMENT );

    if( l_order > 0 ) {
      m_derivativesOffsets[l_order]  =  m_numberOfAlignedBasisFunctions[l_order-1] * NUMBER_OF_QUANTITIES;
      m_derivativesOffsets[l_order] +=  m_derivativesOffsets[l_order-1];
    }
  }

  // intialize the function pointers to the matrix kernels
#define TIME_KERNEL
#include <initialization/bind.h>
#undef TIME_KERNEL
}

void seissol::kernels::Time::computeAder(       double i_timeStepWidth,
                                                real** i_stiffnessMatrices,
                                          const real*  i_degreesOfFreedom,
                                                real   i_starMatrices[3][STAR_NNZ],
                                                real*  o_timeIntegrated,
                                                real*  o_timeDerivatives ) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_degreesOfFreedom)     % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[0]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[1]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[2]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated )      % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeDerivatives)      % ALIGNMENT == 0 || o_timeDerivatives == NULL );

  /*
   * compute ADER scheme.
   */
  // scalars in the taylor-series expansion
  real l_scalar = i_timeStepWidth;

  // temporary result
  real l_temporaryResult[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));
  real l_derivativesBuffer[NUMBER_OF_ALIGNED_DERS] __attribute__((aligned(PAGESIZE_STACK)));

  // initialize time integrated DOFs and derivatives
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof++ ) {
    l_derivativesBuffer[l_dof] = i_degreesOfFreedom[l_dof];
    o_timeIntegrated[l_dof]  = i_degreesOfFreedom[l_dof] * l_scalar;
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
  libxsmm_num_total_flops += NUMBER_OF_ALIGNED_DOFS;
#endif

  for( unsigned int l_dof = NUMBER_OF_ALIGNED_DOFS; l_dof < NUMBER_OF_ALIGNED_DERS; l_dof++ ) {
    l_derivativesBuffer[l_dof] = 0.0;
  }

  // stream out frist derivative (order 0)
  if ( o_timeDerivatives != NULL ) {
    streamstoreFirstDerivative( i_degreesOfFreedom,
                                o_timeDerivatives );
  }

  // compute all derivatives and contributions to the time integrated DOFs
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    // iterate over dimensions 
    for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
      // compute $K_{\xi_c}.Q_k$ and $(K_{\xi_c}.Q_k).A*$
      m_matrixKernels[ (l_derivative-1)*4 + l_c ] ( i_stiffnessMatrices[l_c], l_derivativesBuffer+m_derivativesOffsets[l_derivative-1],  l_temporaryResult,
                                                    NULL,                     NULL,                                                    NULL                                                  ); // These will be be ignored

      m_matrixKernels[ (l_derivative-1)*4 + 3   ] ( l_temporaryResult,        i_starMatrices[l_c],                                     l_derivativesBuffer+m_derivativesOffsets[l_derivative],
                                                    NULL,                     NULL,                                                    NULL                                                  ); // These will be be ignored
    }

    // update scalar for this derivative
    l_scalar *= i_timeStepWidth / real(l_derivative+1);

    // update time integrated DOFs
    integrateInTime( l_derivativesBuffer,
                     l_scalar,
                     l_derivative,
                     o_timeIntegrated,
                     o_timeDerivatives );
  }
}

void seissol::kernels::Time::streamstoreFirstDerivative( const real* i_degreesOfFreedom,
                                                                real* o_derivativesBuffer ) {
#if defined(__AVX512F__) 
#if defined(DOUBLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=8 ) {
    _mm512_stream_pd(&(o_derivativesBuffer[l_dof]), _mm512_load_pd(&(i_degreesOfFreedom[l_dof])));
  }
#elif defined(SINGLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=16 ) {
    _mm512_stream_ps(&(o_derivativesBuffer[l_dof]), _mm512_load_ps(&(i_degreesOfFreedom[l_dof])));
  }
#else
#error no precision was defined 
#endif
#elif defined(__AVX__)
#if defined(DOUBLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=4 ) {
    _mm256_stream_pd(&(o_derivativesBuffer[l_dof]), _mm256_load_pd(&(i_degreesOfFreedom[l_dof])));
  }
#elif defined(SINGLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=8 ) {
    _mm256_stream_ps(&(o_derivativesBuffer[l_dof]), _mm256_load_ps(&(i_degreesOfFreedom[l_dof])));
  }
#else
#error no precision was defined 
#endif
#elif defined(__SSE3__)
#if defined(DOUBLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=2 ) {
    _mm_stream_pd(&(o_derivativesBuffer[l_dof]), _mm_load_pd(&(i_degreesOfFreedom[l_dof])));
  }
#elif defined(SINGLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=4 ) {
    _mm_stream_ps(&(o_derivativesBuffer[l_dof]), _mm_load_ps(&(i_degreesOfFreedom[l_dof])));
  }
#else
#error no precision was defined 
#endif
#elif defined(__MIC__)
#if defined(DOUBLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=8 ) {
    _mm512_storenrngo_pd(&(o_derivativesBuffer[l_dof]), _mm512_load_pd(&(i_degreesOfFreedom[l_dof])));
  }
#elif defined(SINGLE_PRECISION)
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof+=16 ) {
    _mm512_storenrngo_ps(&(o_derivativesBuffer[l_dof]), _mm512_load_ps(&(i_degreesOfFreedom[l_dof])));
  }
#else
#error no precision was defined 
#endif
#else
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof++ ) {
    o_derivativesBuffer[l_dof] = i_degreesOfFreedom[l_dof];
  }
#endif
}

void seissol::kernels::Time::integrateInTime( const real*        i_derivativesBuffer,
                                                    real         i_scalar,
                                                    unsigned int i_derivative,
                                                    real*        o_timeIntegrated,
                                                    real*        o_timeDerivatives ) {
#if defined(DOUBLE_PRECISION)
#if defined(__AVX512F__)
  __m512d l_intrin_scalar = _mm512_broadcast_f64x4(_mm256_broadcast_sd(&i_scalar));
#elif defined(__AVX__)
  __m256d l_intrin_scalar = _mm256_broadcast_sd(&i_scalar);
#elif defined(__SSE3__)
  __m128d l_intrin_scalar = _mm_loaddup_pd(&i_scalar);
#elif defined(__MIC__)
  __m512d l_intrin_scalar = _mm512_extload_pd(&i_scalar, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#else
  real l_scalar = i_scalar;
#endif
#elif defined(SINGLE_PRECISION)
#if defined(__AVX512F__)
  __m512 l_intrin_scalar = _mm512_broadcast_f32x8(_mm256_broadcast_ss(&i_scalar));
#elif defined(__AVX__)
  __m256 l_intrin_scalar = _mm256_broadcast_ss(&i_scalar);
#elif defined(__SSE3__)
  __m128 l_intrin_scalar = _mm_load_ss(&i_scalar);
  l_intrin_scalar = _mm_shuffle_ps(l_intrin_scalar, l_intrin_scalar, 0x00);
#elif defined(__MIC__)
  __m512 l_intrin_scalar = _mm512_extload_ps(&i_scalar, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#else
  real l_scalar = i_scalar;
#endif
#else
#error no precision was defined 
#endif

#ifndef NDEBUG
  long long temp_flops = NUMBER_OF_QUANTITIES*(m_numberOfAlignedBasisFunctions[i_derivative])*2.0;
#ifdef _OPENMP
#pragma omp atomic
#endif
  libxsmm_num_total_flops += temp_flops;
#endif

  if (o_timeDerivatives == NULL) {
    for( unsigned int l_quantity = 0; l_quantity < NUMBER_OF_QUANTITIES; l_quantity++ ) {
      unsigned int l_tint_offset = l_quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS;
      unsigned int l_ders_offset = m_derivativesOffsets[i_derivative] + (l_quantity*m_numberOfAlignedBasisFunctions[i_derivative]);
#if defined(__AVX512F__) || defined(__MIC__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=8 ) {
        __m512d l_source =  _mm512_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm512_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm512_fmadd_pd(l_intrin_scalar, l_source, _mm512_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=16 ) {
        __m512 l_source =  _mm512_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm512_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm512_fmadd_ps(l_intrin_scalar, l_source, _mm512_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#else
#error no precision was defined 
#endif
#elif defined(__AVX2__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=4 ) {
        __m256d l_source = _mm256_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_fmadd_pd(l_intrin_scalar, l_source, _mm256_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=8 ) {
        __m256 l_source = _mm256_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_fmadd_ps(l_intrin_scalar, l_source, _mm256_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#else
#error no precision was defined 
#endif
#elif defined(__AVX__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=4 ) {
        __m256d l_source = _mm256_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_add_pd(_mm256_mul_pd(l_intrin_scalar, l_source), _mm256_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=8 ) {
        __m256 l_source = _mm256_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_add_ps(_mm256_mul_ps(l_intrin_scalar, l_source), _mm256_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#else
#error no precision was defined 
#endif
#elif defined(__SSE3__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=2 ) {
        __m128d l_source = _mm_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm_add_pd(_mm_mul_pd(l_intrin_scalar, l_source), _mm_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=4 ) {
        __m128 l_source = _mm_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm_add_ps(_mm_mul_ps(l_intrin_scalar, l_source), _mm_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
      }
#else
#error no precision was defined 
#endif
#else
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction++ ) {
        o_timeIntegrated[ l_tint_offset + l_basisFunction ] += l_scalar * i_derivativesBuffer[ l_ders_offset + l_basisFunction  ];
      }
#endif
    }
  } else {
    for( unsigned int l_quantity = 0; l_quantity < NUMBER_OF_QUANTITIES; l_quantity++ ) {
      unsigned int l_tint_offset = l_quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS;
      unsigned int l_ders_offset = m_derivativesOffsets[i_derivative] + (l_quantity*m_numberOfAlignedBasisFunctions[i_derivative]);
#if defined(__AVX512F__) || defined(__MIC__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=8 ) {
        __m512d l_source =  _mm512_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm512_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm512_fmadd_pd(l_intrin_scalar, l_source, _mm512_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
#ifdef __AVX512F__
        _mm512_stream_pd(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
#endif
#ifdef __MIC__
        _mm512_storenrngo_pd(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
#endif
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=16 ) {
        __m512 l_source =  _mm512_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm512_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm512_fmadd_ps(l_intrin_scalar, l_source, _mm512_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
#ifdef __AVX512F__
        _mm512_stream_ps(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
#endif
#ifdef __MIC__
        _mm512_storenrngo_ps(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
#endif
      }
#else
#error no precision was defined 
#endif
#elif defined(__AVX2__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=4 ) {
        __m256d l_source = _mm256_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_fmadd_pd(l_intrin_scalar, l_source, _mm256_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
        _mm256_stream_pd(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=8 ) {
        __m256 l_source = _mm256_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_fmadd_ps(l_intrin_scalar, l_source, _mm256_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
        _mm256_stream_ps(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
      }
#else
#error no precision was defined 
#endif
#elif defined(__AVX__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=4 ) {
        __m256d l_source = _mm256_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_add_pd(_mm256_mul_pd(l_intrin_scalar, l_source), _mm256_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
        _mm256_stream_pd(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=8 ) {
        __m256 l_source = _mm256_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm256_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                        _mm256_add_ps(_mm256_mul_ps(l_intrin_scalar, l_source), _mm256_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
        _mm256_stream_ps(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
      }
#else
#error no precision was defined 
#endif
#elif defined(__SSE3__)
#if defined(DOUBLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=2 ) {
        __m128d l_source = _mm_load_pd(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm_store_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                       _mm_add_pd(_mm_mul_pd(l_intrin_scalar, l_source), _mm_load_pd(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
        _mm_stream_pd(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
      }
#elif defined(SINGLE_PRECISION)
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction+=4 ) {
        __m128 l_source = _mm_load_ps(&(i_derivativesBuffer[ l_ders_offset + l_basisFunction  ]));
        _mm_store_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]),
                       _mm_add_ps(_mm_mul_ps(l_intrin_scalar, l_source), _mm_load_ps(&(o_timeIntegrated[ l_tint_offset + l_basisFunction ]))));
        _mm_stream_ps(&(o_timeDerivatives[ l_ders_offset + l_basisFunction ]), l_source);
      }
#else
#error no precision was defined 
#endif
#else
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[i_derivative]; l_basisFunction++ ) {
        o_timeIntegrated[ l_tint_offset + l_basisFunction ] += l_scalar * i_derivativesBuffer[ l_ders_offset + l_basisFunction  ];
        o_timeDerivatives[ l_ders_offset + l_basisFunction ] = i_derivativesBuffer[ l_ders_offset + l_basisFunction  ];
      }
#endif
    }
  }
}

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  // initialization
  o_nonZeroFlops  += NUMBER_OF_DOFS;
  o_hardwareFlops += NUMBER_OF_ALIGNED_DOFS;

  // interate over derivatives
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    // iterate over dimensions
    for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
      o_nonZeroFlops  += m_nonZeroFlops[  (l_derivative-1)*4 + l_c ];
      o_hardwareFlops += m_hardwareFlops[ (l_derivative-1)*4 + l_c ];

      o_nonZeroFlops  += m_nonZeroFlops[  (l_derivative-1)*4 + 3   ];
      o_hardwareFlops += m_hardwareFlops[ (l_derivative-1)*4 + 3   ];
    }

    // update of time integrated DOFs
    o_nonZeroFlops  += seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER - l_derivative ) * NUMBER_OF_QUANTITIES * 2;
    o_hardwareFlops += seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER - l_derivative ) * NUMBER_OF_QUANTITIES * 2;
  }

}

void seissol::kernels::Time::computeExtrapolation(       double i_expansionPoint,
                                                         double i_evaluationPoint,
                                                   const real*  i_timeDerivatives,
                                                         real*  o_timeEvaluated ) {
  /*
   * assert valid input.
   */
  // assert alignments
  assert( ((uintptr_t)i_timeDerivatives) % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeEvaluated)   % ALIGNMENT == 0 );

  // assert that non-backward evaluation in time
  assert( i_evaluationPoint >= i_expansionPoint );

  /*
   * compute extrapolation.
   */
  // copy DOFs (scalar==1)
  memcpy( o_timeEvaluated, i_timeDerivatives, NUMBER_OF_ALIGNED_DOFS*sizeof(real) );

  // initialize scalars in the taylor series expansion (1st derivative)
  real l_deltaT = i_evaluationPoint - i_expansionPoint;
  real l_scalar = 1.0;
 
  // evaluate taylor series expansion at evvaluation point
  for( unsigned int l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    l_scalar *= l_deltaT;
    l_scalar /= (real) l_derivative;

    integrateInTime( i_timeDerivatives,
                     l_scalar,
                     l_derivative,
                     o_timeEvaluated,
                     NULL );
  }
}

void seissol::kernels::Time::computeIntegral(       double i_expansionPoint,
                                                    double i_integrationStart,
                                                    double i_integrationEnd,
                                              const real*  i_timeDerivatives,
                                                    real   o_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] ) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // reset time integrated degrees of freedom
  memset( o_timeIntegrated, 0, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES*sizeof(real) );

  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = (real) 1;
  real l_secondTerm = (real) 1;
  real l_factorial  = (real) 1;
  real l_scalar;
 
  // iterate over time derivatives
  for(int l_derivative = 0; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (real)(l_derivative+1);

    l_scalar  = l_firstTerm - l_secondTerm;
    l_scalar /= l_factorial;

    integrateInTime( i_timeDerivatives,
                     l_scalar,
                     l_derivative,
                     o_timeIntegrated,
                     NULL );
  }
}

void seissol::kernels::Time::computeIntegrals( unsigned short      i_ltsSetup,
                                               const enum faceType i_faceTypes[4],
                                               const double        i_currentTime[5],
                                               double              i_timeStepWidth,
                                               real  *const        i_timeDofs[4],
                                               real                o_integrationBuffer[4][NUMBER_OF_ALIGNED_DOFS],
                                               real  *             o_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] ) {
  /*
   * assert valid input.
   */
  // only lower 10 bits are used for lts encoding
  assert (i_ltsSetup < 2048 );

#ifndef NDEBUG
  // alignment of the time derivatives/integrated dofs and the buffer
  for( int l_dofeighbor = 0; l_dofeighbor < 4; l_dofeighbor++ ) {
    assert( ((uintptr_t)i_timeDofs[l_dofeighbor])          % ALIGNMENT == 0 );
    assert( ((uintptr_t)o_integrationBuffer[l_dofeighbor]) % ALIGNMENT == 0 );
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for( unsigned int l_dofeighbor = 0; l_dofeighbor < 4; l_dofeighbor++ ) {
    // collect information only in the case that neighboring element contributions are required
    if( i_faceTypes[l_dofeighbor] != outflow && i_faceTypes[l_dofeighbor] != dynamicRupture ) {
      // check if the time integration is already done (-> copy pointer)
      if( (i_ltsSetup >> l_dofeighbor ) % 2 == 0 ) {
        o_timeIntegrated[l_dofeighbor] = i_timeDofs[l_dofeighbor];
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        seissol::kernels::Time::computeIntegral(  i_currentTime[    l_dofeighbor+1],
                                                  i_currentTime[    0           ],
                                                  i_currentTime[    0           ] + i_timeStepWidth,
                                                  i_timeDofs[       l_dofeighbor],
                                                  o_integrationBuffer[ l_dofeighbor]                  );

        o_timeIntegrated[l_dofeighbor] = o_integrationBuffer[ l_dofeighbor];
      }
    }
  }
}

void seissol::kernels::Time::computeIntegrals( unsigned short      i_ltsSetup,
                                               const enum faceType i_faceTypes[4],
                                               const double        i_timeStepStart,
                                               const double        i_timeStepWidth,
                                               real * const        i_timeDofs[4],
                                               real                o_integrationBuffer[4][NUMBER_OF_ALIGNED_DOFS],
                                               real *              o_timeIntegrated[4] ) {
  double l_startTimes[5];
  l_startTimes[0] = i_timeStepStart;
  l_startTimes[1] = l_startTimes[2] = l_startTimes[3] = l_startTimes[4] = 0;

  // adjust start times for GTS on derivatives
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    if( (i_ltsSetup >> l_face + 4 ) % 2 ) {
      l_startTimes[l_face+1] = i_timeStepStart;
    }
  }

  // call the more general assembly
  computeIntegrals( i_ltsSetup,
                    i_faceTypes,
                    l_startTimes,
                    i_timeStepWidth,
                    i_timeDofs,
                    o_integrationBuffer,
                    o_timeIntegrated );
}
