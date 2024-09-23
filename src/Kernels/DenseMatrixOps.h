// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPS_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPS_H_

#if defined(__SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#if defined(__AVX512F__)
#include "DenseMatrixOpsAVX512.h"
#elif defined(__MIC__)
#include "DenseMatrixOpsMIC.h"
#elif defined(__AVX2__)
#include "DenseMatrixOpsAVX2.h"
#elif defined(__AVX__)
#include "DenseMatrixOpsAVX.h"
#elif defined(__SSE3__)
#include "DenseMatrixOpsSSE3.h"
#else
#include "DenseMatrixOpsNoarch.h"
#endif

#include <cassert>

extern long long libxsmm_num_total_flops;

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
      assert(numberOfReals % DMO_INCREMENT == 0);
      
      for (unsigned i = 0; i < numberOfReals; i += DMO_INCREMENT) {
        DMO_STREAM(&X[i], &Y[i])
      }
    }
    
    /**
     * Computes Y = scalar * X.
     * The number of rows must be a multiple of the alignment.
     * 
     * @param scalar The value that X shall be scaled with.
     * @param rows Number of rows of X and Y.
     * @param columns Number of columns of X and Y.
     * @param X
     * @param LDX Size of leading dimension of X.
     * @param Y
     * @param LDY Size of leading dimension of Y.
     **/
    inline void SXt(  real scalar,
                      unsigned rows,
                      unsigned columns,
                      real const* X,
                      unsigned LDX,
                      real* Y,
                      unsigned LDY )
    {
      assert(rows % DMO_INCREMENT == 0);
      
      DMO_BROADCAST(&scalar, intrin_scalar)
      
      for (unsigned int column = 0; column < columns; ++column) {
        real const* xCol = &X[column * LDX];
        real* yCol = &Y[column * LDY];
        for (unsigned row = 0; row < rows; row += DMO_INCREMENT) {
          DMO_SXT(intrin_scalar, xCol+row, yCol+row)
        }
      }
            
#ifndef NDEBUG
      long long temp_flops = rows * columns;
#ifdef _OPENMP
#pragma omp atomic
#endif
      libxsmm_num_total_flops += temp_flops;
#endif
    }
    
    /**
     * Computes Y = scalar * X + Y.
     * The number of rows must be a multiple of the alignment.
     * 
     * @param scalar The value that X shall be scaled with.
     * @param rows Number of rows of X and Y.
     * @param columns Number of columns of X and Y.
     * @param X
     * @param LDX Size of leading dimension of X.
     * @param Y
     * @param LDY Size of leading dimension of Y.
     * @param XStore If XStore != NULL a copy of X will be placed at XStore.
     **/
    inline void SXtYp(  real scalar,
                        unsigned rows,
                        unsigned columns,
                        real const* X,
                        unsigned LDX,
                        real* Y,
                        unsigned LDY,
                        real* XStore = NULL )
    {
      assert(rows % DMO_INCREMENT == 0);
      
      DMO_BROADCAST(&scalar, intrin_scalar)
      
      if (XStore == NULL) {
        for (unsigned int column = 0; column < columns; ++column) {
          real const* xCol = &X[column * LDX];
          real* yCol = &Y[column * LDY];
          for (unsigned row = 0; row < rows; row += DMO_INCREMENT) {
            DMO_SXTYP(intrin_scalar, xCol+row, yCol+row)
          }
        }
      } else {
        for (unsigned int column = 0; column < columns; ++column) {
          real const* xCol = &X[column * LDX];
          real* yCol = &Y[column * LDY];
          real* xColStore = &XStore[column * LDX];
          for (unsigned row = 0; row < rows; row += DMO_INCREMENT) {
            DMO_SXTYP(intrin_scalar, xCol+row, yCol+row)
            DMO_STREAM(xCol+row, xColStore+row);
          }
        }
      }
            
#ifndef NDEBUG
      long long temp_flops = rows * columns * 2;
#ifdef _OPENMP
#pragma omp atomic
#endif
      libxsmm_num_total_flops += temp_flops;
#endif
    }
    
    /**
     * Computes Z = (X-Y) * scalar.
     * The number of rows must be a multiple of the alignment.
     * 
     * @param scalar The value that X-Y shall be scaled with.
     * @param rows Number of rows of X, Y and Z.
     * @param columns Number of columns X, Y and Z.
     * @param X
     * @param LDX Size of leading dimension of X.
     * @param Y
     * @param LDY Size of leading dimension of Y.
     * @param Z
     * @param LDZ Size of leading dimension of Z.
     **/
    inline void XYmSt(  real scalar,
                        unsigned rows,
                        unsigned columns,
                        real const* X,
                        unsigned LDX,
                        real const* Y,
                        unsigned LDY,
                        real* Z,
                        unsigned LDZ )
  {
      assert(rows % DMO_INCREMENT == 0);
      
      DMO_BROADCAST(&scalar, intrin_scalar)
      
      for (unsigned int column = 0; column < columns; ++column) {
        real const* xCol = &X[column * LDX];
        real const* yCol = &Y[column * LDY];
        real* zCol = &Z[column * LDZ];
        for (unsigned row = 0; row < rows; row += DMO_INCREMENT) {
          DMO_XYMST(intrin_scalar, xCol+row, yCol+row, zCol+row)
        }
      }

#ifndef NDEBUG
      long long temp_flops = rows * columns * 2;
#ifdef _OPENMP
#pragma omp atomic
#endif
      libxsmm_num_total_flops += temp_flops;
#endif
    }
    
    /**
     * Computes Z = (X-Y) * scalar + Z.
     * The number of rows must be a multiple of the alignment.
     * 
     * @param scalar The value that X-Y shall be scaled with.
     * @param rows Number of rows of X, Y and Z.
     * @param columns Number of columns X, Y and Z.
     * @param X
     * @param LDX Size of leading dimension of X.
     * @param Y
     * @param LDY Size of leading dimension of Y.
     * @param Z
     * @param LDZ Size of leading dimension of Z.
     **/
    inline void XYmStZp(  real scalar,
                          unsigned rows,
                          unsigned columns,
                          real const* X,
                          unsigned LDX,
                          real const* Y,
                          unsigned LDY,
                          real* Z,
                          unsigned LDZ )
    {
      assert(rows % DMO_INCREMENT == 0);
      
      DMO_BROADCAST(&scalar, intrin_scalar)
      
      for (unsigned int column = 0; column < columns; ++column) {
        real const* xCol = &X[column * LDX];
        real const* yCol = &Y[column * LDY];
        real* zCol = &Z[column * LDZ];
        for (unsigned row = 0; row < rows; row += DMO_INCREMENT) {
          DMO_XYMSTZP(intrin_scalar, xCol+row, yCol+row, zCol+row)
        }
      }

#ifndef NDEBUG
      long long temp_flops = rows * columns * 3;
#ifdef _OPENMP
#pragma omp atomic
#endif
      libxsmm_num_total_flops += temp_flops;
#endif
    }
  }
}


#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPS_H_

