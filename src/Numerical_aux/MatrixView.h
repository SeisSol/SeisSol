/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 
#ifndef MATRIXVIEW_H_
#define MATRIXVIEW_H_

#include <Kernels/precision.hpp>
#include <cassert>
#include <cstring>

template<int Stride>
int colMjrIndex(unsigned row, unsigned col)
{
  return row + col * Stride;
}

class MatrixView {
public:
  explicit MatrixView(real* data, unsigned reals, int (*index)(unsigned, unsigned))
    : data(data),
      reals(reals),
      index(index) {}
  
  void setZero() {
    memset(data, 0, reals * sizeof(real));
  }
  
  inline real& operator()(unsigned row, unsigned col) {
    int idx = (*index)(row, col);
    assert(idx != -1);
    return data[idx];
  }
  
  real* data;
  unsigned reals;
  int (*index)(unsigned, unsigned);
};

/**
 * DenseMatrixView allows hassle-free handling of column-major matrix data
 * arrays. The compiler should optimise most of the operations here,
 * however I would not trust it too much and I would not recommend
 * using this class in highly performance-critical parts of the code.
 */
template<unsigned M, unsigned N>
class DenseMatrixView {
public:
  explicit DenseMatrixView(real* data, unsigned stride = M)
    : data(data),
      stride(stride) {}
  
  void setZero() {
    for (unsigned j = 0; j < N; ++j) {
      for (unsigned i = 0; i < M; ++i) {
        operator()(i, j) = 0.0;
      }
    }
  }
  
  unsigned rows() {
    return M;
  }
  
  unsigned cols() {
    return N;
  }

  real& operator()(unsigned i, unsigned j) {
    return data[j * stride + i];
  }
  
  template<unsigned Mb, unsigned Nb>
  DenseMatrixView<Mb, Nb> block(unsigned originI, unsigned originJ) {
    assert(originI + Mb <= M && originJ + Nb <= N);

    return DenseMatrixView<Mb, Nb>(&data[originJ * stride + originI], stride);
  }
  
  real* data;

private:
  unsigned const stride;
};

#endif
