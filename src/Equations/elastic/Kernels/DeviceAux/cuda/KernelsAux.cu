#include <Kernels/precision.hpp>
#include <generated_code/init.h>
#include <generated_code/tensor.h>
#include <vector>
#include <yateto.h>
#include <cuda.h>
#include <cstdio>

#include <iostream>

void checkErr(const std::string &file, int line) {
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    std::cout << '\n'
           << file << ", line " << line << ": " << cudaGetErrorString(error) << " (" << error
           << ")\n" << std::endl;
  }
}

#define CHECK_ERR checkErr(__FILE__, __LINE__)

#include "Equations/elastic/Kernels/DeviceAux/KernelsAux.h"

namespace seissol::kernels::time::aux {
constexpr size_t Quantities = 9;

template<size_t Order> struct FunctionCount;

template<> struct FunctionCount<6> { constexpr static size_t Count = 56; };
template<> struct FunctionCount<5> { constexpr static size_t Count = 35; };
template<> struct FunctionCount<4> { constexpr static size_t Count = 20; };
template<> struct FunctionCount<3> { constexpr static size_t Count = 10; };
template<> struct FunctionCount<2> { constexpr static size_t Count = 4; };
template<> struct FunctionCount<1> { constexpr static size_t Count = 1; };

template<size_t Order>
constexpr size_t Functions = FunctionCount<Order>::Count;

__constant__ real derivative6X[Functions<6> * Functions<5>];
__constant__ real derivative6Y[Functions<6> * Functions<5>];
__constant__ real derivative6Z[Functions<6> * Functions<5>];

constexpr size_t AderMultiple = 1;

#define INDEX(n) (threadIdx.x + Blocksize * (n))
#define OFFSET(n) ((n) * (AderMultiple * blockIdx.x + threadIdx.y))

#define INDEX2(n) (threadIdx.x*2 + Blocksize * (n))

#define INDEX4(n) (threadIdx.x*4 + Blocksize * (n))

#define MATRIX2(I,J,N,M) (INDEX2((I) + (J) * (N)) + OFFSET(Blocksize * (N) * (M)))

#define MATRIX4(I,J,N,M) (INDEX4((I) + (J) * (N)) + OFFSET(Blocksize * (N) * (M)))

#define MATRIX(I,J,N,M) (INDEX((I) + (J) * (N)) + OFFSET(Blocksize * (N) * (M)))
#define VECTOR(I,N) (INDEX(I) + OFFSET(Blocksize * (N)))

template<bool>
__device__ __forceinline__ void assign(real& target, real source);

template<>
__device__ __forceinline__ void assign<false>(real& target, real source) {
  target = source;
}

template<>
__device__ __forceinline__ void assign<true>(real& target, real source) {
  target += source;
}

template<bool>
__device__ __forceinline__ void assign(float4& target, float4 source);

template<>
__device__ __forceinline__ void assign<false>(float4& target, float4 source) {
  target = source;
}

template<>
__device__ __forceinline__ void assign<true>(float4& target, float4 source) {
  float4 target2;
  target2.x = target.x + source.x;
  target2.y = target.y + source.y;
  target2.z = target.z + source.z;
  target2.w = target.w + source.w;
  target = target2;
}

template<size_t TileN1, size_t TileN2, size_t TileK, size_t Order, bool Override>
__device__ __forceinline__ void dgkernel(real* __restrict__ output, const real* __restrict__ dQ, const real* __restrict__ coordinates, const real* __restrict__ derivative,
    real lambda, real mu, real rhoD, int offset) {
    real l2mu = lambda + 2 * mu;
    real n1 = coordinates[VECTOR(offset, 9)];
    real n2 = coordinates[VECTOR(offset+1, 9)];
    real n3 = coordinates[VECTOR(offset+2, 9)];

{
    constexpr size_t TileN = TileN1;
        for (int i = 0; i < Functions<Order>; i += TileN) {
        {
            constexpr size_t TileM = 6;
            constexpr int j = 0;
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                for (int k = 0; k < Functions<Order+1>; k += TileK) {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivative[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<Override>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0] + n2 * result[ii][3] + n3 * result[ii][5]));
                assign<Override>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][3] + n2 * result[ii][1] + n3 * result[ii][4]));
                assign<Override>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][5] + n2 * result[ii][4] + n3 * result[ii][2]));
            }
        }
    }
}
    {
        constexpr size_t TileN = TileN2;
    for (int i = 0; i < Functions<Order>; i += TileN) {
        {
            constexpr size_t TileM = 3;
            constexpr int j = 6;
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                for (int k = 0; k < Functions<Order+1>; k += TileK) {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivative[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<Override>(output[MATRIX(i+ii, 0, Functions<Order>, Quantities)], -n1 * l2mu * result[ii][0] - lambda * (n2 * result[ii][1] + n3 * result[ii][2]));
                assign<Override>(output[MATRIX(i+ii, 1, Functions<Order>, Quantities)], -n2 * l2mu * result[ii][1] - lambda * (n1 * result[ii][0] + n3 * result[ii][2]));
                assign<Override>(output[MATRIX(i+ii, 2, Functions<Order>, Quantities)], -n3 * l2mu * result[ii][2] - lambda * (n1 * result[ii][0] + n2 * result[ii][1]));
                assign<Override>(output[MATRIX(i+ii, 3, Functions<Order>, Quantities)], -mu * (n2 * result[ii][0] + n1 * result[ii][1]));
                assign<Override>(output[MATRIX(i+ii, 4, Functions<Order>, Quantities)], -mu * (n3 * result[ii][1] + n2 * result[ii][2]));
                assign<Override>(output[MATRIX(i+ii, 5, Functions<Order>, Quantities)], -mu * (n1 * result[ii][2] + n3 * result[ii][0]));
            }
        }
    }
    }
}

template<size_t TileN1, size_t TileN2, size_t TileK, size_t Order, bool Override>
__device__ __forceinline__ void dgkernel2b(real* __restrict__ output, const real* __restrict__ dQ, const real* __restrict__ coordinates, const real* __restrict__ derivativeX,
const real* __restrict__ derivativeY, const real* __restrict__ derivativeZ,
    real lambda, real mu, real rhoD, int offset) {
    real l2mu = lambda + 2 * mu;
    real n1 = coordinates[VECTOR(offset, 9)];
    real n2 = coordinates[VECTOR(offset+1, 9)];
    real n3 = coordinates[VECTOR(offset+2, 9)];
    real n4 = coordinates[VECTOR(offset+3, 9)];
    real n5 = coordinates[VECTOR(offset+4, 9)];
    real n6 = coordinates[VECTOR(offset+5, 9)];
    real n7 = coordinates[VECTOR(offset+6, 9)];
    real n8 = coordinates[VECTOR(offset+7, 9)];
    real n9 = coordinates[VECTOR(offset+8, 9)];

{
    constexpr size_t TileN = TileN1;
        for (int i = 0; i < Functions<Order>; i += TileN) {
        {
            constexpr size_t TileM = 3;
            constexpr int j = 0;
            for (int k = 0; k < Functions<Order+1>; k += TileK) {
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeX[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n2 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n3 * result[ii][2]));
            }
            }
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeY[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n4 * result[ii][0]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n5 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n6 * result[ii][2]));
            }
            }
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeZ[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n7 * result[ii][0]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n8 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n9 * result[ii][2]));
            }
            }
        }}
    }
}
{
    constexpr size_t TileN = TileN1;
        for (int i = 0; i < Functions<Order>; i += TileN) {
        {
            constexpr size_t TileM = 3;
            constexpr int j = 3;
            for (int k = 0; k < Functions<Order+1>; k += TileK) {
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeX[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n2 * result[ii][0] + n3 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0] + n3 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][2] + n2 * result[ii][1]));
            }
            }
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeY[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n5 * result[ii][0] + n6 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n4 * result[ii][0] + n6 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n4 * result[ii][2] + n5 * result[ii][1]));
            }
            }
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeZ[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n8 * result[ii][0] + n9 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n7 * result[ii][0] + n9 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n7 * result[ii][2] + n8 * result[ii][1]));
            }
            }
        }}
    }
}
    {
        constexpr size_t TileN = TileN2;
    for (int i = 0; i < Functions<Order>; i += TileN) {
        {
            constexpr size_t TileM = 3;
            constexpr int j = 6;
            for (int k = 0; k < Functions<Order+1>; k += TileK){
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeX[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 0, Functions<Order>, Quantities)], -n1 * l2mu * result[ii][0] - lambda * (n2 * result[ii][1] + n3 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 1, Functions<Order>, Quantities)], -n2 * l2mu * result[ii][1] - lambda * (n1 * result[ii][0] + n3 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 2, Functions<Order>, Quantities)], -n3 * l2mu * result[ii][2] - lambda * (n1 * result[ii][0] + n2 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 3, Functions<Order>, Quantities)], -mu * (n2 * result[ii][0] + n1 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 4, Functions<Order>, Quantities)], -mu * (n3 * result[ii][1] + n2 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 5, Functions<Order>, Quantities)], -mu * (n1 * result[ii][2] + n3 * result[ii][0]));
            }
            }
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeY[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 0, Functions<Order>, Quantities)], -n1 * l2mu * result[ii][0] - lambda * (n5 * result[ii][1] + n6 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 1, Functions<Order>, Quantities)], -n2 * l2mu * result[ii][1] - lambda * (n4 * result[ii][0] + n6 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 2, Functions<Order>, Quantities)], -n3 * l2mu * result[ii][2] - lambda * (n4 * result[ii][0] + n5 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 3, Functions<Order>, Quantities)], -mu * (n5 * result[ii][0] + n4 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 4, Functions<Order>, Quantities)], -mu * (n6 * result[ii][1] + n5 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 5, Functions<Order>, Quantities)], -mu * (n4 * result[ii][2] + n6 * result[ii][0]));
            }
            }
            {
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK; ++kk) {
                                real dQLocal = dQ[MATRIX(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal
                                    * derivativeZ[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                        }
                    }
                }
            }

            // star 0

            #pragma unroll
            for (int ii = 0; ii < TileN; ++ii) {
                assign<false>(output[MATRIX(i+ii, 0, Functions<Order>, Quantities)], -n1 * l2mu * result[ii][0] - lambda * (n8 * result[ii][1] + n9 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 1, Functions<Order>, Quantities)], -n2 * l2mu * result[ii][1] - lambda * (n7 * result[ii][0] + n9 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 2, Functions<Order>, Quantities)], -n3 * l2mu * result[ii][2] - lambda * (n7 * result[ii][0] + n8 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 3, Functions<Order>, Quantities)], -mu * (n8 * result[ii][0] + n7 * result[ii][1]));
                assign<false>(output[MATRIX(i+ii, 4, Functions<Order>, Quantities)], -mu * (n9 * result[ii][1] + n8 * result[ii][2]));
                assign<false>(output[MATRIX(i+ii, 5, Functions<Order>, Quantities)], -mu * (n7 * result[ii][2] + n9 * result[ii][0]));
            }
            }
        }
        }
    }
    }
}

template<size_t Order, size_t Block, bool Override>
__device__ __forceinline__ void sumkernel(real* __restrict__ I, const real* __restrict__ dQ, real scale) {
    // #pragma unroll
    for (int i = 0; i < Functions<Order>; i += Block) {
#pragma unroll
      for (int ii = 0; ii < Block - (Block % 4); ii += 4) {
#pragma unroll
        for (int j = 0; j < Quantities; ++j) {
          float4 preQ = *((float4*)&dQ[MATRIX4(i+ii, j, Functions<Order>, Quantities)]);
          preQ.x *= scale;
          preQ.y *= scale;
          preQ.z *= scale;
          preQ.w *= scale;
            assign<Override>(*((float4*)&I[MATRIX4(i+ii, j, Functions<6>, Quantities)]), preQ);
        }
      }
#pragma unroll
      for (int ii = Block - (Block % 4); ii < Block; ++ii) {
#pragma unroll
        for (int j = 0; j < Quantities; ++j) {
            assign<Override>(I[MATRIX(i+ii, j, Functions<6>, Quantities)], scale * dQ[MATRIX(i+ii, j, Functions<Order>, Quantities)]);
        }
      }
    }
}

template<size_t Order, size_t Block>
__device__ __forceinline__ void copykernel(real* __restrict__ I, const real* __restrict__ dQ) {
  constexpr size_t TotalSize = Functions<Order> * Quantities;
    for (int i = 0; i < TotalSize; i += Block) {
#pragma unroll
      for (int ii = 0; ii < Block; ii+=4) {
            *((float4*)&I[threadIdx.x * 4 + (i + ii) * Blocksize + OFFSET(TotalSize*Blocksize)]) = *((float4*)&dQ[threadIdx.x * 4 + (i + ii) * Blocksize + OFFSET(TotalSize*Blocksize)]);
      }
    }
    int i = (TotalSize / Block) * Block;
#pragma unroll
    for (int ii = 0; ii < (TotalSize % Block) - ((TotalSize % Block) % 4); ii+=4) {
          *((float4*)&I[threadIdx.x * 4 + (i + ii) * Blocksize + OFFSET(TotalSize*Blocksize)]) = *((float4*)&dQ[threadIdx.x * 4 + (i + ii) * Blocksize + OFFSET(TotalSize*Blocksize)]);
    }
    i = (TotalSize / 4) * 4;
#pragma unroll
    for (int ii = 0; ii < TotalSize % 4; ++ii) {
          *((float*)&I[threadIdx.x + (i + ii) * Blocksize + OFFSET(TotalSize*Blocksize)]) = *((float*)&dQ[threadIdx.x + (i + ii) * Blocksize + OFFSET(TotalSize*Blocksize)]);
    }
}

__device__ __forceinline__ void dgkernel56(real* __restrict__ temp, const real* __restrict__ dQ) {
for (int j = 0; j < Quantities; ++j) {
    real column[56];
    #pragma unroll
    for (int i = 0; i < 56; ++i) {
        column[i] = dQ[MATRIX(i, j, 56, Quantities)];
    }

    temp[MATRIX(0, j + Quantities * 0, 35, Quantities * 3)] = 2.f * column[1] + 1.f * column[5] + 1.f * column[7] + .60000000000000000000f * column[10] + .60000000000000000000f * column[12] + .60000000000000000000f * column[15] + .60000000000000000000f * column[17] + .40000000000000000000f * column[21] + .40000000000000000000f * column[23] + .40000000000000000000f * column[25] + .40000000000000000000f * column[27] + .40000000000000000000f * column[30] + .40000000000000000000f * column[32] + .28571428571428571429f * column[35] + .28571428571428571429f * column[37] + .28571428571428571429f * column[39] + .28571428571428571429f * column[42] + .28571428571428571429f * column[44] + .28571428571428571429f * column[46] + .28571428571428571429f * column[48] + .28571428571428571429f * column[51] + .28571428571428571429f * column[53];
    temp[MATRIX(1, j + Quantities * 0, 35, Quantities * 3)] = 6.f * column[4] + 2.f * column[11] + 2.f * column[14] + 2.8571428571428571429f * column[20] + .85714285714285714286f * column[22] + .85714285714285714286f * column[26] + .85714285714285714286f * column[29] + 1.4285714285714285714f * column[36] + .42857142857142857143f * column[38] + 1.4285714285714285714f * column[41] + .42857142857142857143f * column[43] + .42857142857142857143f * column[47] + .42857142857142857143f * column[50];
    temp[MATRIX(2, j + Quantities * 0, 35, Quantities * 3)] = 3.3333333333333333333f * column[5] + -.66666666666666666667f * column[10] + 2.6666666666666666667f * column[12] + 1.1111111111111111111f * column[15] + .57142857142857142857f * column[21] + 2.f * column[23] + -.28571428571428571429f * column[25] + 1.1428571428571428571f * column[27] + .47619047619047619048f * column[30] + -.47619047619047619048f * column[35] + .80952380952380952381f * column[37] + 1.5238095238095238095f * column[39] + .28571428571428571429f * column[42] + 1.f * column[44] + -.14285714285714285714f * column[46] + .57142857142857142857f * column[48] + .23809523809523809524f * column[51];
    temp[MATRIX(3, j + Quantities * 0, 35, Quantities * 3)] = -.33333333333333333333f * column[5] + 3.f * column[7] + -.33333333333333333333f * column[10] + -.33333333333333333333f * column[12] + 1.2222222222222222222f * column[15] + 2.3333333333333333333f * column[17] + -.28571428571428571429f * column[21] + -.28571428571428571429f * column[23] + .57142857142857142857f * column[25] + .57142857142857142857f * column[27] + 1.2380952380952380952f * column[30] + 1.7142857142857142857f * column[32] + -.23809523809523809524f * column[35] + -.23809523809523809524f * column[37] + -.23809523809523809524f * column[39] + .28571428571428571429f * column[42] + .28571428571428571429f * column[44] + .71428571428571428571f * column[46] + .71428571428571428571f * column[48] + 1.0476190476190476190f * column[51] + 1.2857142857142857143f * column[53];
    temp[MATRIX(4, j + Quantities * 0, 35, Quantities * 3)] = 10.f * column[10] + 2.5000000000000000000f * column[21] + 2.5000000000000000000f * column[25] + 5.8333333333333333333f * column[35] + .83333333333333333333f * column[37] + .83333333333333333333f * column[42] + .83333333333333333333f * column[46];
    temp[MATRIX(5, j + Quantities * 0, 35, Quantities * 3)] = 8.4000000000000000000f * column[11] + -1.5000000000000000000f * column[20] + 4.8000000000000000000f * column[22] + 2.1000000000000000000f * column[26] + 3.1666666666666666667f * column[36] + 2.7000000000000000000f * column[38] + -.50000000000000000000f * column[41] + 1.6000000000000000000f * column[43] + .70000000000000000000f * column[47];
    temp[MATRIX(6, j + Quantities * 0, 35, Quantities * 3)] = .20000000000000000000f * column[10] + 4.2000000000000000000f * column[12] + -1.3000000000000000000f * column[21] + 4.2000000000000000000f * column[23] + .50000000000000000000e-1f * column[25] + 1.0500000000000000000f * column[27] + .33333333333333333333f * column[35] + -.66666666666666666667e-1f * column[37] + 3.6000000000000000000f * column[39] + -.43333333333333333333f * column[42] + 1.4000000000000000000f * column[44] + .16666666666666666667e-1f * column[46] + .35000000000000000000f * column[48];
    temp[MATRIX(7, j + Quantities * 0, 35, Quantities * 3)] = -.40000000000000000000f * column[11] + 8.f * column[14] + -1.f * column[20] + -.30000000000000000000f * column[22] + 2.4000000000000000000f * column[26] + 4.5000000000000000000f * column[29] + -.66666666666666666667f * column[36] + -.20000000000000000000f * column[38] + 3.f * column[41] + .90000000000000000000f * column[43] + 1.8000000000000000000f * column[47] + 2.5000000000000000000f * column[50];
    temp[MATRIX(8, j + Quantities * 0, 35, Quantities * 3)] = .13333333333333333333f * column[10] + -.53333333333333333333f * column[12] + 4.4444444444444444444f * column[15] + -.20000000000000000000f * column[21] + -.70000000000000000000f * column[23] + -.80000000000000000000f * column[25] + 3.2000000000000000000f * column[27] + 2.5000000000000000000f * column[30] + .22222222222222222222f * column[35] + -.37777777777777777778f * column[37] + -.71111111111111111111f * column[39] + .60000000000000000000f * column[42] + 2.1000000000000000000f * column[44] + -.60000000000000000000f * column[46] + 2.4000000000000000000f * column[48] + 1.3888888888888888889f * column[51];
    temp[MATRIX(9, j + Quantities * 0, 35, Quantities * 3)] = .66666666666666666667e-1f * column[10] + .66666666666666666667e-1f * column[12] + -.71111111111111111111f * column[15] + 3.7333333333333333333f * column[17] + .10000000000000000000f * column[21] + .10000000000000000000f * column[23] + -.65000000000000000000f * column[25] + -.65000000000000000000f * column[27] + 1.1000000000000000000f * column[30] + 3.6000000000000000000f * column[32] + .11111111111111111111f * column[35] + .11111111111111111111f * column[37] + .11111111111111111111f * column[39] + -.50000000000000000000f * column[42] + -.50000000000000000000f * column[44] + .25000000000000000000f * column[46] + .25000000000000000000f * column[48] + 1.6111111111111111111f * column[51] + 3.f * column[53];
    temp[MATRIX(10, j + Quantities * 0, 35, Quantities * 3)] = 14.f * column[20] + 2.8000000000000000000f * column[36] + 2.8000000000000000000f * column[41];
    temp[MATRIX(11, j + Quantities * 0, 35, Quantities * 3)] = 12.857142857142857143f * column[21] + -2.f * column[35] + 5.7142857142857142857f * column[37] + 2.5714285714285714286f * column[42];
    temp[MATRIX(12, j + Quantities * 0, 35, Quantities * 3)] = .28571428571428571429f * column[20] + 10.285714285714285714f * column[22] + -3.0857142857142857143f * column[36] + 7.7142857142857142857f * column[38] + .57142857142857142857e-1f * column[41] + 2.0571428571428571429f * column[43];
    temp[MATRIX(13, j + Quantities * 0, 35, Quantities * 3)] = .51428571428571428571f * column[21] + 4.8000000000000000000f * column[23] + -.11428571428571428571f * column[35] + -1.6571428571428571429f * column[37] + 5.4857142857142857143f * column[39] + .10285714285714285714f * column[42] + .96000000000000000000f * column[44];
    temp[MATRIX(14, j + Quantities * 0, 35, Quantities * 3)] = -.35714285714285714286f * column[21] + 12.500000000000000000f * column[25] + -1.5000000000000000000f * column[35] + -.21428571428571428571f * column[37] + 2.9285714285714285714f * column[42] + 5.5000000000000000000f * column[46];
    temp[MATRIX(15, j + Quantities * 0, 35, Quantities * 3)] = .21428571428571428571f * column[20] + -.68571428571428571429f * column[22] + 10.500000000000000000f * column[26] + -.81428571428571428571f * column[36] + -.69428571428571428571f * column[38] + -1.7571428571428571429f * column[41] + 5.6228571428571428571f * column[43] + 4.6200000000000000000f * column[47];
    temp[MATRIX(16, j + Quantities * 0, 35, Quantities * 3)] = .18571428571428571429f * column[21] + -.60000000000000000000f * column[23] + .25000000000000000000f * column[25] + 5.2500000000000000000f * column[27] + -.85714285714285714286e-1f * column[35] + .17142857142857142857e-1f * column[37] + -.92571428571428571429f * column[39] + -1.5228571428571428571f * column[42] + 4.9200000000000000000f * column[44] + .11000000000000000000f * column[46] + 2.3100000000000000000f * column[48];
    temp[MATRIX(17, j + Quantities * 0, 35, Quantities * 3)] = .14285714285714285714f * column[20] + .42857142857142857143e-1f * column[22] + -.85714285714285714286f * column[26] + 9.6428571428571428571f * column[29] + .17142857142857142857f * column[36] + .51428571428571428571e-1f * column[38] + -2.0285714285714285714f * column[41] + -.60857142857142857143f * column[43] + 2.4514285714285714286f * column[47] + 7.0714285714285714286f * column[50];
    temp[MATRIX(18, j + Quantities * 0, 35, Quantities * 3)] = .28571428571428571429e-1f * column[21] + .10000000000000000000f * column[23] + .28571428571428571429f * column[25] + -1.1428571428571428571f * column[27] + 5.3571428571428571429f * column[30] + -.57142857142857142857e-1f * column[35] + .97142857142857142857e-1f * column[37] + .18285714285714285714f * column[39] + -.40571428571428571429f * column[42] + -1.4200000000000000000f * column[44] + -.81714285714285714286f * column[46] + 3.2685714285714285714f * column[48] + 3.9285714285714285714f * column[51];
    temp[MATRIX(19, j + Quantities * 0, 35, Quantities * 3)] = -.14285714285714285714e-1f * column[21] + -.14285714285714285714e-1f * column[23] + .17857142857142857143f * column[25] + .17857142857142857143f * column[27] + -1.0714285714285714286f * column[30] + 4.2857142857142857143f * column[32] + -.28571428571428571429e-1f * column[35] + -.28571428571428571429e-1f * column[37] + -.28571428571428571429e-1f * column[39] + .25428571428571428571f * column[42] + .25428571428571428571f * column[44] + -.86428571428571428571f * column[46] + -.86428571428571428571f * column[48] + .78571428571428571429f * column[51] + 4.7142857142857142857f * column[53];
    temp[MATRIX(20, j + Quantities * 0, 35, Quantities * 3)] = 18.f * column[35];
    temp[MATRIX(21, j + Quantities * 0, 35, Quantities * 3)] = 17.111111111111111111f * column[36];
    temp[MATRIX(22, j + Quantities * 0, 35, Quantities * 3)] = .27777777777777777778f * column[35] + 15.277777777777777778f * column[37];
    temp[MATRIX(23, j + Quantities * 0, 35, Quantities * 3)] = .78571428571428571429f * column[36] + 11.785714285714285714f * column[38];
    temp[MATRIX(24, j + Quantities * 0, 35, Quantities * 3)] = .15873015873015873016e-1f * column[35] + .87301587301587301587f * column[37] + 5.2380952380952380952f * column[39];
    temp[MATRIX(25, j + Quantities * 0, 35, Quantities * 3)] = -.31111111111111111111f * column[36] + 16.800000000000000000f * column[41];
    temp[MATRIX(26, j + Quantities * 0, 35, Quantities * 3)] = .22222222222222222222f * column[35] + -.63492063492063492063f * column[37] + 15.428571428571428571f * column[42];
    temp[MATRIX(27, j + Quantities * 0, 35, Quantities * 3)] = .34285714285714285714f * column[36] + -.85714285714285714286f * column[38] + .34285714285714285714f * column[41] + 12.342857142857142857f * column[43];
    temp[MATRIX(28, j + Quantities * 0, 35, Quantities * 3)] = .12698412698412698413e-1f * column[35] + .18412698412698412698f * column[37] + -.60952380952380952381f * column[39] + .61714285714285714286f * column[42] + 5.7600000000000000000f * column[44];
    temp[MATRIX(29, j + Quantities * 0, 35, Quantities * 3)] = .16666666666666666667f * column[35] + .23809523809523809524e-1f * column[37] + -.76190476190476190476f * column[42] + 14.666666666666666667f * column[46];
    temp[MATRIX(30, j + Quantities * 0, 35, Quantities * 3)] = .90476190476190476190e-1f * column[36] + .77142857142857142857e-1f * column[38] + .45714285714285714286f * column[41] + -1.4628571428571428571f * column[43] + 12.320000000000000000f * column[47];
    temp[MATRIX(31, j + Quantities * 0, 35, Quantities * 3)] = .95238095238095238095e-2f * column[35] + -.19047619047619047619e-2f * column[37] + .10285714285714285714f * column[39] + .39619047619047619048f * column[42] + -1.2800000000000000000f * column[44] + .29333333333333333333f * column[46] + 6.1600000000000000000f * column[48];
    temp[MATRIX(32, j + Quantities * 0, 35, Quantities * 3)] = -.19047619047619047619e-1f * column[36] + -.57142857142857142857e-2f * column[38] + .40000000000000000000f * column[41] + .12000000000000000000f * column[43] + -1.3200000000000000000f * column[47] + 11.f * column[50];
    temp[MATRIX(33, j + Quantities * 0, 35, Quantities * 3)] = .63492063492063492063e-2f * column[35] + -.10793650793650793651e-1f * column[37] + -.20317460317460317460e-1f * column[39] + .80000000000000000000e-1f * column[42] + .28000000000000000000f * column[44] + .44000000000000000000f * column[46] + -1.7600000000000000000f * column[48] + 6.1111111111111111111f * column[51];
    temp[MATRIX(34, j + Quantities * 0, 35, Quantities * 3)] = .31746031746031746032e-2f * column[35] + .31746031746031746032e-2f * column[37] + .31746031746031746032e-2f * column[39] + -.45714285714285714286e-1f * column[42] + -.45714285714285714286e-1f * column[44] + .31428571428571428571f * column[46] + .31428571428571428571f * column[48] + -1.3968253968253968254f * column[51] + 4.7142857142857142857f * column[53];
    temp[MATRIX(0, j + Quantities * 1, 35, Quantities * 3)] = 1.f * column[1] + 3.f * column[2] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + -1.f * column[6] + .50000000000000000000f * column[7] + 1.5000000000000000000f * column[8] + .30000000000000000000f * column[10] + .30000000000000000000f * column[11] + .30000000000000000000f * column[12] + 1.5000000000000000000f * column[13] + .30000000000000000000f * column[14] + .30000000000000000000f * column[15] + -.60000000000000000000f * column[16] + .30000000000000000000f * column[17] + .90000000000000000000f * column[18] + .20000000000000000000f * column[20] + .20000000000000000000f * column[21] + .20000000000000000000f * column[22] + .20000000000000000000f * column[23] + -.80000000000000000000f * column[24] + .20000000000000000000f * column[25] + .20000000000000000000f * column[26] + .20000000000000000000f * column[27] + 1.f * column[28] + .20000000000000000000f * column[29] + .20000000000000000000f * column[30] + -.40000000000000000000f * column[31] + .20000000000000000000f * column[32] + .60000000000000000000f * column[33] + .14285714285714285714f * column[35] + .14285714285714285714f * column[36] + .14285714285714285714f * column[37] + .14285714285714285714f * column[38] + .14285714285714285714f * column[39] + 1.f * column[40] + .14285714285714285714f * column[41] + .14285714285714285714f * column[42] + .14285714285714285714f * column[43] + .14285714285714285714f * column[44] + -.57142857142857142857f * column[45] + .14285714285714285714f * column[46] + .14285714285714285714f * column[47] + .14285714285714285714f * column[48] + .71428571428571428571f * column[49] + .14285714285714285714f * column[50] + .14285714285714285714f * column[51] + -.28571428571428571429f * column[52] + .14285714285714285714f * column[53] + .42857142857142857143f * column[54];
    temp[MATRIX(1, j + Quantities * 1, 35, Quantities * 3)] = 3.f * column[4] + 5.f * column[5] + 2.f * column[10] + 1.f * column[11] + -3.f * column[12] + 1.f * column[14] + 1.6666666666666666667f * column[15] + 1.4285714285714285714f * column[20] + .85714285714285714286f * column[21] + .42857142857142857143f * column[22] + 3.f * column[23] + .85714285714285714286f * column[25] + .42857142857142857143f * column[26] + -1.2857142857142857143f * column[27] + .42857142857142857143f * column[29] + .71428571428571428571f * column[30] + 1.0714285714285714286f * column[35] + .71428571428571428571f * column[36] + .42857142857142857143f * column[37] + .21428571428571428571f * column[38] + -2.4285714285714285714f * column[39] + .71428571428571428571f * column[41] + .42857142857142857143f * column[42] + .21428571428571428571f * column[43] + 1.5000000000000000000f * column[44] + .42857142857142857143f * column[46] + .21428571428571428571f * column[47] + -.64285714285714285714f * column[48] + .21428571428571428571f * column[50] + .35714285714285714286f * column[51];
    temp[MATRIX(2, j + Quantities * 1, 35, Quantities * 3)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 6.6666666666666666667f * column[6] + -.33333333333333333333f * column[10] + .66666666666666666667f * column[11] + 1.3333333333333333333f * column[12] + -1.6666666666666666667f * column[13] + -.11111111111111111111f * column[14] + .55555555555555555556f * column[15] + 2.2222222222222222222f * column[16] + -.28571428571428571429f * column[20] + .28571428571428571429f * column[21] + .71428571428571428571f * column[22] + 1.f * column[23] + 4.f * column[24] + -.14285714285714285714f * column[25] + .28571428571428571429f * column[26] + .57142857142857142857f * column[27] + -.71428571428571428571f * column[28] + -.47619047619047619048e-1f * column[29] + .23809523809523809524f * column[30] + .95238095238095238095f * column[31] + -.23809523809523809524f * column[35] + .11904761904761904762f * column[36] + .40476190476190476190f * column[37] + .61904761904761904762f * column[38] + .76190476190476190476f * column[39] + -1.6666666666666666667f * column[40] + -.14285714285714285714f * column[41] + .14285714285714285714f * column[42] + .35714285714285714286f * column[43] + .50000000000000000000f * column[44] + 2.f * column[45] + -.71428571428571428571e-1f * column[46] + .14285714285714285714f * column[47] + .28571428571428571429f * column[48] + -.35714285714285714286f * column[49] + -.23809523809523809524e-1f * column[50] + .11904761904761904762f * column[51] + .47619047619047619048f * column[52];
    temp[MATRIX(3, j + Quantities * 1, 35, Quantities * 3)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + .33333333333333333333f * column[6] + 1.5000000000000000000f * column[7] + 4.5000000000000000000f * column[8] + -.16666666666666666667f * column[10] + -.16666666666666666667f * column[11] + -.16666666666666666667f * column[12] + -.83333333333333333333f * column[13] + .61111111111111111111f * column[14] + .61111111111111111111f * column[15] + -1.2222222222222222222f * column[16] + 1.1666666666666666667f * column[17] + 3.5000000000000000000f * column[18] + -.14285714285714285714f * column[20] + -.14285714285714285714f * column[21] + -.14285714285714285714f * column[22] + -.14285714285714285714f * column[23] + .57142857142857142857f * column[24] + .28571428571428571429f * column[25] + .28571428571428571429f * column[26] + .28571428571428571429f * column[27] + 1.4285714285714285714f * column[28] + .61904761904761904762f * column[29] + .61904761904761904762f * column[30] + -1.2380952380952380952f * column[31] + .85714285714285714286f * column[32] + 2.5714285714285714286f * column[33] + -.11904761904761904762f * column[35] + -.11904761904761904762f * column[36] + -.11904761904761904762f * column[37] + -.11904761904761904762f * column[38] + -.11904761904761904762f * column[39] + -.83333333333333333333f * column[40] + .14285714285714285714f * column[41] + .14285714285714285714f * column[42] + .14285714285714285714f * column[43] + .14285714285714285714f * column[44] + -.57142857142857142857f * column[45] + .35714285714285714286f * column[46] + .35714285714285714286f * column[47] + .35714285714285714286f * column[48] + 1.7857142857142857143f * column[49] + .52380952380952380952f * column[50] + .52380952380952380952f * column[51] + -1.0476190476190476190f * column[52] + .64285714285714285714f * column[53] + 1.9285714285714285714f * column[54];
    temp[MATRIX(4, j + Quantities * 1, 35, Quantities * 3)] = 5.f * column[10] + 7.f * column[11] + 3.7500000000000000000f * column[20] + 1.2500000000000000000f * column[21] + -5.f * column[22] + 1.2500000000000000000f * column[25] + 1.7500000000000000000f * column[26] + 2.9166666666666666667f * column[35] + 1.2500000000000000000f * column[36] + .41666666666666666667f * column[37] + 4.7500000000000000000f * column[38] + 1.2500000000000000000f * column[41] + .41666666666666666667f * column[42] + -1.6666666666666666667f * column[43] + .41666666666666666667f * column[46] + .58333333333333333333f * column[47];
    temp[MATRIX(5, j + Quantities * 1, 35, Quantities * 3)] = -.60000000000000000000f * column[10] + 4.2000000000000000000f * column[11] + 8.4000000000000000000f * column[12] + -.75000000000000000000f * column[20] + 2.5500000000000000000f * column[21] + 2.4000000000000000000f * column[22] + -4.2000000000000000000f * column[23] + -.15000000000000000000f * column[25] + 1.0500000000000000000f * column[26] + 2.1000000000000000000f * column[27] + -.75000000000000000000f * column[35] + 1.5833333333333333333f * column[36] + 1.9500000000000000000f * column[37] + 1.3500000000000000000f * column[38] + 5.2000000000000000000f * column[39] + -.25000000000000000000f * column[41] + .85000000000000000000f * column[42] + .80000000000000000000f * column[43] + -1.4000000000000000000f * column[44] + -.50000000000000000000e-1f * column[46] + .35000000000000000000f * column[47] + .70000000000000000000f * column[48];
    temp[MATRIX(6, j + Quantities * 1, 35, Quantities * 3)] = .10000000000000000000f * column[10] + -.70000000000000000000f * column[11] + 2.1000000000000000000f * column[12] + 10.500000000000000000f * column[13] + .15000000000000000000f * column[20] + -.65000000000000000000f * column[21] + .50000000000000000000f * column[22] + 2.1000000000000000000f * column[23] + -2.1000000000000000000f * column[24] + .25000000000000000000e-1f * column[25] + -.17500000000000000000f * column[26] + .52500000000000000000f * column[27] + 2.6250000000000000000f * column[28] + .16666666666666666667f * column[35] + -.50000000000000000000f * column[36] + -.33333333333333333333e-1f * column[37] + .90000000000000000000f * column[38] + 1.8000000000000000000f * column[39] + 7.f * column[40] + .50000000000000000000e-1f * column[41] + -.21666666666666666667f * column[42] + .16666666666666666667f * column[43] + .70000000000000000000f * column[44] + -.70000000000000000000f * column[45] + .83333333333333333333e-2f * column[46] + -.58333333333333333333e-1f * column[47] + .17500000000000000000f * column[48] + .87500000000000000000f * column[49];
    temp[MATRIX(7, j + Quantities * 1, 35, Quantities * 3)] = -.40000000000000000000f * column[10] + -.20000000000000000000f * column[11] + .60000000000000000000f * column[12] + 4.f * column[14] + 6.6666666666666666667f * column[15] + -.50000000000000000000f * column[20] + -.30000000000000000000f * column[21] + -.15000000000000000000f * column[22] + -1.0500000000000000000f * column[23] + 2.4000000000000000000f * column[25] + 1.2000000000000000000f * column[26] + -3.6000000000000000000f * column[27] + 2.2500000000000000000f * column[29] + 3.7500000000000000000f * column[30] + -.50000000000000000000f * column[35] + -.33333333333333333333f * column[36] + -.20000000000000000000f * column[37] + -.10000000000000000000f * column[38] + 1.1333333333333333333f * column[39] + 1.5000000000000000000f * column[41] + .90000000000000000000f * column[42] + .45000000000000000000f * column[43] + 3.1500000000000000000f * column[44] + 1.8000000000000000000f * column[46] + .90000000000000000000f * column[47] + -2.7000000000000000000f * column[48] + 1.2500000000000000000f * column[50] + 2.0833333333333333333f * column[51];
    temp[MATRIX(8, j + Quantities * 1, 35, Quantities * 3)] = .66666666666666666667e-1f * column[10] + -.13333333333333333333f * column[11] + -.26666666666666666667f * column[12] + .33333333333333333333f * column[13] + -.44444444444444444444f * column[14] + 2.2222222222222222222f * column[15] + 8.8888888888888888889f * column[16] + .10000000000000000000f * column[20] + -.10000000000000000000f * column[21] + -.25000000000000000000f * column[22] + -.35000000000000000000f * column[23] + -1.4000000000000000000f * column[24] + -.40000000000000000000f * column[25] + .80000000000000000000f * column[26] + 1.6000000000000000000f * column[27] + -2.f * column[28] + -.25000000000000000000f * column[29] + 1.2500000000000000000f * column[30] + 5.f * column[31] + .11111111111111111111f * column[35] + -.55555555555555555556e-1f * column[36] + -.18888888888888888889f * column[37] + -.28888888888888888889f * column[38] + -.35555555555555555556f * column[39] + .77777777777777777778f * column[40] + -.30000000000000000000f * column[41] + .30000000000000000000f * column[42] + .75000000000000000000f * column[43] + 1.0500000000000000000f * column[44] + 4.2000000000000000000f * column[45] + -.30000000000000000000f * column[46] + .60000000000000000000f * column[47] + 1.2000000000000000000f * column[48] + -1.5000000000000000000f * column[49] + -.13888888888888888889f * column[50] + .69444444444444444444f * column[51] + 2.7777777777777777778f * column[52];
    temp[MATRIX(9, j + Quantities * 1, 35, Quantities * 3)] = .33333333333333333333e-1f * column[10] + .33333333333333333333e-1f * column[11] + .33333333333333333333e-1f * column[12] + .16666666666666666667f * column[13] + -.35555555555555555556f * column[14] + -.35555555555555555556f * column[15] + .71111111111111111111f * column[16] + 1.8666666666666666667f * column[17] + 5.6000000000000000000f * column[18] + .50000000000000000000e-1f * column[20] + .50000000000000000000e-1f * column[21] + .50000000000000000000e-1f * column[22] + .50000000000000000000e-1f * column[23] + -.20000000000000000000f * column[24] + -.32500000000000000000f * column[25] + -.32500000000000000000f * column[26] + -.32500000000000000000f * column[27] + -1.6250000000000000000f * column[28] + .55000000000000000000f * column[29] + .55000000000000000000f * column[30] + -1.1000000000000000000f * column[31] + 1.8000000000000000000f * column[32] + 5.4000000000000000000f * column[33] + .55555555555555555556e-1f * column[35] + .55555555555555555556e-1f * column[36] + .55555555555555555556e-1f * column[37] + .55555555555555555556e-1f * column[38] + .55555555555555555556e-1f * column[39] + .38888888888888888889f * column[40] + -.25000000000000000000f * column[41] + -.25000000000000000000f * column[42] + -.25000000000000000000f * column[43] + -.25000000000000000000f * column[44] + 1.f * column[45] + .12500000000000000000f * column[46] + .12500000000000000000f * column[47] + .12500000000000000000f * column[48] + .62500000000000000000f * column[49] + .80555555555555555556f * column[50] + .80555555555555555556f * column[51] + -1.6111111111111111111f * column[52] + 1.5000000000000000000f * column[53] + 4.5000000000000000000f * column[54];
    temp[MATRIX(10, j + Quantities * 1, 35, Quantities * 3)] = 7.f * column[20] + 9.f * column[21] + 5.6000000000000000000f * column[35] + 1.4000000000000000000f * column[36] + -7.f * column[37] + 1.4000000000000000000f * column[41] + 1.8000000000000000000f * column[42];
    temp[MATRIX(11, j + Quantities * 1, 35, Quantities * 3)] = -.71428571428571428571f * column[20] + 6.4285714285714285714f * column[21] + 10.285714285714285714f * column[22] + -1.f * column[35] + 4.5714285714285714286f * column[36] + 2.8571428571428571429f * column[37] + -6.4285714285714285714f * column[38] + -.14285714285714285714f * column[41] + 1.2857142857142857143f * column[42] + 2.0571428571428571429f * column[43];
    temp[MATRIX(12, j + Quantities * 1, 35, Quantities * 3)] = .14285714285714285714f * column[20] + -1.2857142857142857143f * column[21] + 5.1428571428571428571f * column[22] + 12.f * column[23] + .25714285714285714286f * column[35] + -1.5428571428571428571f * column[36] + 2.5714285714285714286f * column[37] + 3.8571428571428571429f * column[38] + -5.1428571428571428571f * column[39] + .28571428571428571429e-1f * column[41] + -.25714285714285714286f * column[42] + 1.0285714285714285714f * column[43] + 2.4000000000000000000f * column[44];
    temp[MATRIX(13, j + Quantities * 1, 35, Quantities * 3)] = -.28571428571428571429e-1f * column[20] + .25714285714285714286f * column[21] + -1.0285714285714285714f * column[22] + 2.4000000000000000000f * column[23] + 14.400000000000000000f * column[24] + -.57142857142857142857e-1f * column[35] + .37142857142857142857f * column[36] + -.82857142857142857143f * column[37] + .17142857142857142857f * column[38] + 2.7428571428571428571f * column[39] + -2.4000000000000000000f * column[40] + -.57142857142857142857e-2f * column[41] + .51428571428571428571e-1f * column[42] + -.20571428571428571429f * column[43] + .48000000000000000000f * column[44] + 2.8800000000000000000f * column[45];
    temp[MATRIX(14, j + Quantities * 1, 35, Quantities * 3)] = -.53571428571428571429f * column[20] + -.17857142857142857143f * column[21] + .71428571428571428571f * column[22] + 6.2500000000000000000f * column[25] + 8.7500000000000000000f * column[26] + -.75000000000000000000f * column[35] + -.32142857142857142857f * column[36] + -.10714285714285714286f * column[37] + -1.2214285714285714286f * column[38] + 4.3928571428571428571f * column[41] + 1.4642857142857142857f * column[42] + -5.8571428571428571429f * column[43] + 2.7500000000000000000f * column[46] + 3.8500000000000000000f * column[47];
    temp[MATRIX(15, j + Quantities * 1, 35, Quantities * 3)] = .10714285714285714286f * column[20] + -.36428571428571428571f * column[21] + -.34285714285714285714f * column[22] + .60000000000000000000f * column[23] + -.75000000000000000000f * column[25] + 5.2500000000000000000f * column[26] + 10.500000000000000000f * column[27] + .19285714285714285714f * column[35] + -.40714285714285714286f * column[36] + -.50142857142857142857f * column[37] + -.34714285714285714286f * column[38] + -1.3371428571428571429f * column[39] + -.87857142857142857143f * column[41] + 2.9871428571428571429f * column[42] + 2.8114285714285714286f * column[43] + -4.9200000000000000000f * column[44] + -.33000000000000000000f * column[46] + 2.3100000000000000000f * column[47] + 4.6200000000000000000f * column[48];
    temp[MATRIX(16, j + Quantities * 1, 35, Quantities * 3)] = -.21428571428571428571e-1f * column[20] + .92857142857142857143e-1f * column[21] + -.71428571428571428571e-1f * column[22] + -.30000000000000000000f * column[23] + .30000000000000000000f * column[24] + .12500000000000000000f * column[25] + -.87500000000000000000f * column[26] + 2.6250000000000000000f * column[27] + 13.125000000000000000f * column[28] + -.42857142857142857143e-1f * column[35] + .12857142857142857143f * column[36] + .85714285714285714286e-2f * column[37] + -.23142857142857142857f * column[38] + -.46285714285714285714f * column[39] + -1.8000000000000000000f * column[40] + .17571428571428571429f * column[41] + -.76142857142857142857f * column[42] + .58571428571428571429f * column[43] + 2.4600000000000000000f * column[44] + -2.4600000000000000000f * column[45] + .55000000000000000000e-1f * column[46] + -.38500000000000000000f * column[47] + 1.1550000000000000000f * column[48] + 5.7750000000000000000f * column[49];
    temp[MATRIX(17, j + Quantities * 1, 35, Quantities * 3)] = .71428571428571428571e-1f * column[20] + .42857142857142857143e-1f * column[21] + .21428571428571428571e-1f * column[22] + .15000000000000000000f * column[23] + -.85714285714285714286f * column[25] + -.42857142857142857143f * column[26] + 1.2857142857142857143f * column[27] + 4.8214285714285714286f * column[29] + 8.0357142857142857143f * column[30] + .12857142857142857143f * column[35] + .85714285714285714286e-1f * column[36] + .51428571428571428571e-1f * column[37] + .25714285714285714286e-1f * column[38] + -.29142857142857142857f * column[39] + -1.0142857142857142857f * column[41] + -.60857142857142857143f * column[42] + -.30428571428571428571f * column[43] + -2.1300000000000000000f * column[44] + 2.4514285714285714286f * column[46] + 1.2257142857142857143f * column[47] + -3.6771428571428571429f * column[48] + 3.5357142857142857143f * column[50] + 5.8928571428571428571f * column[51];
    temp[MATRIX(18, j + Quantities * 1, 35, Quantities * 3)] = -.14285714285714285714e-1f * column[20] + .14285714285714285714e-1f * column[21] + .35714285714285714286e-1f * column[22] + .50000000000000000000e-1f * column[23] + .20000000000000000000f * column[24] + .14285714285714285714f * column[25] + -.28571428571428571429f * column[26] + -.57142857142857142857f * column[27] + .71428571428571428571f * column[28] + -.53571428571428571429f * column[29] + 2.6785714285714285714f * column[30] + 10.714285714285714286f * column[31] + -.28571428571428571429e-1f * column[35] + .14285714285714285714e-1f * column[36] + .48571428571428571429e-1f * column[37] + .74285714285714285714e-1f * column[38] + .91428571428571428571e-1f * column[39] + -.20000000000000000000f * column[40] + .20285714285714285714f * column[41] + -.20285714285714285714f * column[42] + -.50714285714285714286f * column[43] + -.71000000000000000000f * column[44] + -2.8400000000000000000f * column[45] + -.40857142857142857143f * column[46] + .81714285714285714286f * column[47] + 1.6342857142857142857f * column[48] + -2.0428571428571428571f * column[49] + -.39285714285714285714f * column[50] + 1.9642857142857142857f * column[51] + 7.8571428571428571429f * column[52];
    temp[MATRIX(19, j + Quantities * 1, 35, Quantities * 3)] = -.71428571428571428571e-2f * column[20] + -.71428571428571428571e-2f * column[21] + -.71428571428571428571e-2f * column[22] + -.71428571428571428571e-2f * column[23] + .28571428571428571429e-1f * column[24] + .89285714285714285714e-1f * column[25] + .89285714285714285714e-1f * column[26] + .89285714285714285714e-1f * column[27] + .44642857142857142857f * column[28] + -.53571428571428571429f * column[29] + -.53571428571428571429f * column[30] + 1.0714285714285714286f * column[31] + 2.1428571428571428571f * column[32] + 6.4285714285714285714f * column[33] + -.14285714285714285714e-1f * column[35] + -.14285714285714285714e-1f * column[36] + -.14285714285714285714e-1f * column[37] + -.14285714285714285714e-1f * column[38] + -.14285714285714285714e-1f * column[39] + -.10000000000000000000f * column[40] + .12714285714285714286f * column[41] + .12714285714285714286f * column[42] + .12714285714285714286f * column[43] + .12714285714285714286f * column[44] + -.50857142857142857143f * column[45] + -.43214285714285714286f * column[46] + -.43214285714285714286f * column[47] + -.43214285714285714286f * column[48] + -2.1607142857142857143f * column[49] + .39285714285714285714f * column[50] + .39285714285714285714f * column[51] + -.78571428571428571429f * column[52] + 2.3571428571428571429f * column[53] + 7.0714285714285714286f * column[54];
    temp[MATRIX(20, j + Quantities * 1, 35, Quantities * 3)] = 9.f * column[35] + 11.f * column[36];
    temp[MATRIX(21, j + Quantities * 1, 35, Quantities * 3)] = -.77777777777777777778f * column[35] + 8.5555555555555555556f * column[36] + 12.222222222222222222f * column[37];
    temp[MATRIX(22, j + Quantities * 1, 35, Quantities * 3)] = .13888888888888888889f * column[35] + -1.5277777777777777778f * column[36] + 7.6388888888888888889f * column[37] + 13.750000000000000000f * column[38];
    temp[MATRIX(23, j + Quantities * 1, 35, Quantities * 3)] = -.35714285714285714286e-1f * column[35] + .39285714285714285714f * column[36] + -1.9642857142857142857f * column[37] + 5.8928571428571428571f * column[38] + 15.714285714285714286f * column[39];
    temp[MATRIX(24, j + Quantities * 1, 35, Quantities * 3)] = .79365079365079365079e-2f * column[35] + -.87301587301587301587e-1f * column[36] + .43650793650793650794f * column[37] + -1.3095238095238095238f * column[38] + 2.6190476190476190476f * column[39] + 18.333333333333333333f * column[40];
    temp[MATRIX(25, j + Quantities * 1, 35, Quantities * 3)] = -.62222222222222222222f * column[35] + -.15555555555555555556f * column[36] + .77777777777777777778f * column[37] + 8.4000000000000000000f * column[41] + 10.800000000000000000f * column[42];
    temp[MATRIX(26, j + Quantities * 1, 35, Quantities * 3)] = .11111111111111111111f * column[35] + -.50793650793650793651f * column[36] + -.31746031746031746032f * column[37] + .71428571428571428571f * column[38] + -.85714285714285714286f * column[41] + 7.7142857142857142857f * column[42] + 12.342857142857142857f * column[43];
    temp[MATRIX(27, j + Quantities * 1, 35, Quantities * 3)] = -.28571428571428571429e-1f * column[35] + .17142857142857142857f * column[36] + -.28571428571428571429f * column[37] + -.42857142857142857143f * column[38] + .57142857142857142857f * column[39] + .17142857142857142857f * column[41] + -1.5428571428571428571f * column[42] + 6.1714285714285714286f * column[43] + 14.400000000000000000f * column[44];
    temp[MATRIX(28, j + Quantities * 1, 35, Quantities * 3)] = .63492063492063492063e-2f * column[35] + -.41269841269841269841e-1f * column[36] + .92063492063492063492e-1f * column[37] + -.19047619047619047619e-1f * column[38] + -.30476190476190476190f * column[39] + .26666666666666666667f * column[40] + -.34285714285714285714e-1f * column[41] + .30857142857142857143f * column[42] + -1.2342857142857142857f * column[43] + 2.8800000000000000000f * column[44] + 17.280000000000000000f * column[45];
    temp[MATRIX(29, j + Quantities * 1, 35, Quantities * 3)] = .83333333333333333333e-1f * column[35] + .35714285714285714286e-1f * column[36] + .11904761904761904762e-1f * column[37] + .13571428571428571429f * column[38] + -1.1428571428571428571f * column[41] + -.38095238095238095238f * column[42] + 1.5238095238095238095f * column[43] + 7.3333333333333333333f * column[46] + 10.266666666666666667f * column[47];
    temp[MATRIX(30, j + Quantities * 1, 35, Quantities * 3)] = -.21428571428571428571e-1f * column[35] + .45238095238095238095e-1f * column[36] + .55714285714285714286e-1f * column[37] + .38571428571428571429e-1f * column[38] + .14857142857142857143f * column[39] + .22857142857142857143f * column[41] + -.77714285714285714286f * column[42] + -.73142857142857142857f * column[43] + 1.2800000000000000000f * column[44] + -.88000000000000000000f * column[46] + 6.1600000000000000000f * column[47] + 12.320000000000000000f * column[48];
    temp[MATRIX(31, j + Quantities * 1, 35, Quantities * 3)] = .47619047619047619048e-2f * column[35] + -.14285714285714285714e-1f * column[36] + -.95238095238095238095e-3f * column[37] + .25714285714285714286e-1f * column[38] + .51428571428571428571e-1f * column[39] + .20000000000000000000f * column[40] + -.45714285714285714286e-1f * column[41] + .19809523809523809524f * column[42] + -.15238095238095238095f * column[43] + -.64000000000000000000f * column[44] + .64000000000000000000f * column[45] + .14666666666666666667f * column[46] + -1.0266666666666666667f * column[47] + 3.0800000000000000000f * column[48] + 15.400000000000000000f * column[49];
    temp[MATRIX(32, j + Quantities * 1, 35, Quantities * 3)] = -.14285714285714285714e-1f * column[35] + -.95238095238095238095e-2f * column[36] + -.57142857142857142857e-2f * column[37] + -.28571428571428571429e-2f * column[38] + .32380952380952380952e-1f * column[39] + .20000000000000000000f * column[41] + .12000000000000000000f * column[42] + .60000000000000000000e-1f * column[43] + .42000000000000000000f * column[44] + -1.3200000000000000000f * column[46] + -.66000000000000000000f * column[47] + 1.9800000000000000000f * column[48] + 5.5000000000000000000f * column[50] + 9.1666666666666666667f * column[51];
    temp[MATRIX(33, j + Quantities * 1, 35, Quantities * 3)] = .31746031746031746032e-2f * column[35] + -.15873015873015873016e-2f * column[36] + -.53968253968253968254e-2f * column[37] + -.82539682539682539683e-2f * column[38] + -.10158730158730158730e-1f * column[39] + .22222222222222222222e-1f * column[40] + -.40000000000000000000e-1f * column[41] + .40000000000000000000e-1f * column[42] + .10000000000000000000f * column[43] + .14000000000000000000f * column[44] + .56000000000000000000f * column[45] + .22000000000000000000f * column[46] + -.44000000000000000000f * column[47] + -.88000000000000000000f * column[48] + 1.1000000000000000000f * column[49] + -.61111111111111111111f * column[50] + 3.0555555555555555556f * column[51] + 12.222222222222222222f * column[52];
    temp[MATRIX(34, j + Quantities * 1, 35, Quantities * 3)] = .15873015873015873016e-2f * column[35] + .15873015873015873016e-2f * column[36] + .15873015873015873016e-2f * column[37] + .15873015873015873016e-2f * column[38] + .15873015873015873016e-2f * column[39] + .11111111111111111111e-1f * column[40] + -.22857142857142857143e-1f * column[41] + -.22857142857142857143e-1f * column[42] + -.22857142857142857143e-1f * column[43] + -.22857142857142857143e-1f * column[44] + .91428571428571428571e-1f * column[45] + .15714285714285714286f * column[46] + .15714285714285714286f * column[47] + .15714285714285714286f * column[48] + .78571428571428571429f * column[49] + -.69841269841269841270f * column[50] + -.69841269841269841270f * column[51] + 1.3968253968253968254f * column[52] + 2.3571428571428571429f * column[53] + 7.0714285714285714286f * column[54];
    temp[MATRIX(0, j + Quantities * 2, 35, Quantities * 3)] = 1.f * column[1] + 1.f * column[2] + 4.f * column[3] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + .50000000000000000000f * column[6] + .50000000000000000000f * column[7] + .50000000000000000000f * column[8] + -2.5000000000000000000f * column[9] + .30000000000000000000f * column[10] + .30000000000000000000f * column[11] + .30000000000000000000f * column[12] + .30000000000000000000f * column[13] + .30000000000000000000f * column[14] + .30000000000000000000f * column[15] + .30000000000000000000f * column[16] + .30000000000000000000f * column[17] + .30000000000000000000f * column[18] + 3.3000000000000000000f * column[19] + .20000000000000000000f * column[20] + .20000000000000000000f * column[21] + .20000000000000000000f * column[22] + .20000000000000000000f * column[23] + .20000000000000000000f * column[24] + .20000000000000000000f * column[25] + .20000000000000000000f * column[26] + .20000000000000000000f * column[27] + .20000000000000000000f * column[28] + .20000000000000000000f * column[29] + .20000000000000000000f * column[30] + .20000000000000000000f * column[31] + .20000000000000000000f * column[32] + .20000000000000000000f * column[33] + -2.8000000000000000000f * column[34] + .14285714285714285714f * column[35] + .14285714285714285714f * column[36] + .14285714285714285714f * column[37] + .14285714285714285714f * column[38] + .14285714285714285714f * column[39] + .14285714285714285714f * column[40] + .14285714285714285714f * column[41] + .14285714285714285714f * column[42] + .14285714285714285714f * column[43] + .14285714285714285714f * column[44] + .14285714285714285714f * column[45] + .14285714285714285714f * column[46] + .14285714285714285714f * column[47] + .14285714285714285714f * column[48] + .14285714285714285714f * column[49] + .14285714285714285714f * column[50] + .14285714285714285714f * column[51] + .14285714285714285714f * column[52] + .14285714285714285714f * column[53] + .14285714285714285714f * column[54] + 3.1428571428571428571f * column[55];
    temp[MATRIX(1, j + Quantities * 2, 35, Quantities * 3)] = 3.f * column[4] + 1.f * column[5] + 6.f * column[7] + 2.f * column[10] + 1.f * column[11] + .33333333333333333333f * column[12] + 1.f * column[14] + .33333333333333333333f * column[15] + -4.6666666666666666667f * column[17] + 1.4285714285714285714f * column[20] + .85714285714285714286f * column[21] + .42857142857142857143f * column[22] + .14285714285714285714f * column[23] + .85714285714285714286f * column[25] + .42857142857142857143f * column[26] + .14285714285714285714f * column[27] + .42857142857142857143f * column[29] + .14285714285714285714f * column[30] + 5.1428571428571428571f * column[32] + 1.0714285714285714286f * column[35] + .71428571428571428571f * column[36] + .42857142857142857143f * column[37] + .21428571428571428571f * column[38] + .71428571428571428571e-1f * column[39] + .71428571428571428571f * column[41] + .42857142857142857143f * column[42] + .21428571428571428571f * column[43] + .71428571428571428571e-1f * column[44] + .42857142857142857143f * column[46] + .21428571428571428571f * column[47] + .71428571428571428571e-1f * column[48] + .21428571428571428571f * column[50] + .71428571428571428571e-1f * column[51] + -4.9285714285714285714f * column[53];
    temp[MATRIX(2, j + Quantities * 2, 35, Quantities * 3)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 2.6666666666666666667f * column[6] + 6.f * column[8] + -.33333333333333333333f * column[10] + .66666666666666666667f * column[11] + 1.3333333333333333333f * column[12] + 1.6666666666666666667f * column[13] + -.11111111111111111111f * column[14] + .55555555555555555556f * column[15] + .88888888888888888889f * column[16] + -4.6666666666666666667f * column[18] + -.28571428571428571429f * column[20] + .28571428571428571429f * column[21] + .71428571428571428571f * column[22] + 1.f * column[23] + 1.1428571428571428571f * column[24] + -.14285714285714285714f * column[25] + .28571428571428571429f * column[26] + .57142857142857142857f * column[27] + .71428571428571428571f * column[28] + -.47619047619047619048e-1f * column[29] + .23809523809523809524f * column[30] + .38095238095238095238f * column[31] + 5.1428571428571428571f * column[33] + -.23809523809523809524f * column[35] + .11904761904761904762f * column[36] + .40476190476190476190f * column[37] + .61904761904761904762f * column[38] + .76190476190476190476f * column[39] + .83333333333333333333f * column[40] + -.14285714285714285714f * column[41] + .14285714285714285714f * column[42] + .35714285714285714286f * column[43] + .50000000000000000000f * column[44] + .57142857142857142857f * column[45] + -.71428571428571428571e-1f * column[46] + .14285714285714285714f * column[47] + .28571428571428571429f * column[48] + .35714285714285714286f * column[49] + -.23809523809523809524e-1f * column[50] + .11904761904761904762f * column[51] + .19047619047619047619f * column[52] + -4.9285714285714285714f * column[54];
    temp[MATRIX(3, j + Quantities * 2, 35, Quantities * 3)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + -.16666666666666666667f * column[6] + 1.5000000000000000000f * column[7] + 1.5000000000000000000f * column[8] + 7.5000000000000000000f * column[9] + -.16666666666666666667f * column[10] + -.16666666666666666667f * column[11] + -.16666666666666666667f * column[12] + -.16666666666666666667f * column[13] + .61111111111111111111f * column[14] + .61111111111111111111f * column[15] + .61111111111111111111f * column[16] + 1.1666666666666666667f * column[17] + 1.1666666666666666667f * column[18] + -3.5000000000000000000f * column[19] + -.14285714285714285714f * column[20] + -.14285714285714285714f * column[21] + -.14285714285714285714f * column[22] + -.14285714285714285714f * column[23] + -.14285714285714285714f * column[24] + .28571428571428571429f * column[25] + .28571428571428571429f * column[26] + .28571428571428571429f * column[27] + .28571428571428571429f * column[28] + .61904761904761904762f * column[29] + .61904761904761904762f * column[30] + .61904761904761904762f * column[31] + .85714285714285714286f * column[32] + .85714285714285714286f * column[33] + 6.f * column[34] + -.11904761904761904762f * column[35] + -.11904761904761904762f * column[36] + -.11904761904761904762f * column[37] + -.11904761904761904762f * column[38] + -.11904761904761904762f * column[39] + -.11904761904761904762f * column[40] + .14285714285714285714f * column[41] + .14285714285714285714f * column[42] + .14285714285714285714f * column[43] + .14285714285714285714f * column[44] + .14285714285714285714f * column[45] + .35714285714285714286f * column[46] + .35714285714285714286f * column[47] + .35714285714285714286f * column[48] + .35714285714285714286f * column[49] + .52380952380952380952f * column[50] + .52380952380952380952f * column[51] + .52380952380952380952f * column[52] + .64285714285714285714f * column[53] + .64285714285714285714f * column[54] + -4.2857142857142857143f * column[55];
    temp[MATRIX(4, j + Quantities * 2, 35, Quantities * 3)] = 5.f * column[10] + 1.f * column[11] + 8.f * column[14] + 3.7500000000000000000f * column[20] + 1.2500000000000000000f * column[21] + .25000000000000000000f * column[22] + 1.2500000000000000000f * column[25] + .25000000000000000000f * column[26] + -6.7500000000000000000f * column[29] + 2.9166666666666666667f * column[35] + 1.2500000000000000000f * column[36] + .41666666666666666667f * column[37] + .83333333333333333333e-1f * column[38] + 1.2500000000000000000f * column[41] + .41666666666666666667f * column[42] + .83333333333333333333e-1f * column[43] + .41666666666666666667f * column[46] + .83333333333333333333e-1f * column[47] + 7.0833333333333333333f * column[50];
    temp[MATRIX(5, j + Quantities * 2, 35, Quantities * 3)] = -.60000000000000000000f * column[10] + 4.2000000000000000000f * column[11] + 2.4000000000000000000f * column[12] + 8.f * column[15] + -.75000000000000000000f * column[20] + 2.5500000000000000000f * column[21] + 2.4000000000000000000f * column[22] + 1.0500000000000000000f * column[23] + -.15000000000000000000f * column[25] + 1.0500000000000000000f * column[26] + .60000000000000000000f * column[27] + -6.7500000000000000000f * column[30] + -.75000000000000000000f * column[35] + 1.5833333333333333333f * column[36] + 1.9500000000000000000f * column[37] + 1.3500000000000000000f * column[38] + .53333333333333333333f * column[39] + -.25000000000000000000f * column[41] + .85000000000000000000f * column[42] + .80000000000000000000f * column[43] + .35000000000000000000f * column[44] + -.50000000000000000000e-1f * column[46] + .35000000000000000000f * column[47] + .20000000000000000000f * column[48] + 7.0833333333333333333f * column[51];
    temp[MATRIX(6, j + Quantities * 2, 35, Quantities * 3)] = .10000000000000000000f * column[10] + -.70000000000000000000f * column[11] + 2.1000000000000000000f * column[12] + 4.5000000000000000000f * column[13] + 8.f * column[16] + .15000000000000000000f * column[20] + -.65000000000000000000f * column[21] + .50000000000000000000f * column[22] + 2.1000000000000000000f * column[23] + 3.1500000000000000000f * column[24] + .25000000000000000000e-1f * column[25] + -.17500000000000000000f * column[26] + .52500000000000000000f * column[27] + 1.1250000000000000000f * column[28] + -6.7500000000000000000f * column[31] + .16666666666666666667f * column[35] + -.50000000000000000000f * column[36] + -.33333333333333333333e-1f * column[37] + .90000000000000000000f * column[38] + 1.8000000000000000000f * column[39] + 2.3333333333333333333f * column[40] + .50000000000000000000e-1f * column[41] + -.21666666666666666667f * column[42] + .16666666666666666667f * column[43] + .70000000000000000000f * column[44] + 1.0500000000000000000f * column[45] + .83333333333333333333e-2f * column[46] + -.58333333333333333333e-1f * column[47] + .17500000000000000000f * column[48] + .37500000000000000000f * column[49] + 7.0833333333333333333f * column[52];
    temp[MATRIX(7, j + Quantities * 2, 35, Quantities * 3)] = -.40000000000000000000f * column[10] + -.20000000000000000000f * column[11] + -.66666666666666666667e-1f * column[12] + 4.f * column[14] + 1.3333333333333333333f * column[15] + 9.3333333333333333333f * column[17] + -.50000000000000000000f * column[20] + -.30000000000000000000f * column[21] + -.15000000000000000000f * column[22] + -.50000000000000000000e-1f * column[23] + 2.4000000000000000000f * column[25] + 1.2000000000000000000f * column[26] + .40000000000000000000f * column[27] + 2.2500000000000000000f * column[29] + .75000000000000000000f * column[30] + -6.f * column[32] + -.50000000000000000000f * column[35] + -.33333333333333333333f * column[36] + -.20000000000000000000f * column[37] + -.10000000000000000000f * column[38] + -.33333333333333333333e-1f * column[39] + 1.5000000000000000000f * column[41] + .90000000000000000000f * column[42] + .45000000000000000000f * column[43] + .15000000000000000000f * column[44] + 1.8000000000000000000f * column[46] + .90000000000000000000f * column[47] + .30000000000000000000f * column[48] + 1.2500000000000000000f * column[50] + .41666666666666666667f * column[51] + 7.5000000000000000000f * column[53];
    temp[MATRIX(8, j + Quantities * 2, 35, Quantities * 3)] = .66666666666666666667e-1f * column[10] + -.13333333333333333333f * column[11] + -.26666666666666666667f * column[12] + -.33333333333333333333f * column[13] + -.44444444444444444444f * column[14] + 2.2222222222222222222f * column[15] + 3.5555555555555555556f * column[16] + 9.3333333333333333333f * column[18] + .10000000000000000000f * column[20] + -.10000000000000000000f * column[21] + -.25000000000000000000f * column[22] + -.35000000000000000000f * column[23] + -.40000000000000000000f * column[24] + -.40000000000000000000f * column[25] + .80000000000000000000f * column[26] + 1.6000000000000000000f * column[27] + 2.f * column[28] + -.25000000000000000000f * column[29] + 1.2500000000000000000f * column[30] + 2.f * column[31] + -6.f * column[33] + .11111111111111111111f * column[35] + -.55555555555555555556e-1f * column[36] + -.18888888888888888889f * column[37] + -.28888888888888888889f * column[38] + -.35555555555555555556f * column[39] + -.38888888888888888889f * column[40] + -.30000000000000000000f * column[41] + .30000000000000000000f * column[42] + .75000000000000000000f * column[43] + 1.0500000000000000000f * column[44] + 1.2000000000000000000f * column[45] + -.30000000000000000000f * column[46] + .60000000000000000000f * column[47] + 1.2000000000000000000f * column[48] + 1.5000000000000000000f * column[49] + -.13888888888888888889f * column[50] + .69444444444444444444f * column[51] + 1.1111111111111111111f * column[52] + 7.5000000000000000000f * column[54];
    temp[MATRIX(9, j + Quantities * 2, 35, Quantities * 3)] = .33333333333333333333e-1f * column[10] + .33333333333333333333e-1f * column[11] + .33333333333333333333e-1f * column[12] + .33333333333333333333e-1f * column[13] + -.35555555555555555556f * column[14] + -.35555555555555555556f * column[15] + -.35555555555555555556f * column[16] + 1.8666666666666666667f * column[17] + 1.8666666666666666667f * column[18] + 11.200000000000000000f * column[19] + .50000000000000000000e-1f * column[20] + .50000000000000000000e-1f * column[21] + .50000000000000000000e-1f * column[22] + .50000000000000000000e-1f * column[23] + .50000000000000000000e-1f * column[24] + -.32500000000000000000f * column[25] + -.32500000000000000000f * column[26] + -.32500000000000000000f * column[27] + -.32500000000000000000f * column[28] + .55000000000000000000f * column[29] + .55000000000000000000f * column[30] + .55000000000000000000f * column[31] + 1.8000000000000000000f * column[32] + 1.8000000000000000000f * column[33] + -4.2000000000000000000f * column[34] + .55555555555555555556e-1f * column[35] + .55555555555555555556e-1f * column[36] + .55555555555555555556e-1f * column[37] + .55555555555555555556e-1f * column[38] + .55555555555555555556e-1f * column[39] + .55555555555555555556e-1f * column[40] + -.25000000000000000000f * column[41] + -.25000000000000000000f * column[42] + -.25000000000000000000f * column[43] + -.25000000000000000000f * column[44] + -.25000000000000000000f * column[45] + .12500000000000000000f * column[46] + .12500000000000000000f * column[47] + .12500000000000000000f * column[48] + .12500000000000000000f * column[49] + .80555555555555555556f * column[50] + .80555555555555555556f * column[51] + .80555555555555555556f * column[52] + 1.5000000000000000000f * column[53] + 1.5000000000000000000f * column[54] + 9.f * column[55];
    temp[MATRIX(10, j + Quantities * 2, 35, Quantities * 3)] = 7.f * column[20] + 1.f * column[21] + 10.f * column[25] + 5.6000000000000000000f * column[35] + 1.4000000000000000000f * column[36] + .20000000000000000000f * column[37] + 1.4000000000000000000f * column[41] + .20000000000000000000f * column[42] + -8.8000000000000000000f * column[46];
    temp[MATRIX(11, j + Quantities * 2, 35, Quantities * 3)] = -.71428571428571428571f * column[20] + 6.4285714285714285714f * column[21] + 2.2857142857142857143f * column[22] + 10.f * column[26] + -1.f * column[35] + 4.5714285714285714286f * column[36] + 2.8571428571428571429f * column[37] + .77142857142857142857f * column[38] + -.14285714285714285714f * column[41] + 1.2857142857142857143f * column[42] + .45714285714285714286f * column[43] + -8.8000000000000000000f * column[47];
    temp[MATRIX(12, j + Quantities * 2, 35, Quantities * 3)] = .14285714285714285714f * column[20] + -1.2857142857142857143f * column[21] + 5.1428571428571428571f * column[22] + 4.f * column[23] + 10.f * column[27] + .25714285714285714286f * column[35] + -1.5428571428571428571f * column[36] + 2.5714285714285714286f * column[37] + 3.8571428571428571429f * column[38] + 2.0571428571428571429f * column[39] + .28571428571428571429e-1f * column[41] + -.25714285714285714286f * column[42] + 1.0285714285714285714f * column[43] + .80000000000000000000f * column[44] + -8.8000000000000000000f * column[48];
    temp[MATRIX(13, j + Quantities * 2, 35, Quantities * 3)] = -.28571428571428571429e-1f * column[20] + .25714285714285714286f * column[21] + -1.0285714285714285714f * column[22] + 2.4000000000000000000f * column[23] + 6.4000000000000000000f * column[24] + 10.f * column[28] + -.57142857142857142857e-1f * column[35] + .37142857142857142857f * column[36] + -.82857142857142857143f * column[37] + .17142857142857142857f * column[38] + 2.7428571428571428571f * column[39] + 4.8000000000000000000f * column[40] + -.57142857142857142857e-2f * column[41] + .51428571428571428571e-1f * column[42] + -.20571428571428571429f * column[43] + .48000000000000000000f * column[44] + 1.2800000000000000000f * column[45] + -8.8000000000000000000f * column[49];
    temp[MATRIX(14, j + Quantities * 2, 35, Quantities * 3)] = -.53571428571428571429f * column[20] + -.17857142857142857143f * column[21] + -.35714285714285714286e-1f * column[22] + 6.2500000000000000000f * column[25] + 1.2500000000000000000f * column[26] + 11.250000000000000000f * column[29] + -.75000000000000000000f * column[35] + -.32142857142857142857f * column[36] + -.10714285714285714286f * column[37] + -.21428571428571428571e-1f * column[38] + 4.3928571428571428571f * column[41] + 1.4642857142857142857f * column[42] + .29285714285714285714f * column[43] + 2.7500000000000000000f * column[46] + .55000000000000000000f * column[47] + -8.2500000000000000000f * column[50];
    temp[MATRIX(15, j + Quantities * 2, 35, Quantities * 3)] = .10714285714285714286f * column[20] + -.36428571428571428571f * column[21] + -.34285714285714285714f * column[22] + -.15000000000000000000f * column[23] + -.75000000000000000000f * column[25] + 5.2500000000000000000f * column[26] + 3.f * column[27] + 11.250000000000000000f * column[30] + .19285714285714285714f * column[35] + -.40714285714285714286f * column[36] + -.50142857142857142857f * column[37] + -.34714285714285714286f * column[38] + -.13714285714285714286f * column[39] + -.87857142857142857143f * column[41] + 2.9871428571428571429f * column[42] + 2.8114285714285714286f * column[43] + 1.2300000000000000000f * column[44] + -.33000000000000000000f * column[46] + 2.3100000000000000000f * column[47] + 1.3200000000000000000f * column[48] + -8.2500000000000000000f * column[51];
    temp[MATRIX(16, j + Quantities * 2, 35, Quantities * 3)] = -.21428571428571428571e-1f * column[20] + .92857142857142857143e-1f * column[21] + -.71428571428571428571e-1f * column[22] + -.30000000000000000000f * column[23] + -.45000000000000000000f * column[24] + .12500000000000000000f * column[25] + -.87500000000000000000f * column[26] + 2.6250000000000000000f * column[27] + 5.6250000000000000000f * column[28] + 11.250000000000000000f * column[31] + -.42857142857142857143e-1f * column[35] + .12857142857142857143f * column[36] + .85714285714285714286e-2f * column[37] + -.23142857142857142857f * column[38] + -.46285714285714285714f * column[39] + -.60000000000000000000f * column[40] + .17571428571428571429f * column[41] + -.76142857142857142857f * column[42] + .58571428571428571429f * column[43] + 2.4600000000000000000f * column[44] + 3.6900000000000000000f * column[45] + .55000000000000000000e-1f * column[46] + -.38500000000000000000f * column[47] + 1.1550000000000000000f * column[48] + 2.4750000000000000000f * column[49] + -8.2500000000000000000f * column[52];
    temp[MATRIX(17, j + Quantities * 2, 35, Quantities * 3)] = .71428571428571428571e-1f * column[20] + .42857142857142857143e-1f * column[21] + .21428571428571428571e-1f * column[22] + .71428571428571428571e-2f * column[23] + -.85714285714285714286f * column[25] + -.42857142857142857143f * column[26] + -.14285714285714285714f * column[27] + 4.8214285714285714286f * column[29] + 1.6071428571428571429f * column[30] + 12.857142857142857143f * column[32] + .12857142857142857143f * column[35] + .85714285714285714286e-1f * column[36] + .51428571428571428571e-1f * column[37] + .25714285714285714286e-1f * column[38] + .85714285714285714286e-2f * column[39] + -1.0142857142857142857f * column[41] + -.60857142857142857143f * column[42] + -.30428571428571428571f * column[43] + -.10142857142857142857f * column[44] + 2.4514285714285714286f * column[46] + 1.2257142857142857143f * column[47] + .40857142857142857143f * column[48] + 3.5357142857142857143f * column[50] + 1.1785714285714285714f * column[51] + -7.0714285714285714286f * column[53];
    temp[MATRIX(18, j + Quantities * 2, 35, Quantities * 3)] = -.14285714285714285714e-1f * column[20] + .14285714285714285714e-1f * column[21] + .35714285714285714286e-1f * column[22] + .50000000000000000000e-1f * column[23] + .57142857142857142857e-1f * column[24] + .14285714285714285714f * column[25] + -.28571428571428571429f * column[26] + -.57142857142857142857f * column[27] + -.71428571428571428571f * column[28] + -.53571428571428571429f * column[29] + 2.6785714285714285714f * column[30] + 4.2857142857142857143f * column[31] + 12.857142857142857143f * column[33] + -.28571428571428571429e-1f * column[35] + .14285714285714285714e-1f * column[36] + .48571428571428571429e-1f * column[37] + .74285714285714285714e-1f * column[38] + .91428571428571428571e-1f * column[39] + .10000000000000000000f * column[40] + .20285714285714285714f * column[41] + -.20285714285714285714f * column[42] + -.50714285714285714286f * column[43] + -.71000000000000000000f * column[44] + -.81142857142857142857f * column[45] + -.40857142857142857143f * column[46] + .81714285714285714286f * column[47] + 1.6342857142857142857f * column[48] + 2.0428571428571428571f * column[49] + -.39285714285714285714f * column[50] + 1.9642857142857142857f * column[51] + 3.1428571428571428571f * column[52] + -7.0714285714285714286f * column[54];
    temp[MATRIX(19, j + Quantities * 2, 35, Quantities * 3)] = -.71428571428571428571e-2f * column[20] + -.71428571428571428571e-2f * column[21] + -.71428571428571428571e-2f * column[22] + -.71428571428571428571e-2f * column[23] + -.71428571428571428571e-2f * column[24] + .89285714285714285714e-1f * column[25] + .89285714285714285714e-1f * column[26] + .89285714285714285714e-1f * column[27] + .89285714285714285714e-1f * column[28] + -.53571428571428571429f * column[29] + -.53571428571428571429f * column[30] + -.53571428571428571429f * column[31] + 2.1428571428571428571f * column[32] + 2.1428571428571428571f * column[33] + 15.f * column[34] + -.14285714285714285714e-1f * column[35] + -.14285714285714285714e-1f * column[36] + -.14285714285714285714e-1f * column[37] + -.14285714285714285714e-1f * column[38] + -.14285714285714285714e-1f * column[39] + -.14285714285714285714e-1f * column[40] + .12714285714285714286f * column[41] + .12714285714285714286f * column[42] + .12714285714285714286f * column[43] + .12714285714285714286f * column[44] + .12714285714285714286f * column[45] + -.43214285714285714286f * column[46] + -.43214285714285714286f * column[47] + -.43214285714285714286f * column[48] + -.43214285714285714286f * column[49] + .39285714285714285714f * column[50] + .39285714285714285714f * column[51] + .39285714285714285714f * column[52] + 2.3571428571428571429f * column[53] + 2.3571428571428571429f * column[54] + -4.7142857142857142857f * column[55];
    temp[MATRIX(20, j + Quantities * 2, 35, Quantities * 3)] = 9.f * column[35] + 1.f * column[36] + 12.f * column[41];
    temp[MATRIX(21, j + Quantities * 2, 35, Quantities * 3)] = -.77777777777777777778f * column[35] + 8.5555555555555555556f * column[36] + 2.2222222222222222222f * column[37] + 12.f * column[42];
    temp[MATRIX(22, j + Quantities * 2, 35, Quantities * 3)] = .13888888888888888889f * column[35] + -1.5277777777777777778f * column[36] + 7.6388888888888888889f * column[37] + 3.7500000000000000000f * column[38] + 12.f * column[43];
    temp[MATRIX(23, j + Quantities * 2, 35, Quantities * 3)] = -.35714285714285714286e-1f * column[35] + .39285714285714285714f * column[36] + -1.9642857142857142857f * column[37] + 5.8928571428571428571f * column[38] + 5.7142857142857142857f * column[39] + 12.f * column[44];
    temp[MATRIX(24, j + Quantities * 2, 35, Quantities * 3)] = .79365079365079365079e-2f * column[35] + -.87301587301587301587e-1f * column[36] + .43650793650793650794f * column[37] + -1.3095238095238095238f * column[38] + 2.6190476190476190476f * column[39] + 8.3333333333333333333f * column[40] + 12.f * column[45];
    temp[MATRIX(25, j + Quantities * 2, 35, Quantities * 3)] = -.62222222222222222222f * column[35] + -.15555555555555555556f * column[36] + -.22222222222222222222e-1f * column[37] + 8.4000000000000000000f * column[41] + 1.2000000000000000000f * column[42] + 13.200000000000000000f * column[46];
    temp[MATRIX(26, j + Quantities * 2, 35, Quantities * 3)] = .11111111111111111111f * column[35] + -.50793650793650793651f * column[36] + -.31746031746031746032f * column[37] + -.85714285714285714286e-1f * column[38] + -.85714285714285714286f * column[41] + 7.7142857142857142857f * column[42] + 2.7428571428571428571f * column[43] + 13.200000000000000000f * column[47];
    temp[MATRIX(27, j + Quantities * 2, 35, Quantities * 3)] = -.28571428571428571429e-1f * column[35] + .17142857142857142857f * column[36] + -.28571428571428571429f * column[37] + -.42857142857142857143f * column[38] + -.22857142857142857143f * column[39] + .17142857142857142857f * column[41] + -1.5428571428571428571f * column[42] + 6.1714285714285714286f * column[43] + 4.8000000000000000000f * column[44] + 13.200000000000000000f * column[48];
    temp[MATRIX(28, j + Quantities * 2, 35, Quantities * 3)] = .63492063492063492063e-2f * column[35] + -.41269841269841269841e-1f * column[36] + .92063492063492063492e-1f * column[37] + -.19047619047619047619e-1f * column[38] + -.30476190476190476190f * column[39] + -.53333333333333333333f * column[40] + -.34285714285714285714e-1f * column[41] + .30857142857142857143f * column[42] + -1.2342857142857142857f * column[43] + 2.8800000000000000000f * column[44] + 7.6800000000000000000f * column[45] + 13.200000000000000000f * column[49];
    temp[MATRIX(29, j + Quantities * 2, 35, Quantities * 3)] = .83333333333333333333e-1f * column[35] + .35714285714285714286e-1f * column[36] + .11904761904761904762e-1f * column[37] + .23809523809523809524e-2f * column[38] + -1.1428571428571428571f * column[41] + -.38095238095238095238f * column[42] + -.76190476190476190476e-1f * column[43] + 7.3333333333333333333f * column[46] + 1.4666666666666666667f * column[47] + 14.666666666666666667f * column[50];
    temp[MATRIX(30, j + Quantities * 2, 35, Quantities * 3)] = -.21428571428571428571e-1f * column[35] + .45238095238095238095e-1f * column[36] + .55714285714285714286e-1f * column[37] + .38571428571428571429e-1f * column[38] + .15238095238095238095e-1f * column[39] + .22857142857142857143f * column[41] + -.77714285714285714286f * column[42] + -.73142857142857142857f * column[43] + -.32000000000000000000f * column[44] + -.88000000000000000000f * column[46] + 6.1600000000000000000f * column[47] + 3.5200000000000000000f * column[48] + 14.666666666666666667f * column[51];
    temp[MATRIX(31, j + Quantities * 2, 35, Quantities * 3)] = .47619047619047619048e-2f * column[35] + -.14285714285714285714e-1f * column[36] + -.95238095238095238095e-3f * column[37] + .25714285714285714286e-1f * column[38] + .51428571428571428571e-1f * column[39] + .66666666666666666667e-1f * column[40] + -.45714285714285714286e-1f * column[41] + .19809523809523809524f * column[42] + -.15238095238095238095f * column[43] + -.64000000000000000000f * column[44] + -.96000000000000000000f * column[45] + .14666666666666666667f * column[46] + -1.0266666666666666667f * column[47] + 3.0800000000000000000f * column[48] + 6.6000000000000000000f * column[49] + 14.666666666666666667f * column[52];
    temp[MATRIX(32, j + Quantities * 2, 35, Quantities * 3)] = -.14285714285714285714e-1f * column[35] + -.95238095238095238095e-2f * column[36] + -.57142857142857142857e-2f * column[37] + -.28571428571428571429e-2f * column[38] + -.95238095238095238095e-3f * column[39] + .20000000000000000000f * column[41] + .12000000000000000000f * column[42] + .60000000000000000000e-1f * column[43] + .20000000000000000000e-1f * column[44] + -1.3200000000000000000f * column[46] + -.66000000000000000000f * column[47] + -.22000000000000000000f * column[48] + 5.5000000000000000000f * column[50] + 1.8333333333333333333f * column[51] + 16.500000000000000000f * column[53];
    temp[MATRIX(33, j + Quantities * 2, 35, Quantities * 3)] = .31746031746031746032e-2f * column[35] + -.15873015873015873016e-2f * column[36] + -.53968253968253968254e-2f * column[37] + -.82539682539682539683e-2f * column[38] + -.10158730158730158730e-1f * column[39] + -.11111111111111111111e-1f * column[40] + -.40000000000000000000e-1f * column[41] + .40000000000000000000e-1f * column[42] + .10000000000000000000f * column[43] + .14000000000000000000f * column[44] + .16000000000000000000f * column[45] + .22000000000000000000f * column[46] + -.44000000000000000000f * column[47] + -.88000000000000000000f * column[48] + -1.1000000000000000000f * column[49] + -.61111111111111111111f * column[50] + 3.0555555555555555556f * column[51] + 4.8888888888888888889f * column[52] + 16.500000000000000000f * column[54];
    temp[MATRIX(34, j + Quantities * 2, 35, Quantities * 3)] = .15873015873015873016e-2f * column[35] + .15873015873015873016e-2f * column[36] + .15873015873015873016e-2f * column[37] + .15873015873015873016e-2f * column[38] + .15873015873015873016e-2f * column[39] + .15873015873015873016e-2f * column[40] + -.22857142857142857143e-1f * column[41] + -.22857142857142857143e-1f * column[42] + -.22857142857142857143e-1f * column[43] + -.22857142857142857143e-1f * column[44] + -.22857142857142857143e-1f * column[45] + .15714285714285714286f * column[46] + .15714285714285714286f * column[47] + .15714285714285714286f * column[48] + .15714285714285714286f * column[49] + -.69841269841269841270f * column[50] + -.69841269841269841270f * column[51] + -.69841269841269841270f * column[52] + 2.3571428571428571429f * column[53] + 2.3571428571428571429f * column[54] + 18.857142857142857143f * column[55];
    }
}
__device__ __forceinline__ void dgkernel35(real* __restrict__ temp, const real* __restrict__ dQ) {
for (int j = 0; j < Quantities; ++j) {
    real column[35];
    #pragma unroll
    for (int i = 0; i < 35; ++i) {
        column[i] = dQ[MATRIX(i, j, 35, Quantities)];
    }

    temp[MATRIX(0, j + Quantities * 0, 20, Quantities * 3)] = 2.f * column[1] + 1.f * column[5] + 1.f * column[7] + .60000000000000000000f * column[10] + .60000000000000000000f * column[12] + .60000000000000000000f * column[15] + .60000000000000000000f * column[17] + .40000000000000000000f * column[21] + .40000000000000000000f * column[23] + .40000000000000000000f * column[25] + .40000000000000000000f * column[27] + .40000000000000000000f * column[30] + .40000000000000000000f * column[32];
    temp[MATRIX(1, j + Quantities * 0, 20, Quantities * 3)] = 6.f * column[4] + 2.f * column[11] + 2.f * column[14] + 2.8571428571428571429f * column[20] + .85714285714285714286f * column[22] + .85714285714285714286f * column[26] + .85714285714285714286f * column[29];
    temp[MATRIX(2, j + Quantities * 0, 20, Quantities * 3)] = 3.3333333333333333333f * column[5] + -.66666666666666666667f * column[10] + 2.6666666666666666667f * column[12] + 1.1111111111111111111f * column[15] + .57142857142857142857f * column[21] + 2.f * column[23] + -.28571428571428571429f * column[25] + 1.1428571428571428571f * column[27] + .47619047619047619048f * column[30];
    temp[MATRIX(3, j + Quantities * 0, 20, Quantities * 3)] = -.33333333333333333333f * column[5] + 3.f * column[7] + -.33333333333333333333f * column[10] + -.33333333333333333333f * column[12] + 1.2222222222222222222f * column[15] + 2.3333333333333333333f * column[17] + -.28571428571428571429f * column[21] + -.28571428571428571429f * column[23] + .57142857142857142857f * column[25] + .57142857142857142857f * column[27] + 1.2380952380952380952f * column[30] + 1.7142857142857142857f * column[32];
    temp[MATRIX(4, j + Quantities * 0, 20, Quantities * 3)] = 10.f * column[10] + 2.5000000000000000000f * column[21] + 2.5000000000000000000f * column[25];
    temp[MATRIX(5, j + Quantities * 0, 20, Quantities * 3)] = 8.4000000000000000000f * column[11] + -1.5000000000000000000f * column[20] + 4.8000000000000000000f * column[22] + 2.1000000000000000000f * column[26];
    temp[MATRIX(6, j + Quantities * 0, 20, Quantities * 3)] = .20000000000000000000f * column[10] + 4.2000000000000000000f * column[12] + -1.3000000000000000000f * column[21] + 4.2000000000000000000f * column[23] + .50000000000000000000e-1f * column[25] + 1.0500000000000000000f * column[27];
    temp[MATRIX(7, j + Quantities * 0, 20, Quantities * 3)] = -.40000000000000000000f * column[11] + 8.f * column[14] + -1.f * column[20] + -.30000000000000000000f * column[22] + 2.4000000000000000000f * column[26] + 4.5000000000000000000f * column[29];
    temp[MATRIX(8, j + Quantities * 0, 20, Quantities * 3)] = .13333333333333333333f * column[10] + -.53333333333333333333f * column[12] + 4.4444444444444444444f * column[15] + -.20000000000000000000f * column[21] + -.70000000000000000000f * column[23] + -.80000000000000000000f * column[25] + 3.2000000000000000000f * column[27] + 2.5000000000000000000f * column[30];
    temp[MATRIX(9, j + Quantities * 0, 20, Quantities * 3)] = .66666666666666666667e-1f * column[10] + .66666666666666666667e-1f * column[12] + -.71111111111111111111f * column[15] + 3.7333333333333333333f * column[17] + .10000000000000000000f * column[21] + .10000000000000000000f * column[23] + -.65000000000000000000f * column[25] + -.65000000000000000000f * column[27] + 1.1000000000000000000f * column[30] + 3.6000000000000000000f * column[32];
    temp[MATRIX(10, j + Quantities * 0, 20, Quantities * 3)] = 14.f * column[20];
    temp[MATRIX(11, j + Quantities * 0, 20, Quantities * 3)] = 12.857142857142857143f * column[21];
    temp[MATRIX(12, j + Quantities * 0, 20, Quantities * 3)] = .28571428571428571429f * column[20] + 10.285714285714285714f * column[22];
    temp[MATRIX(13, j + Quantities * 0, 20, Quantities * 3)] = .51428571428571428571f * column[21] + 4.8000000000000000000f * column[23];
    temp[MATRIX(14, j + Quantities * 0, 20, Quantities * 3)] = -.35714285714285714286f * column[21] + 12.500000000000000000f * column[25];
    temp[MATRIX(15, j + Quantities * 0, 20, Quantities * 3)] = .21428571428571428571f * column[20] + -.68571428571428571429f * column[22] + 10.500000000000000000f * column[26];
    temp[MATRIX(16, j + Quantities * 0, 20, Quantities * 3)] = .18571428571428571429f * column[21] + -.60000000000000000000f * column[23] + .25000000000000000000f * column[25] + 5.2500000000000000000f * column[27];
    temp[MATRIX(17, j + Quantities * 0, 20, Quantities * 3)] = .14285714285714285714f * column[20] + .42857142857142857143e-1f * column[22] + -.85714285714285714286f * column[26] + 9.6428571428571428571f * column[29];
    temp[MATRIX(18, j + Quantities * 0, 20, Quantities * 3)] = .28571428571428571429e-1f * column[21] + .10000000000000000000f * column[23] + .28571428571428571429f * column[25] + -1.1428571428571428571f * column[27] + 5.3571428571428571429f * column[30];
    temp[MATRIX(19, j + Quantities * 0, 20, Quantities * 3)] = -.14285714285714285714e-1f * column[21] + -.14285714285714285714e-1f * column[23] + .17857142857142857143f * column[25] + .17857142857142857143f * column[27] + -1.0714285714285714286f * column[30] + 4.2857142857142857143f * column[32];
    temp[MATRIX(0, j + Quantities * 1, 20, Quantities * 3)] = 1.f * column[1] + 3.f * column[2] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + -1.f * column[6] + .50000000000000000000f * column[7] + 1.5000000000000000000f * column[8] + .30000000000000000000f * column[10] + .30000000000000000000f * column[11] + .30000000000000000000f * column[12] + 1.5000000000000000000f * column[13] + .30000000000000000000f * column[14] + .30000000000000000000f * column[15] + -.60000000000000000000f * column[16] + .30000000000000000000f * column[17] + .90000000000000000000f * column[18] + .20000000000000000000f * column[20] + .20000000000000000000f * column[21] + .20000000000000000000f * column[22] + .20000000000000000000f * column[23] + -.80000000000000000000f * column[24] + .20000000000000000000f * column[25] + .20000000000000000000f * column[26] + .20000000000000000000f * column[27] + 1.f * column[28] + .20000000000000000000f * column[29] + .20000000000000000000f * column[30] + -.40000000000000000000f * column[31] + .20000000000000000000f * column[32] + .60000000000000000000f * column[33];
    temp[MATRIX(1, j + Quantities * 1, 20, Quantities * 3)] = 3.f * column[4] + 5.f * column[5] + 2.f * column[10] + 1.f * column[11] + -3.f * column[12] + 1.f * column[14] + 1.6666666666666666667f * column[15] + 1.4285714285714285714f * column[20] + .85714285714285714286f * column[21] + .42857142857142857143f * column[22] + 3.f * column[23] + .85714285714285714286f * column[25] + .42857142857142857143f * column[26] + -1.2857142857142857143f * column[27] + .42857142857142857143f * column[29] + .71428571428571428571f * column[30];
    temp[MATRIX(2, j + Quantities * 1, 20, Quantities * 3)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 6.6666666666666666667f * column[6] + -.33333333333333333333f * column[10] + .66666666666666666667f * column[11] + 1.3333333333333333333f * column[12] + -1.6666666666666666667f * column[13] + -.11111111111111111111f * column[14] + .55555555555555555556f * column[15] + 2.2222222222222222222f * column[16] + -.28571428571428571429f * column[20] + .28571428571428571429f * column[21] + .71428571428571428571f * column[22] + 1.f * column[23] + 4.f * column[24] + -.14285714285714285714f * column[25] + .28571428571428571429f * column[26] + .57142857142857142857f * column[27] + -.71428571428571428571f * column[28] + -.47619047619047619048e-1f * column[29] + .23809523809523809524f * column[30] + .95238095238095238095f * column[31];
    temp[MATRIX(3, j + Quantities * 1, 20, Quantities * 3)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + .33333333333333333333f * column[6] + 1.5000000000000000000f * column[7] + 4.5000000000000000000f * column[8] + -.16666666666666666667f * column[10] + -.16666666666666666667f * column[11] + -.16666666666666666667f * column[12] + -.83333333333333333333f * column[13] + .61111111111111111111f * column[14] + .61111111111111111111f * column[15] + -1.2222222222222222222f * column[16] + 1.1666666666666666667f * column[17] + 3.5000000000000000000f * column[18] + -.14285714285714285714f * column[20] + -.14285714285714285714f * column[21] + -.14285714285714285714f * column[22] + -.14285714285714285714f * column[23] + .57142857142857142857f * column[24] + .28571428571428571429f * column[25] + .28571428571428571429f * column[26] + .28571428571428571429f * column[27] + 1.4285714285714285714f * column[28] + .61904761904761904762f * column[29] + .61904761904761904762f * column[30] + -1.2380952380952380952f * column[31] + .85714285714285714286f * column[32] + 2.5714285714285714286f * column[33];
    temp[MATRIX(4, j + Quantities * 1, 20, Quantities * 3)] = 5.f * column[10] + 7.f * column[11] + 3.7500000000000000000f * column[20] + 1.2500000000000000000f * column[21] + -5.f * column[22] + 1.2500000000000000000f * column[25] + 1.7500000000000000000f * column[26];
    temp[MATRIX(5, j + Quantities * 1, 20, Quantities * 3)] = -.60000000000000000000f * column[10] + 4.2000000000000000000f * column[11] + 8.4000000000000000000f * column[12] + -.75000000000000000000f * column[20] + 2.5500000000000000000f * column[21] + 2.4000000000000000000f * column[22] + -4.2000000000000000000f * column[23] + -.15000000000000000000f * column[25] + 1.0500000000000000000f * column[26] + 2.1000000000000000000f * column[27];
    temp[MATRIX(6, j + Quantities * 1, 20, Quantities * 3)] = .10000000000000000000f * column[10] + -.70000000000000000000f * column[11] + 2.1000000000000000000f * column[12] + 10.500000000000000000f * column[13] + .15000000000000000000f * column[20] + -.65000000000000000000f * column[21] + .50000000000000000000f * column[22] + 2.1000000000000000000f * column[23] + -2.1000000000000000000f * column[24] + .25000000000000000000e-1f * column[25] + -.17500000000000000000f * column[26] + .52500000000000000000f * column[27] + 2.6250000000000000000f * column[28];
    temp[MATRIX(7, j + Quantities * 1, 20, Quantities * 3)] = -.40000000000000000000f * column[10] + -.20000000000000000000f * column[11] + .60000000000000000000f * column[12] + 4.f * column[14] + 6.6666666666666666667f * column[15] + -.50000000000000000000f * column[20] + -.30000000000000000000f * column[21] + -.15000000000000000000f * column[22] + -1.0500000000000000000f * column[23] + 2.4000000000000000000f * column[25] + 1.2000000000000000000f * column[26] + -3.6000000000000000000f * column[27] + 2.2500000000000000000f * column[29] + 3.7500000000000000000f * column[30];
    temp[MATRIX(8, j + Quantities * 1, 20, Quantities * 3)] = .66666666666666666667e-1f * column[10] + -.13333333333333333333f * column[11] + -.26666666666666666667f * column[12] + .33333333333333333333f * column[13] + -.44444444444444444444f * column[14] + 2.2222222222222222222f * column[15] + 8.8888888888888888889f * column[16] + .10000000000000000000f * column[20] + -.10000000000000000000f * column[21] + -.25000000000000000000f * column[22] + -.35000000000000000000f * column[23] + -1.4000000000000000000f * column[24] + -.40000000000000000000f * column[25] + .80000000000000000000f * column[26] + 1.6000000000000000000f * column[27] + -2.f * column[28] + -.25000000000000000000f * column[29] + 1.2500000000000000000f * column[30] + 5.f * column[31];
    temp[MATRIX(9, j + Quantities * 1, 20, Quantities * 3)] = .33333333333333333333e-1f * column[10] + .33333333333333333333e-1f * column[11] + .33333333333333333333e-1f * column[12] + .16666666666666666667f * column[13] + -.35555555555555555556f * column[14] + -.35555555555555555556f * column[15] + .71111111111111111111f * column[16] + 1.8666666666666666667f * column[17] + 5.6000000000000000000f * column[18] + .50000000000000000000e-1f * column[20] + .50000000000000000000e-1f * column[21] + .50000000000000000000e-1f * column[22] + .50000000000000000000e-1f * column[23] + -.20000000000000000000f * column[24] + -.32500000000000000000f * column[25] + -.32500000000000000000f * column[26] + -.32500000000000000000f * column[27] + -1.6250000000000000000f * column[28] + .55000000000000000000f * column[29] + .55000000000000000000f * column[30] + -1.1000000000000000000f * column[31] + 1.8000000000000000000f * column[32] + 5.4000000000000000000f * column[33];
    temp[MATRIX(10, j + Quantities * 1, 20, Quantities * 3)] = 7.f * column[20] + 9.f * column[21];
    temp[MATRIX(11, j + Quantities * 1, 20, Quantities * 3)] = -.71428571428571428571f * column[20] + 6.4285714285714285714f * column[21] + 10.285714285714285714f * column[22];
    temp[MATRIX(12, j + Quantities * 1, 20, Quantities * 3)] = .14285714285714285714f * column[20] + -1.2857142857142857143f * column[21] + 5.1428571428571428571f * column[22] + 12.f * column[23];
    temp[MATRIX(13, j + Quantities * 1, 20, Quantities * 3)] = -.28571428571428571429e-1f * column[20] + .25714285714285714286f * column[21] + -1.0285714285714285714f * column[22] + 2.4000000000000000000f * column[23] + 14.400000000000000000f * column[24];
    temp[MATRIX(14, j + Quantities * 1, 20, Quantities * 3)] = -.53571428571428571429f * column[20] + -.17857142857142857143f * column[21] + .71428571428571428571f * column[22] + 6.2500000000000000000f * column[25] + 8.7500000000000000000f * column[26];
    temp[MATRIX(15, j + Quantities * 1, 20, Quantities * 3)] = .10714285714285714286f * column[20] + -.36428571428571428571f * column[21] + -.34285714285714285714f * column[22] + .60000000000000000000f * column[23] + -.75000000000000000000f * column[25] + 5.2500000000000000000f * column[26] + 10.500000000000000000f * column[27];
    temp[MATRIX(16, j + Quantities * 1, 20, Quantities * 3)] = -.21428571428571428571e-1f * column[20] + .92857142857142857143e-1f * column[21] + -.71428571428571428571e-1f * column[22] + -.30000000000000000000f * column[23] + .30000000000000000000f * column[24] + .12500000000000000000f * column[25] + -.87500000000000000000f * column[26] + 2.6250000000000000000f * column[27] + 13.125000000000000000f * column[28];
    temp[MATRIX(17, j + Quantities * 1, 20, Quantities * 3)] = .71428571428571428571e-1f * column[20] + .42857142857142857143e-1f * column[21] + .21428571428571428571e-1f * column[22] + .15000000000000000000f * column[23] + -.85714285714285714286f * column[25] + -.42857142857142857143f * column[26] + 1.2857142857142857143f * column[27] + 4.8214285714285714286f * column[29] + 8.0357142857142857143f * column[30];
    temp[MATRIX(18, j + Quantities * 1, 20, Quantities * 3)] = -.14285714285714285714e-1f * column[20] + .14285714285714285714e-1f * column[21] + .35714285714285714286e-1f * column[22] + .50000000000000000000e-1f * column[23] + .20000000000000000000f * column[24] + .14285714285714285714f * column[25] + -.28571428571428571429f * column[26] + -.57142857142857142857f * column[27] + .71428571428571428571f * column[28] + -.53571428571428571429f * column[29] + 2.6785714285714285714f * column[30] + 10.714285714285714286f * column[31];
    temp[MATRIX(19, j + Quantities * 1, 20, Quantities * 3)] = -.71428571428571428571e-2f * column[20] + -.71428571428571428571e-2f * column[21] + -.71428571428571428571e-2f * column[22] + -.71428571428571428571e-2f * column[23] + .28571428571428571429e-1f * column[24] + .89285714285714285714e-1f * column[25] + .89285714285714285714e-1f * column[26] + .89285714285714285714e-1f * column[27] + .44642857142857142857f * column[28] + -.53571428571428571429f * column[29] + -.53571428571428571429f * column[30] + 1.0714285714285714286f * column[31] + 2.1428571428571428571f * column[32] + 6.4285714285714285714f * column[33];
    temp[MATRIX(0, j + Quantities * 2, 20, Quantities * 3)] = 1.f * column[1] + 1.f * column[2] + 4.f * column[3] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + .50000000000000000000f * column[6] + .50000000000000000000f * column[7] + .50000000000000000000f * column[8] + -2.5000000000000000000f * column[9] + .30000000000000000000f * column[10] + .30000000000000000000f * column[11] + .30000000000000000000f * column[12] + .30000000000000000000f * column[13] + .30000000000000000000f * column[14] + .30000000000000000000f * column[15] + .30000000000000000000f * column[16] + .30000000000000000000f * column[17] + .30000000000000000000f * column[18] + 3.3000000000000000000f * column[19] + .20000000000000000000f * column[20] + .20000000000000000000f * column[21] + .20000000000000000000f * column[22] + .20000000000000000000f * column[23] + .20000000000000000000f * column[24] + .20000000000000000000f * column[25] + .20000000000000000000f * column[26] + .20000000000000000000f * column[27] + .20000000000000000000f * column[28] + .20000000000000000000f * column[29] + .20000000000000000000f * column[30] + .20000000000000000000f * column[31] + .20000000000000000000f * column[32] + .20000000000000000000f * column[33] + -2.8000000000000000000f * column[34];
    temp[MATRIX(1, j + Quantities * 2, 20, Quantities * 3)] = 3.f * column[4] + 1.f * column[5] + 6.f * column[7] + 2.f * column[10] + 1.f * column[11] + .33333333333333333333f * column[12] + 1.f * column[14] + .33333333333333333333f * column[15] + -4.6666666666666666667f * column[17] + 1.4285714285714285714f * column[20] + .85714285714285714286f * column[21] + .42857142857142857143f * column[22] + .14285714285714285714f * column[23] + .85714285714285714286f * column[25] + .42857142857142857143f * column[26] + .14285714285714285714f * column[27] + .42857142857142857143f * column[29] + .14285714285714285714f * column[30] + 5.1428571428571428571f * column[32];
    temp[MATRIX(2, j + Quantities * 2, 20, Quantities * 3)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 2.6666666666666666667f * column[6] + 6.f * column[8] + -.33333333333333333333f * column[10] + .66666666666666666667f * column[11] + 1.3333333333333333333f * column[12] + 1.6666666666666666667f * column[13] + -.11111111111111111111f * column[14] + .55555555555555555556f * column[15] + .88888888888888888889f * column[16] + -4.6666666666666666667f * column[18] + -.28571428571428571429f * column[20] + .28571428571428571429f * column[21] + .71428571428571428571f * column[22] + 1.f * column[23] + 1.1428571428571428571f * column[24] + -.14285714285714285714f * column[25] + .28571428571428571429f * column[26] + .57142857142857142857f * column[27] + .71428571428571428571f * column[28] + -.47619047619047619048e-1f * column[29] + .23809523809523809524f * column[30] + .38095238095238095238f * column[31] + 5.1428571428571428571f * column[33];
    temp[MATRIX(3, j + Quantities * 2, 20, Quantities * 3)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + -.16666666666666666667f * column[6] + 1.5000000000000000000f * column[7] + 1.5000000000000000000f * column[8] + 7.5000000000000000000f * column[9] + -.16666666666666666667f * column[10] + -.16666666666666666667f * column[11] + -.16666666666666666667f * column[12] + -.16666666666666666667f * column[13] + .61111111111111111111f * column[14] + .61111111111111111111f * column[15] + .61111111111111111111f * column[16] + 1.1666666666666666667f * column[17] + 1.1666666666666666667f * column[18] + -3.5000000000000000000f * column[19] + -.14285714285714285714f * column[20] + -.14285714285714285714f * column[21] + -.14285714285714285714f * column[22] + -.14285714285714285714f * column[23] + -.14285714285714285714f * column[24] + .28571428571428571429f * column[25] + .28571428571428571429f * column[26] + .28571428571428571429f * column[27] + .28571428571428571429f * column[28] + .61904761904761904762f * column[29] + .61904761904761904762f * column[30] + .61904761904761904762f * column[31] + .85714285714285714286f * column[32] + .85714285714285714286f * column[33] + 6.f * column[34];
    temp[MATRIX(4, j + Quantities * 2, 20, Quantities * 3)] = 5.f * column[10] + 1.f * column[11] + 8.f * column[14] + 3.7500000000000000000f * column[20] + 1.2500000000000000000f * column[21] + .25000000000000000000f * column[22] + 1.2500000000000000000f * column[25] + .25000000000000000000f * column[26] + -6.7500000000000000000f * column[29];
    temp[MATRIX(5, j + Quantities * 2, 20, Quantities * 3)] = -.60000000000000000000f * column[10] + 4.2000000000000000000f * column[11] + 2.4000000000000000000f * column[12] + 8.f * column[15] + -.75000000000000000000f * column[20] + 2.5500000000000000000f * column[21] + 2.4000000000000000000f * column[22] + 1.0500000000000000000f * column[23] + -.15000000000000000000f * column[25] + 1.0500000000000000000f * column[26] + .60000000000000000000f * column[27] + -6.7500000000000000000f * column[30];
    temp[MATRIX(6, j + Quantities * 2, 20, Quantities * 3)] = .10000000000000000000f * column[10] + -.70000000000000000000f * column[11] + 2.1000000000000000000f * column[12] + 4.5000000000000000000f * column[13] + 8.f * column[16] + .15000000000000000000f * column[20] + -.65000000000000000000f * column[21] + .50000000000000000000f * column[22] + 2.1000000000000000000f * column[23] + 3.1500000000000000000f * column[24] + .25000000000000000000e-1f * column[25] + -.17500000000000000000f * column[26] + .52500000000000000000f * column[27] + 1.1250000000000000000f * column[28] + -6.7500000000000000000f * column[31];
    temp[MATRIX(7, j + Quantities * 2, 20, Quantities * 3)] = -.40000000000000000000f * column[10] + -.20000000000000000000f * column[11] + -.66666666666666666667e-1f * column[12] + 4.f * column[14] + 1.3333333333333333333f * column[15] + 9.3333333333333333333f * column[17] + -.50000000000000000000f * column[20] + -.30000000000000000000f * column[21] + -.15000000000000000000f * column[22] + -.50000000000000000000e-1f * column[23] + 2.4000000000000000000f * column[25] + 1.2000000000000000000f * column[26] + .40000000000000000000f * column[27] + 2.2500000000000000000f * column[29] + .75000000000000000000f * column[30] + -6.f * column[32];
    temp[MATRIX(8, j + Quantities * 2, 20, Quantities * 3)] = .66666666666666666667e-1f * column[10] + -.13333333333333333333f * column[11] + -.26666666666666666667f * column[12] + -.33333333333333333333f * column[13] + -.44444444444444444444f * column[14] + 2.2222222222222222222f * column[15] + 3.5555555555555555556f * column[16] + 9.3333333333333333333f * column[18] + .10000000000000000000f * column[20] + -.10000000000000000000f * column[21] + -.25000000000000000000f * column[22] + -.35000000000000000000f * column[23] + -.40000000000000000000f * column[24] + -.40000000000000000000f * column[25] + .80000000000000000000f * column[26] + 1.6000000000000000000f * column[27] + 2.f * column[28] + -.25000000000000000000f * column[29] + 1.2500000000000000000f * column[30] + 2.f * column[31] + -6.f * column[33];
    temp[MATRIX(9, j + Quantities * 2, 20, Quantities * 3)] = .33333333333333333333e-1f * column[10] + .33333333333333333333e-1f * column[11] + .33333333333333333333e-1f * column[12] + .33333333333333333333e-1f * column[13] + -.35555555555555555556f * column[14] + -.35555555555555555556f * column[15] + -.35555555555555555556f * column[16] + 1.8666666666666666667f * column[17] + 1.8666666666666666667f * column[18] + 11.200000000000000000f * column[19] + .50000000000000000000e-1f * column[20] + .50000000000000000000e-1f * column[21] + .50000000000000000000e-1f * column[22] + .50000000000000000000e-1f * column[23] + .50000000000000000000e-1f * column[24] + -.32500000000000000000f * column[25] + -.32500000000000000000f * column[26] + -.32500000000000000000f * column[27] + -.32500000000000000000f * column[28] + .55000000000000000000f * column[29] + .55000000000000000000f * column[30] + .55000000000000000000f * column[31] + 1.8000000000000000000f * column[32] + 1.8000000000000000000f * column[33] + -4.2000000000000000000f * column[34];
    temp[MATRIX(10, j + Quantities * 2, 20, Quantities * 3)] = 7.f * column[20] + 1.f * column[21] + 10.f * column[25];
    temp[MATRIX(11, j + Quantities * 2, 20, Quantities * 3)] = -.71428571428571428571f * column[20] + 6.4285714285714285714f * column[21] + 2.2857142857142857143f * column[22] + 10.f * column[26];
    temp[MATRIX(12, j + Quantities * 2, 20, Quantities * 3)] = .14285714285714285714f * column[20] + -1.2857142857142857143f * column[21] + 5.1428571428571428571f * column[22] + 4.f * column[23] + 10.f * column[27];
    temp[MATRIX(13, j + Quantities * 2, 20, Quantities * 3)] = -.28571428571428571429e-1f * column[20] + .25714285714285714286f * column[21] + -1.0285714285714285714f * column[22] + 2.4000000000000000000f * column[23] + 6.4000000000000000000f * column[24] + 10.f * column[28];
    temp[MATRIX(14, j + Quantities * 2, 20, Quantities * 3)] = -.53571428571428571429f * column[20] + -.17857142857142857143f * column[21] + -.35714285714285714286e-1f * column[22] + 6.2500000000000000000f * column[25] + 1.2500000000000000000f * column[26] + 11.250000000000000000f * column[29];
    temp[MATRIX(15, j + Quantities * 2, 20, Quantities * 3)] = .10714285714285714286f * column[20] + -.36428571428571428571f * column[21] + -.34285714285714285714f * column[22] + -.15000000000000000000f * column[23] + -.75000000000000000000f * column[25] + 5.2500000000000000000f * column[26] + 3.f * column[27] + 11.250000000000000000f * column[30];
    temp[MATRIX(16, j + Quantities * 2, 20, Quantities * 3)] = -.21428571428571428571e-1f * column[20] + .92857142857142857143e-1f * column[21] + -.71428571428571428571e-1f * column[22] + -.30000000000000000000f * column[23] + -.45000000000000000000f * column[24] + .12500000000000000000f * column[25] + -.87500000000000000000f * column[26] + 2.6250000000000000000f * column[27] + 5.6250000000000000000f * column[28] + 11.250000000000000000f * column[31];
    temp[MATRIX(17, j + Quantities * 2, 20, Quantities * 3)] = .71428571428571428571e-1f * column[20] + .42857142857142857143e-1f * column[21] + .21428571428571428571e-1f * column[22] + .71428571428571428571e-2f * column[23] + -.85714285714285714286f * column[25] + -.42857142857142857143f * column[26] + -.14285714285714285714f * column[27] + 4.8214285714285714286f * column[29] + 1.6071428571428571429f * column[30] + 12.857142857142857143f * column[32];
    temp[MATRIX(18, j + Quantities * 2, 20, Quantities * 3)] = -.14285714285714285714e-1f * column[20] + .14285714285714285714e-1f * column[21] + .35714285714285714286e-1f * column[22] + .50000000000000000000e-1f * column[23] + .57142857142857142857e-1f * column[24] + .14285714285714285714f * column[25] + -.28571428571428571429f * column[26] + -.57142857142857142857f * column[27] + -.71428571428571428571f * column[28] + -.53571428571428571429f * column[29] + 2.6785714285714285714f * column[30] + 4.2857142857142857143f * column[31] + 12.857142857142857143f * column[33];
    temp[MATRIX(19, j + Quantities * 2, 20, Quantities * 3)] = -.71428571428571428571e-2f * column[20] + -.71428571428571428571e-2f * column[21] + -.71428571428571428571e-2f * column[22] + -.71428571428571428571e-2f * column[23] + -.71428571428571428571e-2f * column[24] + .89285714285714285714e-1f * column[25] + .89285714285714285714e-1f * column[26] + .89285714285714285714e-1f * column[27] + .89285714285714285714e-1f * column[28] + -.53571428571428571429f * column[29] + -.53571428571428571429f * column[30] + -.53571428571428571429f * column[31] + 2.1428571428571428571f * column[32] + 2.1428571428571428571f * column[33] + 15.f * column[34];
    }
}
__device__ __forceinline__ void dgkernel20(real* __restrict__ temp, const real* __restrict__ dQ) {
    #pragma unroll
for (int j = 0; j < Quantities; ++j) {
    real column[20];
    #pragma unroll
    for (int i = 0; i < 20; ++i) {
        column[i] = dQ[MATRIX(i, j, 20, Quantities)];
    }

    temp[MATRIX(0, j + Quantities * 0, 10, Quantities * 3)] = 2.f * column[1] + 1.f * column[5] + 1.f * column[7] + .60000000000000000000f * column[10] + .60000000000000000000f * column[12] + .60000000000000000000f * column[15] + .60000000000000000000f * column[17];
    temp[MATRIX(1, j + Quantities * 0, 10, Quantities * 3)] = 6.f * column[4] + 2.f * column[11] + 2.f * column[14];
    temp[MATRIX(2, j + Quantities * 0, 10, Quantities * 3)] = 3.3333333333333333333f * column[5] + -.66666666666666666667f * column[10] + 2.6666666666666666667f * column[12] + 1.1111111111111111111f * column[15];
    temp[MATRIX(3, j + Quantities * 0, 10, Quantities * 3)] = -.33333333333333333333f * column[5] + 3.f * column[7] + -.33333333333333333333f * column[10] + -.33333333333333333333f * column[12] + 1.2222222222222222222f * column[15] + 2.3333333333333333333f * column[17];
    temp[MATRIX(4, j + Quantities * 0, 10, Quantities * 3)] = 10.f * column[10];
    temp[MATRIX(5, j + Quantities * 0, 10, Quantities * 3)] = 8.4000000000000000000f * column[11];
    temp[MATRIX(6, j + Quantities * 0, 10, Quantities * 3)] = .20000000000000000000f * column[10] + 4.2000000000000000000f * column[12];
    temp[MATRIX(7, j + Quantities * 0, 10, Quantities * 3)] = -.40000000000000000000f * column[11] + 8.f * column[14];
    temp[MATRIX(8, j + Quantities * 0, 10, Quantities * 3)] = .13333333333333333333f * column[10] + -.53333333333333333333f * column[12] + 4.4444444444444444444f * column[15];
    temp[MATRIX(9, j + Quantities * 0, 10, Quantities * 3)] = .66666666666666666667e-1f * column[10] + .66666666666666666667e-1f * column[12] + -.71111111111111111111f * column[15] + 3.7333333333333333333f * column[17];
    temp[MATRIX(0, j + Quantities * 1, 10, Quantities * 3)] = 1.f * column[1] + 3.f * column[2] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + -1.f * column[6] + .50000000000000000000f * column[7] + 1.5000000000000000000f * column[8] + .30000000000000000000f * column[10] + .30000000000000000000f * column[11] + .30000000000000000000f * column[12] + 1.5000000000000000000f * column[13] + .30000000000000000000f * column[14] + .30000000000000000000f * column[15] + -.60000000000000000000f * column[16] + .30000000000000000000f * column[17] + .90000000000000000000f * column[18];
    temp[MATRIX(1, j + Quantities * 1, 10, Quantities * 3)] = 3.f * column[4] + 5.f * column[5] + 2.f * column[10] + 1.f * column[11] + -3.f * column[12] + 1.f * column[14] + 1.6666666666666666667f * column[15];
    temp[MATRIX(2, j + Quantities * 1, 10, Quantities * 3)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 6.6666666666666666667f * column[6] + -.33333333333333333333f * column[10] + .66666666666666666667f * column[11] + 1.3333333333333333333f * column[12] + -1.6666666666666666667f * column[13] + -.11111111111111111111f * column[14] + .55555555555555555556f * column[15] + 2.2222222222222222222f * column[16];
    temp[MATRIX(3, j + Quantities * 1, 10, Quantities * 3)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + .33333333333333333333f * column[6] + 1.5000000000000000000f * column[7] + 4.5000000000000000000f * column[8] + -.16666666666666666667f * column[10] + -.16666666666666666667f * column[11] + -.16666666666666666667f * column[12] + -.83333333333333333333f * column[13] + .61111111111111111111f * column[14] + .61111111111111111111f * column[15] + -1.2222222222222222222f * column[16] + 1.1666666666666666667f * column[17] + 3.5000000000000000000f * column[18];
    temp[MATRIX(4, j + Quantities * 1, 10, Quantities * 3)] = 5.f * column[10] + 7.f * column[11];
    temp[MATRIX(5, j + Quantities * 1, 10, Quantities * 3)] = -.60000000000000000000f * column[10] + 4.2000000000000000000f * column[11] + 8.4000000000000000000f * column[12];
    temp[MATRIX(6, j + Quantities * 1, 10, Quantities * 3)] = .10000000000000000000f * column[10] + -.70000000000000000000f * column[11] + 2.1000000000000000000f * column[12] + 10.500000000000000000f * column[13];
    temp[MATRIX(7, j + Quantities * 1, 10, Quantities * 3)] = -.40000000000000000000f * column[10] + -.20000000000000000000f * column[11] + .60000000000000000000f * column[12] + 4.f * column[14] + 6.6666666666666666667f * column[15];
    temp[MATRIX(8, j + Quantities * 1, 10, Quantities * 3)] = .66666666666666666667e-1f * column[10] + -.13333333333333333333f * column[11] + -.26666666666666666667f * column[12] + .33333333333333333333f * column[13] + -.44444444444444444444f * column[14] + 2.2222222222222222222f * column[15] + 8.8888888888888888889f * column[16];
    temp[MATRIX(9, j + Quantities * 1, 10, Quantities * 3)] = .33333333333333333333e-1f * column[10] + .33333333333333333333e-1f * column[11] + .33333333333333333333e-1f * column[12] + .16666666666666666667f * column[13] + -.35555555555555555556f * column[14] + -.35555555555555555556f * column[15] + .71111111111111111111f * column[16] + 1.8666666666666666667f * column[17] + 5.6000000000000000000f * column[18];
    temp[MATRIX(0, j + Quantities * 2, 10, Quantities * 3)] = 1.f * column[1] + 1.f * column[2] + 4.f * column[3] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + .50000000000000000000f * column[6] + .50000000000000000000f * column[7] + .50000000000000000000f * column[8] + -2.5000000000000000000f * column[9] + .30000000000000000000f * column[10] + .30000000000000000000f * column[11] + .30000000000000000000f * column[12] + .30000000000000000000f * column[13] + .30000000000000000000f * column[14] + .30000000000000000000f * column[15] + .30000000000000000000f * column[16] + .30000000000000000000f * column[17] + .30000000000000000000f * column[18] + 3.3000000000000000000f * column[19];
    temp[MATRIX(1, j + Quantities * 2, 10, Quantities * 3)] = 3.f * column[4] + 1.f * column[5] + 6.f * column[7] + 2.f * column[10] + 1.f * column[11] + .33333333333333333333f * column[12] + 1.f * column[14] + .33333333333333333333f * column[15] + -4.6666666666666666667f * column[17];
    temp[MATRIX(2, j + Quantities * 2, 10, Quantities * 3)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 2.6666666666666666667f * column[6] + 6.f * column[8] + -.33333333333333333333f * column[10] + .66666666666666666667f * column[11] + 1.3333333333333333333f * column[12] + 1.6666666666666666667f * column[13] + -.11111111111111111111f * column[14] + .55555555555555555556f * column[15] + .88888888888888888889f * column[16] + -4.6666666666666666667f * column[18];
    temp[MATRIX(3, j + Quantities * 2, 10, Quantities * 3)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + -.16666666666666666667f * column[6] + 1.5000000000000000000f * column[7] + 1.5000000000000000000f * column[8] + 7.5000000000000000000f * column[9] + -.16666666666666666667f * column[10] + -.16666666666666666667f * column[11] + -.16666666666666666667f * column[12] + -.16666666666666666667f * column[13] + .61111111111111111111f * column[14] + .61111111111111111111f * column[15] + .61111111111111111111f * column[16] + 1.1666666666666666667f * column[17] + 1.1666666666666666667f * column[18] + -3.5000000000000000000f * column[19];
    temp[MATRIX(4, j + Quantities * 2, 10, Quantities * 3)] = 5.f * column[10] + 1.f * column[11] + 8.f * column[14];
    temp[MATRIX(5, j + Quantities * 2, 10, Quantities * 3)] = -.60000000000000000000f * column[10] + 4.2000000000000000000f * column[11] + 2.4000000000000000000f * column[12] + 8.f * column[15];
    temp[MATRIX(6, j + Quantities * 2, 10, Quantities * 3)] = .10000000000000000000f * column[10] + -.70000000000000000000f * column[11] + 2.1000000000000000000f * column[12] + 4.5000000000000000000f * column[13] + 8.f * column[16];
    temp[MATRIX(7, j + Quantities * 2, 10, Quantities * 3)] = -.40000000000000000000f * column[10] + -.20000000000000000000f * column[11] + -.66666666666666666667e-1f * column[12] + 4.f * column[14] + 1.3333333333333333333f * column[15] + 9.3333333333333333333f * column[17];
    temp[MATRIX(8, j + Quantities * 2, 10, Quantities * 3)] = .66666666666666666667e-1f * column[10] + -.13333333333333333333f * column[11] + -.26666666666666666667f * column[12] + -.33333333333333333333f * column[13] + -.44444444444444444444f * column[14] + 2.2222222222222222222f * column[15] + 3.5555555555555555556f * column[16] + 9.3333333333333333333f * column[18];
    temp[MATRIX(9, j + Quantities * 2, 10, Quantities * 3)] = .33333333333333333333e-1f * column[10] + .33333333333333333333e-1f * column[11] + .33333333333333333333e-1f * column[12] + .33333333333333333333e-1f * column[13] + -.35555555555555555556f * column[14] + -.35555555555555555556f * column[15] + -.35555555555555555556f * column[16] + 1.8666666666666666667f * column[17] + 1.8666666666666666667f * column[18] + 11.200000000000000000f * column[19];
    }
}
__device__ __forceinline__ void dgkernel10(real* __restrict__ temp, const real* __restrict__ dQ) {
    #pragma unroll
for (int j = 0; j < Quantities; ++j) {
    real column[10];
    #pragma unroll
    for (int i = 0; i < 10; ++i) {
        column[i] = dQ[MATRIX(i, j, 10, Quantities)];
    }

    temp[0 * Quantities * 3 + (j + Quantities * 0)] = 2.f * column[1] + 1.f * column[5] + 1.f * column[7];
    temp[1 * Quantities * 3 + (j + Quantities * 0)] = 6.f * column[4];
    temp[2 * Quantities * 3 + (j + Quantities * 0)] = 3.3333333333333333333f * column[5];
    temp[3 * Quantities * 3 + (j + Quantities * 0)] = -.33333333333333333333f * column[5] + 3.f * column[7];
    temp[0 * Quantities * 3 + (j + Quantities * 1)] = 1.f * column[1] + 3.f * column[2] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + -1.f * column[6] + .50000000000000000000f * column[7] + 1.5000000000000000000f * column[8];
    temp[1 * Quantities * 3 + (j + Quantities * 1)] = 3.f * column[4] + 5.f * column[5];
    temp[2 * Quantities * 3 + (j + Quantities * 1)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 6.6666666666666666667f * column[6];
    temp[3 * Quantities * 3 + (j + Quantities * 1)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + .33333333333333333333f * column[6] + 1.5000000000000000000f * column[7] + 4.5000000000000000000f * column[8];
    temp[0 * Quantities * 3 + (j + Quantities * 2)] = 1.f * column[1] + 1.f * column[2] + 4.f * column[3] + .50000000000000000000f * column[4] + .50000000000000000000f * column[5] + .50000000000000000000f * column[6] + .50000000000000000000f * column[7] + .50000000000000000000f * column[8] + -2.5000000000000000000f * column[9];
    temp[1 * Quantities * 3 + (j + Quantities * 2)] = 3.f * column[4] + 1.f * column[5] + 6.f * column[7];
    temp[2 * Quantities * 3 + (j + Quantities * 2)] = -.33333333333333333333f * column[4] + 1.6666666666666666667f * column[5] + 2.6666666666666666667f * column[6] + 6.f * column[8];
    temp[3 * Quantities * 3 + (j + Quantities * 2)] = -.16666666666666666667f * column[4] + -.16666666666666666667f * column[5] + -.16666666666666666667f * column[6] + 1.5000000000000000000f * column[7] + 1.5000000000000000000f * column[8] + 7.5000000000000000000f * column[9];
    }
}
__device__ __forceinline__ void dgkernel4(real* __restrict__ temp, const real* __restrict__ dQ) {
    #pragma unroll
for (int j = 0; j < Quantities; ++j) {
    real column[4];
    #pragma unroll
    for (int i = 0; i < 4; ++i) {
        column[i] = dQ[MATRIX(i, j, 4, Quantities)];
    }

    temp[0 * Quantities * 3 + (j + Quantities * 0)] = 2.f * column[1];
    temp[0 * Quantities * 3 + (j + Quantities * 1)] = 1.f * column[1] + 3.f * column[2];
    temp[0 * Quantities * 3 + (j + Quantities * 2)] = 1.f * column[1] + 1.f * column[2] + 4.f * column[3];
    }
}



template<size_t N, size_t N1, size_t Nmax>
__device__ __forceinline__ void dgkernelPart1(real* __restrict__ temp, const real* __restrict__ dQ, const real* __restrict__ derivativeX,
const real* __restrict__ derivativeY, const real* __restrict__ derivativeZ) {
    constexpr size_t TileK = 12;
    for (int j = 0; j < Quantities; ++j) {
        float column[N1];
        #pragma unroll
        for (int i = 0; i < N1; ++i) {
            column[i] = dQ[MATRIX(i, j, N1, Quantities)];
        }
        for (int i = 0; i < N; ++i) {
            float outX = 0, outY = 0, outZ = 0;
            constexpr int kr = N1 - (N1 % TileK);
            /*
            for (int k = 0; k < kr; k += TileK) {
                #pragma unroll
                for (int kk = 0; kk < TileK; ++kk) {
                    outX += derivativeX[i * Nmax + k + kk] * column[k + kk];
                }
            }
            #pragma unroll
            for (int kk = 0; kk < N1 % TileK; ++kk) {
                outX += derivativeX[i * Nmax + kr + kk] * column[kr + kk];
            }
            for (int k = 0; k < kr; k += TileK) {
                #pragma unroll
                for (int kk = 0; kk < TileK; ++kk) {
                    outY += derivativeY[i * Nmax + k + kk] * column[k + kk];
                }
            }
            #pragma unroll
            for (int kk = 0; kk < N1 % TileK; ++kk) {
                outY += derivativeY[i * Nmax + kr + kk] * column[kr + kk];
            }
            for (int k = 0; k < kr; k += TileK) {
                #pragma unroll
                for (int kk = 0; kk < TileK; ++kk) {
                    outZ += derivativeZ[i * Nmax + k + kk] * column[k + kk];
                }
            }
            #pragma unroll
            for (int kk = 0; kk < N1 % TileK; ++kk) {
                outZ += derivativeZ[i * Nmax + kr + kk] * column[kr + kk];
            }
            */
            for (int k = 0; k < kr; k += TileK) {
                /*#pragma unroll
                for (int kk = 0; kk < TileK; ++kk) {
                    outX += derivativeX[i * Nmax + k + kk] * column[k + kk];
                    outY += derivativeY[i * Nmax + k + kk] * column[k + kk];
                    outZ += derivativeZ[i * Nmax + k + kk] * column[k + kk];
                }*/
                #pragma unroll
                for (int kk = 0; kk < TileK; kk += 4) {
                    float4 dx = *(float4*)&derivativeX[i * Nmax + k + kk];
                    float4 dy = *(float4*)&derivativeY[i * Nmax + k + kk];
                    float4 dz = *(float4*)&derivativeZ[i * Nmax + k + kk];
                    outX += dx.x * column[k + kk + 0];
                    outY += dy.x * column[k + kk + 0];
                    outZ += dx.x * column[k + kk + 0];
                    outX += dx.y * column[k + kk + 1];
                    outY += dy.y * column[k + kk + 1];
                    outZ += dx.y * column[k + kk + 1];
                    outX += dx.z * column[k + kk + 2];
                    outY += dy.z * column[k + kk + 2];
                    outZ += dx.z * column[k + kk + 2];
                    outX += dx.w * column[k + kk + 3];
                    outY += dy.w * column[k + kk + 3];
                    outZ += dx.z * column[k + kk + 3];
                }
            }
            #pragma unroll
            for (int kk = 0; kk < N1 % TileK; ++kk) {
                outX += derivativeX[i * Nmax + kr + kk] * column[kr + kk];
                outY += derivativeY[i * Nmax + kr + kk] * column[kr + kk];
                outZ += derivativeZ[i * Nmax + kr + kk] * column[kr + kk];
            }
            temp[MATRIX(i, j, N, Quantities * 3)] = outX;
            temp[MATRIX(i, j + Quantities, N, Quantities * 3)] = outY;
            temp[MATRIX(i, j + Quantities * 2, N, Quantities * 3)] = outZ;
        }
    }
}

template<size_t N, size_t Nmax>
__device__ __forceinline__ void dgkernelPart2(real scale, real* __restrict__ output, real* __restrict__ acc, const real* __restrict__ temp, const real* __restrict__ coordinates,
    real lambda, real mu, real rhoD) {
        real l2mu = lambda + 2 * mu;
        #pragma unroll
    for (int i = 0; i < N; ++i) {
        real outRow[Quantities] = {0};
        real accRow[Quantities] = {0};
        #pragma unroll
        for (int j = 0; j < Quantities; ++j) {
            accRow[j] = acc[MATRIX(i,j, Nmax, Quantities)];
        }
        #pragma unroll
        for (int d = 0; d < 3; ++d) {
            real inRow[Quantities];
            #pragma unroll
            for (int j = 0; j < Quantities; ++j) {
                if constexpr (N >= 10) {
                    inRow[j] = temp[MATRIX(i, j * Quantities * d, N, Quantities * 3)];
                }
                else {
                    inRow[j] = temp[i * Quantities * 3 + (j + Quantities * d)];
                }
            }
            real n1 = coordinates[VECTOR(0 + d * 3, 9)];
            real n2 = coordinates[VECTOR(1 + d * 3, 9)];
            real n3 = coordinates[VECTOR(2 + d * 3, 9)];
            outRow[0] += -n1 * l2mu * inRow[0] - lambda * (n2 * inRow[1] + n3 * inRow[2]);
            outRow[1] += -n2 * l2mu * inRow[1] - lambda * (n1 * inRow[0] + n3 * inRow[2]);
            outRow[2] += -n3 * l2mu * inRow[2] - lambda * (n1 * inRow[0] + n2 * inRow[1]);
            outRow[3] += -mu * (n2 * inRow[0] + n1 * inRow[1]);
            outRow[4] += -mu * (n3 * inRow[1] + n2 * inRow[2]);
            outRow[5] += -mu * (n1 * inRow[2] + n3 * inRow[0]);
            outRow[6] += -rhoD * (n1 * inRow[0] + n2 * inRow[3] + n3 * inRow[5]);
            outRow[7] += -rhoD * (n1 * inRow[3] + n2 * inRow[1] + n3 * inRow[4]);
            outRow[8] += -rhoD * (n1 * inRow[5] + n2 * inRow[4] + n3 * inRow[2]);
        }
        #pragma unroll
        for (int j = 0; j < Quantities; ++j) {
            output[MATRIX(i,j, N, Quantities)] = outRow[j];
            acc[MATRIX(i,j, Nmax, Quantities)] = scale * accRow[j] + outRow[j];
        }
    }
}

template<size_t Nmax>
__device__ __forceinline__ void dgkernelInit(real scale, real* __restrict__ acc, const real* __restrict__ dofs) {
    // TODO: move into first kernel
        #pragma unroll
        for (int i = 0; i < Nmax; ++i) {
            for (int j = 0; j < Quantities; ++j) {
                acc[MATRIX(i,j, Nmax, Quantities)] = scale * dofs[MATRIX(i,j, Nmax, Quantities)];
            }
        }
    }

__global__ __launch_bounds__(AderMultiple*Blocksize) void dgkernelFull(std::size_t blockCount, real scale, real* __restrict__ I, const real* __restrict__ Q, real* __restrict__ dQ6, real* __restrict__ dQ5, real* __restrict__ dQ4, real* __restrict__ dQ3, real* __restrict__ dQ2, real* __restrict__ dQ1, const real* __restrict__ stardata, const real* __restrict__ coordinates, real* __restrict__ temp) {
    if (threadIdx.y + AderMultiple * blockIdx.x < blockCount) {
        real rhoD = stardata[VECTOR(0, 3)];
        real mu = stardata[VECTOR(1, 3)];
        real lambda = stardata[VECTOR(2, 3)];
        /*dgkernelPart1<56, 35, 56>(temp, dQ6, derivative6X, derivative6Y, derivative6Z);
        dgkernelPart2<35, 56>(dQ5, I,  temp, coordinates, lambda, mu, rhoD);
        dgkernelPart1<35, 20, 56>(temp, dQ5, derivative6X, derivative6Y, derivative6Z);
        dgkernelPart2<20, 56>(dQ4, I,  temp, coordinates, lambda, mu, rhoD);
        dgkernelPart1<20, 10, 56>(temp, dQ4, derivative6X, derivative6Y, derivative6Z);
        dgkernelPart2<10, 56>(dQ3, I,  temp, coordinates, lambda, mu, rhoD);
        dgkernelPart1<10, 4, 56>(temp, dQ3, derivative6X, derivative6Y, derivative6Z);
        dgkernelPart2<4, 56>(dQ2, I,  temp, coordinates, lambda, mu, rhoD);
        dgkernelPart1<4, 1, 56>(temp, dQ2, derivative6X, derivative6Y, derivative6Z);
        dgkernelPart2<1, 56>(dQ1, I,  temp, coordinates, lambda, mu, rhoD);*/
        real tempreg[4 * Quantities * 3];
        real coeff = scale;
        dgkernelInit<56>(coeff, I, Q);
        dgkernel56(temp, Q);
        dgkernelPart2<35, 56>(coeff, dQ5, I,  temp, coordinates, lambda, mu, rhoD);
        dgkernel35(temp, dQ5);
        dgkernelPart2<20, 56>(coeff, dQ4, I,  temp, coordinates, lambda, mu, rhoD);
        dgkernel20(temp, dQ4);
        dgkernelPart2<10, 56>(coeff, dQ3, I,  temp, coordinates, lambda, mu, rhoD);
        dgkernel10(tempreg, dQ3);
        dgkernelPart2<4, 56>(coeff, dQ2, I,  tempreg, coordinates, lambda, mu, rhoD);
        dgkernel4(tempreg, dQ2);
        dgkernelPart2<1, 56>(coeff, dQ1, I,  tempreg, coordinates, lambda, mu, rhoD);
    }
}

__global__ __launch_bounds__(AderMultiple*Blocksize) void kernelXTiled3x32(std::size_t blockCount, real scale, real* __restrict__ I, const real* __restrict__ Q, real* __restrict__ dQ6, real* __restrict__ dQ5, real* __restrict__ dQ4, real* __restrict__ dQ3, real* __restrict__ dQ2, real* __restrict__ dQ1, const real* __restrict__ stardata, const real* __restrict__ coordinates) {
    if (threadIdx.y + AderMultiple * blockIdx.x < blockCount) {
      real rhoD = stardata[VECTOR(0, 3)];
      real mu = stardata[VECTOR(1, 3)];
      real lambda = stardata[VECTOR(2, 3)];
      
      real scale1 = scale;
      real scale2 = 1;
      // TODO: zero init

      copykernel<6, 32>(dQ6, Q);

      // 56x9 -> 36x9
      dgkernel2b<12, 12, 14, 5, true>(dQ5, Q, coordinates, derivative6X, derivative6Y, derivative6Z, lambda, mu, rhoD, 0);
      sumkernel<5, 12, true>(I, dQ5, scale1);
      // sumkernel<5, true>(I, dQ5, scale1);
      scale1 *= scale;
      scale2 *= scale;

      // 36x9 -> 20x9
      dgkernel2b<10, 10, 12, 4, true>(dQ4, dQ5, coordinates, derivative6X, derivative6Y, derivative6Z, lambda, mu, rhoD, 0);
      sumkernel<4, 10, false>(I, dQ4, scale1);
      scale1 *= scale / 2;
      scale2 *= scale;

      // 20x9 -> 10x9
      dgkernel<10, 10, 10, 3, true>(dQ3, dQ3, coordinates, derivative6X, lambda, mu, rhoD, 0);
      dgkernel<10, 10, 10, 3, false>(dQ3, dQ4, coordinates, derivative6Y, lambda, mu, rhoD, 3);
      dgkernel<10, 10, 10, 3, false>(dQ3, dQ4, coordinates, derivative6Z, lambda, mu, rhoD, 6);
      sumkernel<3, 10, false>(I, dQ3, scale1);
      scale1 *= scale / 3;
      scale2 *= scale / 2;

      // 10x9 -> 4x9
      dgkernel<4, 4, 10, 2, true>(dQ2, dQ3, coordinates, derivative6X, lambda, mu, rhoD, 0);
      dgkernel<4, 4, 10, 2, false>(dQ2, dQ3, coordinates, derivative6Y, lambda, mu, rhoD, 3);
      dgkernel<4, 4, 10, 2, false>(dQ2, dQ3, coordinates, derivative6Z, lambda, mu, rhoD, 6);
      sumkernel<2, 4, false>(I, dQ2, scale1);
      scale1 *= scale / 4;
      scale2 *= scale / 3;

      // 4x9 -> 1x9
      dgkernel<1, 1, 4, 1, true>(dQ1, dQ2, coordinates, derivative6X, lambda, mu, rhoD, 0);
      dgkernel<1, 1, 4, 1, false>(dQ1, dQ2, coordinates, derivative6Y, lambda, mu, rhoD, 3);
      dgkernel<1, 1, 4, 1, false>(dQ1, dQ2, coordinates, derivative6Z, lambda, mu, rhoD, 6);
      sumkernel<1, 1, false>(I, dQ1, scale1);
    }
}
void aderLauncher(std::size_t count, real timestep, const real* dofs, real* buffers, real* derivatives, const real* stardata, const real* coordinates, real* temp, void* stream) {
  std::size_t blocks = (count + Blocksize - 1) / Blocksize;
  std::size_t launchSize = (blocks + AderMultiple - 1) / AderMultiple;

  dim3 grid(launchSize, 1, 1);
  dim3 block(Blocksize, AderMultiple, 1);
  std::vector<std::size_t> dataOffsets(6);
  dataOffsets[0] = 0;
  dataOffsets[1] = dataOffsets[0] + Functions<6> * Blocksize * blocks;
  dataOffsets[2] = dataOffsets[1] + Functions<5> * Blocksize * blocks;
  dataOffsets[3] = dataOffsets[2] + Functions<4> * Blocksize * blocks;
  dataOffsets[4] = dataOffsets[3] + Functions<3> * Blocksize * blocks;
  dataOffsets[5] = dataOffsets[4] + Functions<2> * Blocksize * blocks;
  cudaStream_t streamObject = reinterpret_cast<cudaStream_t>(stream);
  cudaMemcpyAsync(derivatives, dofs, sizeof(real) * blocks * Quantities * Functions<6>, cudaMemcpyDeviceToDevice, streamObject);
  dgkernelFull <<<grid, block, 0, streamObject>>> (blocks, timestep, buffers, dofs, derivatives, derivatives + dataOffsets[1], derivatives + dataOffsets[2], derivatives + dataOffsets[3], derivatives + dataOffsets[4], derivatives + dataOffsets[5], stardata, coordinates, temp);
  CHECK_ERR;
}

template<typename T>
__device__ __forceinline__ T& iacc1(T* __restrict__ data, size_t size, size_t inblock, size_t block) {
    return data[block * size + inblock + threadIdx.x];
}

constexpr size_t InterleaveMultiple = 3;

template<typename T>
__global__ __launch_bounds__(Blocksize * InterleaveMultiple) void interleave(const T** source, T* target, size_t realdim, size_t paddeddim, size_t size, size_t count) {
    __shared__ T swap[Blocksize * (Blocksize * InterleaveMultiple + 1)];
    const size_t block = blockIdx.x * InterleaveMultiple + threadIdx.y;
    if (block < (count + Blocksize - 1) / Blocksize) {
      using TPtr = const T*;
      TPtr sourcePtrs[Blocksize];
      bool notZero = false;
    #pragma unroll
      for (size_t j = 0; j < Blocksize; ++j) {
        sourcePtrs[j] = source[block * Blocksize + j];
        notZero |= sourcePtrs[j] != nullptr;
      }
      if (notZero) {
        for (size_t j = 0; j < size; j += Blocksize) {
          #pragma unroll
            for (size_t i = 0; i < Blocksize; ++i) {
              if (block * Blocksize + i < count && j + threadIdx.x < size && sourcePtrs[i] != nullptr) {
                swap[i * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + threadIdx.x] = sourcePtrs[i][j + threadIdx.x];
              }
              else {
                swap[i * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + threadIdx.x] = 0;
              }
            }
            __syncwarp();
          #pragma unroll
            for (size_t i = 0; i < Blocksize; ++i) {
              size_t index = i + j;
              if (i + j < size && index % paddeddim < realdim) {
                size_t realindex = (index / paddeddim) * realdim + (index % paddeddim);
                target[size * Blocksize * block + realindex * Blocksize + threadIdx.x] = swap[threadIdx.x * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + i];
              }
            }
        }
      }
    }
}

template<typename T>
__global__ __launch_bounds__(Blocksize * InterleaveMultiple) void deinterleave(const T* source, T** target, size_t realdim, size_t paddeddim, size_t size, size_t count) {
    __shared__ T swap[Blocksize * (Blocksize * InterleaveMultiple + 1)];
    const size_t block = blockIdx.x * InterleaveMultiple + threadIdx.y;
    if (block < (count + Blocksize - 1) / Blocksize) {
      using TPtr = T*;
      TPtr targetPtrs[Blocksize];
      bool notZero = false;
    #pragma unroll
      for (size_t j = 0; j < Blocksize; ++j) {
        targetPtrs[j] = target[block * Blocksize + j];
        notZero |= targetPtrs[j] != nullptr;
      }
      if (notZero) {
        for (size_t j = 0; j < size; j += Blocksize) {
          #pragma unroll
            for (size_t i = 0; i < Blocksize; ++i) {
              size_t index = i + j;
              if (index < size && index % paddeddim < realdim) {
                size_t realindex = (index / paddeddim) * realdim + (index % paddeddim);
                swap[i * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + threadIdx.x] = source[size * Blocksize * block + realindex * Blocksize + threadIdx.x];
              }
              else {
                swap[i * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + threadIdx.x] = 0;
              }
            }
            __syncwarp();
          #pragma unroll
            for (size_t i = 0; i < Blocksize; ++i) {
              if (block * Blocksize + i < count && j + threadIdx.x < size && targetPtrs[i] != nullptr) {
                targetPtrs[i][j + threadIdx.x] = swap[threadIdx.x * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + i];
              }
            }
        }
      }
    }
}

void interleaveLauncher(std::size_t count, std::size_t size, std::size_t realdim, std::size_t paddeddim, const real** indata, real* outdata, void* stream) {
  cudaStream_t streamObject = reinterpret_cast<cudaStream_t>(stream);
  dim3 grid((count + Blocksize * InterleaveMultiple - 1) / (Blocksize * InterleaveMultiple), 1, 1);
  dim3 block(Blocksize, InterleaveMultiple, 1);
  interleave<<<grid, block, 0, streamObject>>>(indata, outdata, realdim, paddeddim, size, count);
  CHECK_ERR;
}
void deinterleaveLauncher(std::size_t count, std::size_t size, std::size_t realdim, std::size_t paddeddim, const real* indata, real** outdata, void* stream) {
  cudaStream_t streamObject = reinterpret_cast<cudaStream_t>(stream);
  dim3 grid((count + Blocksize * InterleaveMultiple - 1) / (Blocksize * InterleaveMultiple), 1, 1);
  dim3 block(Blocksize, InterleaveMultiple, 1);
  deinterleave<<<grid, block, 0, streamObject>>>(indata, outdata, realdim, paddeddim, size, count);
  CHECK_ERR;
}

void setConstantData(const real* startptr) {
  cudaMemcpyToSymbol(derivative6X, startptr, sizeof(derivative6X));
  startptr += sizeof(derivative6X) / sizeof(real);
  cudaMemcpyToSymbol(derivative6Y, startptr, sizeof(derivative6Y));
  startptr += sizeof(derivative6Y) / sizeof(real);
  cudaMemcpyToSymbol(derivative6Z, startptr, sizeof(derivative6Z));
}

}

namespace seissol::kernels::local_flux::aux::details {

__global__ void kernelFreeSurfaceGravity(
  real** dofsFaceBoundaryNodalPtrs,
  real** displacementDataPtrs,
  double* rhos,
  double g,
  size_t numElements) {

  const int tid = threadIdx.x;
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    const double rho = rhos[elementId];
    real* elementBoundaryDofs = dofsFaceBoundaryNodalPtrs[elementId];
    real* elementDisplacement = displacementDataPtrs[elementId];

    constexpr auto numNodes = seissol::nodal::tensor::nodes2D::Shape[0];
    if (tid < numNodes) {
      constexpr auto ldINodal = yateto::leadDim<seissol::init::INodal>();

      const auto pressureAtBnd = static_cast<real>(-1.0) * rho * g * elementDisplacement[tid];

      #pragma unroll
      for (int component{0}; component < 3; ++component) {
        elementBoundaryDofs[tid + component * ldINodal] =
          2.0 * pressureAtBnd - elementBoundaryDofs[tid + component * ldINodal];
      }
    }
  }
}

void launchFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                              real** displacementDataPtrs,
                              double* rhos,
                              double g,
                              size_t numElements,
                              void* deviceStream) {
  dim3 block(yateto::leadDim<seissol::nodal::init::nodes2D>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelFreeSurfaceGravity<<<grid, block, 0, stream>>>(dofsFaceBoundaryNodalPtrs,
                                                       displacementDataPtrs,
                                                       rhos,
                                                       g,
                                                       numElements);
}


__global__ void  kernelEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                                    real** easiBoundaryMapPtrs,
                                    real** easiBoundaryConstantPtrs,
                                    size_t numElements) {

  const int tid = threadIdx.x;
  const int elementId = blockIdx.x;

  constexpr auto ldINodalDim = yateto::leadDim<seissol::init::INodal>();
  constexpr auto iNodalDim0 = seissol::tensor::INodal::Shape[0];
  constexpr auto iNodalDim1 = seissol::tensor::INodal::Shape[1];
  __shared__ __align__(8) real resultTerm[iNodalDim1][iNodalDim0];

  constexpr auto ldConstantDim = yateto::leadDim<seissol::init::easiBoundaryConstant>();
  constexpr auto constantDim0 = seissol::tensor::easiBoundaryConstant::Shape[0];
  constexpr auto constantDim1 = seissol::tensor::easiBoundaryConstant::Shape[1];
  __shared__ __align__(8) real rightTerm[iNodalDim1][ldConstantDim];

  constexpr auto ldMapDim = yateto::leadDim<seissol::init::easiBoundaryMap>();
  constexpr auto mapDim0 = seissol::tensor::easiBoundaryMap::Shape[0];
  constexpr auto mapDim1 = seissol::tensor::easiBoundaryMap::Shape[1];
  constexpr auto mapDim2 = seissol::tensor::easiBoundaryMap::Shape[2];
  __shared__ __align__(8) real leftTerm[mapDim0][mapDim2];

  static_assert(iNodalDim1 == constantDim0, "supposed to be equal");
  static_assert(iNodalDim1 == mapDim0, "supposed to be equal");


  if (elementId < numElements) {
    real* dofsFaceBoundaryNodal = dofsFaceBoundaryNodalPtrs[elementId];
    real* easiBoundaryMap = easiBoundaryMapPtrs[elementId];
    auto easiBoundaryConstant = easiBoundaryConstantPtrs[elementId];

    for (int i = tid; i < (ldConstantDim * constantDim1); i += blockDim.x) {
      const auto b = i % ldConstantDim;
      const auto l = i / ldConstantDim;
      rightTerm[b][l] = easiBoundaryConstant[i];
    }
    __syncthreads();

    for (int i = 0; i < iNodalDim1; ++i) {
      if (tid < iNodalDim0) resultTerm[i][tid] = 0.0;
    }
    __syncthreads();

    for (int b = 0; b < mapDim1; ++b) {
      for (int l = 0; l < mapDim2; ++l) {
        if (tid < mapDim0) {
          leftTerm[tid][l] = easiBoundaryMap[tid + ldMapDim * (b + l * mapDim1)];
        }
      }
      __syncthreads();

      if (tid < mapDim2) {
        const real col = dofsFaceBoundaryNodal[tid + b * ldINodalDim];
        for (int a = 0; a < mapDim0; ++a) {
          resultTerm[a][tid] += leftTerm[a][tid] * col;
        }
      }
      __syncthreads();
    }

    if (tid < iNodalDim0) {
      for (int a = 0; a < iNodalDim1; ++a) {
        dofsFaceBoundaryNodal[tid + a * ldINodalDim] = resultTerm[a][tid] + rightTerm[a][tid];
      }
    }
  }
}


void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream) {

  dim3 block(yateto::leadDim<seissol::init::INodal>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelEasiBoundary<<<grid, block, 0, stream>>>(dofsFaceBoundaryNodalPtrs,
                                                 easiBoundaryMapPtrs,
                                                 easiBoundaryConstantPtrs,
                                                 numElements);
}
} // seissol::kernels::local_flux::aux::details



namespace seissol::kernels::time::aux {
__global__ void kernelextractRotationMatrices(real** displacementToFaceNormalPtrs,
                                              real** displacementToGlobalDataPtrs,
                                              real** TPtrs,
                                              real** TinvPtrs,
                                              size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    real* displacementToFaceNormal = displacementToFaceNormalPtrs[elementId];
    real* displacementToGlobalData = displacementToGlobalDataPtrs[elementId];
    auto* T = TPtrs[elementId];
    auto* Tinv = TinvPtrs[elementId];

    constexpr auto ldTinv = yateto::leadDim<seissol::init::Tinv>();
    constexpr auto ldT = yateto::leadDim<seissol::init::T>();
    constexpr auto ldDisplacement = yateto::leadDim<seissol::init::displacementRotationMatrix>();

    const int i = threadIdx.x;
    const int j = threadIdx.y;

    displacementToFaceNormal[i + j * ldDisplacement] = Tinv[(i + 6) + (j + 6) * ldTinv];
    displacementToGlobalData[i + j * ldDisplacement] = T[(i + 6) + (j + 6) * ldT];
  }
}

void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                             real** displacementToGlobalDataPtrs,
                             real** TPtrs,
                             real** TinvPtrs,
                             size_t numElements,
                             void* deviceStream) {
  dim3 block(3, 3, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelextractRotationMatrices<<<grid, block, 0, stream>>>(
    displacementToFaceNormalPtrs,
    displacementToGlobalDataPtrs,
    TPtrs,
    TinvPtrs,
    numElements);
}

__global__  void kernelInitializeTaylorSeriesForGravitationalBoundary(
  real** prevCoefficientsPtrs,
  real** integratedDisplacementNodalPtrs,
  real** rotatedFaceDisplacementPtrs,
  double deltaTInt,
  size_t numElements) {

  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    auto* prevCoefficients = prevCoefficientsPtrs[elementId];
    auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
    const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

    assert(nodal::tensor::nodes2D::Shape[0] <= yateto::leadDim<seissol::init::rotatedFaceDisplacement>());

    const int tid = threadIdx.x;
    constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];
    if (tid < num2dNodes) {
      prevCoefficients[tid] = rotatedFaceDisplacement[tid];
      integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
    }
  }
}

void initializeTaylorSeriesForGravitationalBoundary(
  real** prevCoefficientsPtrs,
  real** integratedDisplacementNodalPtrs,
  real** rotatedFaceDisplacementPtrs,
  double deltaTInt,
  size_t numElements,
  void* deviceStream) {

  dim3 block(yateto::leadDim<seissol::nodal::init::nodes2D>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelInitializeTaylorSeriesForGravitationalBoundary<<<grid, block, 0, stream>>>(
    prevCoefficientsPtrs,
    integratedDisplacementNodalPtrs,
    rotatedFaceDisplacementPtrs,
    deltaTInt,
    numElements);
}

__global__ void kernelComputeInvAcousticImpedance(double* invImpedances,
                                                  double* rhos,
                                                  double* lambdas,
                                                  size_t numElements) {

  size_t index = threadIdx.x + blockIdx.x * blockDim.x;
  if (index < numElements) {
    invImpedances[index] = 1.0 / std::sqrt(lambdas[index] * rhos[index]);
  }
}

void computeInvAcousticImpedance(double* invImpedances,
                                 double* rhos,
                                 double* lambdas,
                                 size_t numElements,
                                 void* deviceStream) {
  constexpr size_t blockSize{256};
  dim3 block(blockSize, 1, 1);
  dim3 grid((numElements + blockSize - 1) / blockSize, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelComputeInvAcousticImpedance<<<grid, block, 0, stream>>>(invImpedances,
                                                                rhos,
                                                                lambdas,
                                                                numElements);
}

__global__ void kernelUpdateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                                    real** prevCoefficientsPtrs,
                                                    real** integratedDisplacementNodalPtrs,
                                                    real** dofsFaceNodalPtrs,
                                                    double* invImpedances,
                                                    double* rhos,
                                                    double g,
                                                    double factorEvaluated,
                                                    double factorInt,
                                                    size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    constexpr int pIdx = 0;
    constexpr int uIdx = 6;
    constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];

    const int tid = threadIdx.x;
    if (tid < num2dNodes) {

      real* dofsFaceNodal = dofsFaceNodalPtrs[elementId];
      constexpr auto ldINodal = yateto::leadDim<seissol::init::INodal>();

      const auto uInside = dofsFaceNodal[tid + (uIdx + 0) * ldINodal];
      const auto vInside = dofsFaceNodal[tid + (uIdx + 1) * ldINodal];
      const auto wInside = dofsFaceNodal[tid + (uIdx + 2) * ldINodal];
      const auto pressureInside = dofsFaceNodal[tid + pIdx * ldINodal];

      real* prevCoefficients = prevCoefficientsPtrs[elementId];
#ifdef USE_ELASTIC
      const auto rho = rhos[elementId];
      const auto invImpedance = invImpedances[elementId];

      const double curCoeff = uInside - invImpedance * (rho * g * prevCoefficients[tid] + pressureInside);
#else
      const double curCoeff = uInside;
#endif
      prevCoefficients[tid] = curCoeff;

      constexpr auto ldFaceDisplacement = yateto::leadDim<seissol::init::faceDisplacement>();
      static_assert(num2dNodes <= ldFaceDisplacement, "");

      real* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];
      rotatedFaceDisplacement[tid + 0 * ldFaceDisplacement] += factorEvaluated * curCoeff;
      rotatedFaceDisplacement[tid + 1 * ldFaceDisplacement] += factorEvaluated * vInside;
      rotatedFaceDisplacement[tid + 2 * ldFaceDisplacement] += factorEvaluated * wInside;

      constexpr auto ldIntegratedFaceDisplacement = yateto::leadDim<seissol::init::averageNormalDisplacement>();
      static_assert(num2dNodes <= ldIntegratedFaceDisplacement, "");

      real* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
      integratedDisplacementNodal[tid] += factorInt * curCoeff;
    }
  }
}

void updateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                   real** prevCoefficientsPtrs,
                                   real** integratedDisplacementNodalPtrs,
                                   real** dofsFaceNodalPtrs,
                                   double* invImpedances,
                                   double* rhos,
                                   double g,
                                   double factorEvaluated,
                                   double factorInt,
                                   size_t numElements,
                                   void* deviceStream) {
  dim3 block(yateto::leadDim<seissol::nodal::init::nodes2D>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelUpdateRotatedFaceDisplacement<<<grid, block, 0, stream>>>(
    rotatedFaceDisplacementPtrs,
    prevCoefficientsPtrs,
    integratedDisplacementNodalPtrs,
    dofsFaceNodalPtrs,
    invImpedances,
    rhos,
    g,
    factorEvaluated,
    factorInt,
    numElements);
}
} // namespace seissol::kernels::time::aux