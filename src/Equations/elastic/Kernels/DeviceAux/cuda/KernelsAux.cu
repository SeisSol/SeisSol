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
template<> struct FunctionCount<5> { constexpr static size_t Count = 36; };
template<> struct FunctionCount<4> { constexpr static size_t Count = 20; };
template<> struct FunctionCount<3> { constexpr static size_t Count = 10; };
template<> struct FunctionCount<2> { constexpr static size_t Count = 4; };
template<> struct FunctionCount<1> { constexpr static size_t Count = 1; };

template<size_t Order>
constexpr size_t Functions = FunctionCount<Order>::Count;

__constant__ real derivative6X[Functions<6> * Functions<5>];
__constant__ real derivative6Y[Functions<6> * Functions<5>];
__constant__ real derivative6Z[Functions<6> * Functions<5>];

constexpr size_t AderMultiple = 3;

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
__device__ __forceinline__ void dgkernelb(real* __restrict__ output, const real* __restrict__ dQ, const real* __restrict__ coordinates, const real* __restrict__ derivative,
    real lambda, real mu, real rhoD, int offset) {
    real l2mu = lambda + 2 * mu;
    real n1 = coordinates[VECTOR(offset, 9)];
    real n2 = coordinates[VECTOR(offset+1, 9)];
    real n3 = coordinates[VECTOR(offset+2, 9)];

{
    constexpr size_t TileN = TileN1;
        for (int i = 0; i < Functions<Order>; i += TileN) {
        {
            constexpr size_t TileM = 3;
            constexpr int j = 0;
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                for (int k = 0; k < Functions<Order+1>; k += TileK) {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK - (TileK % 4); kk += 4) {
                                float4 dQLocal = *(float4*)&dQ[MATRIX4(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.x
                                    * derivative[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.y
                                    * derivative[k * TileN + (kk+1) * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.z
                                    * derivative[k * TileN + (kk+2) * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.w
                                    * derivative[k * TileN + (kk+3) * TileN + ii + i * Functions<6>];
                            }
                            }
                            #pragma unroll
                            for (int kk = TileK - (TileK % 4); kk < TileK; ++kk) {
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
                assign<Override>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0]));
                assign<Override>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n2 * result[ii][1]));
                assign<Override>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n3 * result[ii][2]));
            }
        }
    }
}
{
    constexpr size_t TileN = TileN1;
        for (int i = 0; i < Functions<Order>; i += TileN) {
        {
            constexpr size_t TileM = 3;
            constexpr int j = 3;
            real result[TileN][TileM] = { 0 };
            {

            // derivative * dQ
    // #pragma unroll
                for (int k = 0; k < Functions<Order+1>; k += TileK) {
    #pragma unroll
                        for (int jj = 0; jj < TileM; ++jj) {
    #pragma unroll
                            for (int kk = 0; kk < TileK - (TileK % 4); kk += 4) {
                                float4 dQLocal = *(float4*)&dQ[MATRIX4(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.x
                                    * derivative[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.y
                                    * derivative[k * TileN + (kk+1) * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.z
                                    * derivative[k * TileN + (kk+2) * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.w
                                    * derivative[k * TileN + (kk+3) * TileN + ii + i * Functions<6>];
                            }
                            }
                            #pragma unroll
                            for (int kk = TileK - (TileK % 4); kk < TileK; ++kk) {
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
                assign<Override>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n2 * result[ii][0] + n3 * result[ii][2]));
                assign<Override>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0] + n3 * result[ii][1]));
                assign<Override>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][2] + n2 * result[ii][1]));
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
                            for (int kk = 0; kk < TileK - (TileK % 4); kk += 4) {
                                float4 dQLocal = *(float4*)&dQ[MATRIX4(k+kk, j+jj, Functions<Order+1>, Quantities)];
    #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.x
                                    * derivative[k * TileN + kk * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.y
                                    * derivative[k * TileN + (kk+1) * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.z
                                    * derivative[k * TileN + (kk+2) * TileN + ii + i * Functions<6>];
                            }
                                #pragma unroll
                                for (int ii = 0; ii < TileN; ++ii) {
                                result[ii][jj]
                                    +=
                                    dQLocal.w
                                    * derivative[k * TileN + (kk+3) * TileN + ii + i * Functions<6>];
                            }
                            }
                            #pragma unroll
                            for (int kk = TileK - (TileK % 4); kk < TileK; ++kk) {
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
__device__ __forceinline__ void dgkernel2(real* __restrict__ output, const real* __restrict__ dQ, const real* __restrict__ coordinates, const real* __restrict__ derivativeX,
const real* __restrict__ derivativeY, const real* __restrict__ derivativeZ,
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
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0] + n2 * result[ii][3] + n3 * result[ii][5]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][3] + n2 * result[ii][1] + n3 * result[ii][4]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][5] + n2 * result[ii][4] + n3 * result[ii][2]));
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
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0] + n2 * result[ii][3] + n3 * result[ii][5]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][3] + n2 * result[ii][1] + n3 * result[ii][4]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][5] + n2 * result[ii][4] + n3 * result[ii][2]));
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
                assign<false>(output[MATRIX(i+ii, 6, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][0] + n2 * result[ii][3] + n3 * result[ii][5]));
                assign<false>(output[MATRIX(i+ii, 7, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][3] + n2 * result[ii][1] + n3 * result[ii][4]));
                assign<false>(output[MATRIX(i+ii, 8, Functions<Order>, Quantities)], -rhoD * (n1 * result[ii][5] + n2 * result[ii][4] + n3 * result[ii][2]));
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
                                    * derivativeZ[k * TileN + kk * TileN + ii + i * Functions<6>];
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
      // dgkernelb<36, 36, 7, 5, true>(dQ5, Q, coordinates, derivative6X, lambda, mu, rhoD, 0);
      // dgkernelb<36, 36, 7, 5, false>(dQ5, Q, coordinates, derivative6Y, lambda, mu, rhoD, 3);
      // dgkernelb<36, 36, 7, 5, false>(dQ5, Q, coordinates, derivative6Z, lambda, mu, rhoD, 6);
      sumkernel<5, 12, true>(I, dQ5, scale1);
      // sumkernel<5, true>(I, dQ5, scale1);
      scale1 *= scale;
      scale2 *= scale;

      // 36x9 -> 20x9
      dgkernel2b<10, 10, 12, 4, true>(dQ4, dQ5, coordinates, derivative6X, derivative6Y, derivative6Z, lambda, mu, rhoD, 0);
      // dgkernelb<20, 20, 6, 4, true>(dQ4, dQ5, coordinates, derivative6X, lambda, mu, rhoD, 0);
      // dgkernelb<20, 20, 6, 4, false>(dQ4, dQ5, coordinates, derivative6Y, lambda, mu, rhoD, 3);
      // dgkernelb<20, 20, 6, 4, false>(dQ4, dQ5, coordinates, derivative6Z, lambda, mu, rhoD, 6);
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
void aderLauncher(std::size_t count, real timestep, const real* dofs, real* buffers, real* derivatives, const real* stardata, const real* coordinates, void* stream) {
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
  kernelXTiled3x32 <<<grid, block, 0, streamObject>>> (blocks, timestep, buffers, dofs, derivatives, derivatives + dataOffsets[1], derivatives + dataOffsets[2], derivatives + dataOffsets[3], derivatives + dataOffsets[4], derivatives + dataOffsets[5], stardata, coordinates);
  CHECK_ERR;
}

template<typename T>
__device__ __forceinline__ T& iacc1(T* __restrict__ data, size_t size, size_t inblock, size_t block) {
    return data[block * size + inblock + threadIdx.x];
}

constexpr size_t InterleaveMultiple = 3;

template<typename T>
__global__ __launch_bounds__(Blocksize * InterleaveMultiple) void interleave(const T** source, T* target, size_t size, size_t count) {
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
              if (i + j < size) {
                target[size * Blocksize * block + (i + j) * Blocksize + threadIdx.x] = swap[threadIdx.x * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + i];
              }
            }
        }
      }
    }
}

template<typename T>
__global__ __launch_bounds__(Blocksize * InterleaveMultiple) void deinterleave(const T* source, T** target, size_t size, size_t count) {
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
              if (i + j < size) {
                swap[i * (Blocksize * InterleaveMultiple + 1) + threadIdx.y * Blocksize + threadIdx.x] = source[size * Blocksize * block + (i + j) * Blocksize + threadIdx.x];
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

void interleaveLauncher(std::size_t count, std::size_t size, const real** indata, real* outdata, void* stream) {
  cudaStream_t streamObject = reinterpret_cast<cudaStream_t>(stream);
  dim3 grid((count + Blocksize * InterleaveMultiple - 1) / (Blocksize * InterleaveMultiple), 1, 1);
  dim3 block(Blocksize, InterleaveMultiple, 1);
  interleave<<<grid, block, 0, streamObject>>>(indata, outdata, size, count);
  CHECK_ERR;
}
void deinterleaveLauncher(std::size_t count, std::size_t size, const real* indata, real** outdata, void* stream) {
  cudaStream_t streamObject = reinterpret_cast<cudaStream_t>(stream);
  dim3 grid((count + Blocksize * InterleaveMultiple - 1) / (Blocksize * InterleaveMultiple), 1, 1);
  dim3 block(Blocksize, InterleaveMultiple, 1);
  deinterleave<<<grid, block, 0, streamObject>>>(indata, outdata, size, count);
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