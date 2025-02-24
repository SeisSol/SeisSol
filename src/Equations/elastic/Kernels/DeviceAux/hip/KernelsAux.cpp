// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include <cstdio>
#include <hip/hip_runtime.h>
#include <init.h>
#include <tensor.h>
#include <yateto.h>

#ifdef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
namespace {
constexpr std::size_t Blocksize = 64;

template <typename SourceRealT>
constexpr std::size_t RestFunctions(std::size_t SourceOrder, std::size_t ThisOrder) {
  std::size_t total = 0;
  for (std::size_t j = ThisOrder; j < SourceOrder; ++j) {
    total +=
        seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder - ThisOrder);
  }
  return total;
}

template <bool Integral,
          std::size_t Quantities,
          std::size_t ThisOrder,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          std::size_t Offset,
          std::size_t SharedOffset>
static __device__ __forceinline__ void taylorSumInner(TargetRealT* const __restrict__ target,
                                                      const SourceRealT* const __restrict__ source,
                                                      TargetRealT start,
                                                      TargetRealT end,
                                                      TargetRealT startCoeff,
                                                      TargetRealT endCoeff,
                                                      SourceRealT* const __restrict__ shmem,
                                                      TargetRealT reg[Quantities]) {
  constexpr std::size_t MemorySize =
      seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder) * Quantities;
  constexpr std::size_t SourceStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder - ThisOrder);
  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);
  constexpr std::size_t RestMemSize =
      SharedOffset > 0 ? 0 : RestFunctions<SourceRealT>(SourceOrder, ThisOrder) * Quantities;
  constexpr std::size_t SourceMemSize = SourceStride * Quantities;
  constexpr bool UseShared = MemorySize >= RestMemSize;
  constexpr std::size_t LoadSize = UseShared ? RestMemSize : SourceMemSize;
  constexpr TargetRealT DivisionCoefficient =
      Integral ? static_cast<TargetRealT>(ThisOrder + 2) : static_cast<TargetRealT>(ThisOrder + 1);

  static_assert(seissol::tensor::dQ::size(ThisOrder) == SourceMemSize,
                "Tensor size mismatch in explicit kernel.");

  if constexpr (LoadSize > 0) {
    __syncthreads();
    constexpr std::size_t Rounds = LoadSize / Blocksize;
    constexpr std::size_t Rest = LoadSize % Blocksize;
#pragma unroll
    for (std::size_t j = 0; j < Rounds; ++j) {
      shmem[j * Blocksize + threadIdx.x] = source[Offset + j * Blocksize + threadIdx.x];
    }
    if constexpr (Rest > 0) {
      if (threadIdx.x < Rest) {
        shmem[Rounds * Blocksize + threadIdx.x] = source[Offset + Rounds * Blocksize + threadIdx.x];
      }
    }
    __syncthreads();
  }

  const TargetRealT coeff = endCoeff - startCoeff;
  constexpr std::size_t BasisFunctionsSize = std::min(SourceStride, TargetStride);
  if (threadIdx.x < BasisFunctionsSize) {
#pragma unroll
    for (std::size_t j = 0; j < Quantities; ++j) {
      // FIXME: non-optimal warp utilization (as before... But now it's in registers)
      reg[j] +=
          coeff * static_cast<TargetRealT>(shmem[SharedOffset + SourceStride * j + threadIdx.x]);
    }
  }
  // TODO(David): are we sure about this? Or is ThisOrder + 1 < SourceOrder enough?
  if constexpr (ThisOrder + 1 < std::min(SourceOrder, TargetOrder)) {
    constexpr std::size_t SharedPosition = UseShared ? SharedOffset + SourceMemSize : 0;
    const TargetRealT newStartCoeff = startCoeff * start / DivisionCoefficient;
    const TargetRealT newEndCoeff = endCoeff * end / DivisionCoefficient;
    taylorSumInner<Integral,
                   Quantities,
                   ThisOrder + 1,
                   SourceRealT,
                   TargetRealT,
                   SourceOrder,
                   TargetOrder,
                   Offset + LoadSize,
                   SharedPosition>(
        target, source, start, end, newStartCoeff, newEndCoeff, shmem, reg);
  }
}

template <bool Integral,
          std::size_t Quantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder>
static __global__ __launch_bounds__(Blocksize) void taylorSumKernel(TargetRealT** targetBatch,
                                                                    const SourceRealT** sourceBatch,
                                                                    TargetRealT start,
                                                                    TargetRealT end) {
  int batchId = blockIdx.x;

  __shared__ SourceRealT
      shmem[seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder) *
            Quantities];
  TargetRealT reg[Quantities] = {0};

  const SourceRealT* const __restrict__ source =
      const_cast<const SourceRealT*>(sourceBatch[batchId]);
  TargetRealT* const __restrict__ target = targetBatch[batchId];

  const TargetRealT startCoeff = Integral ? start : 0;
  const TargetRealT endCoeff = Integral ? end : 1;

  taylorSumInner<Integral, Quantities, 0, SourceRealT, TargetRealT, SourceOrder, TargetOrder, 0, 0>(
      target, source, start, end, startCoeff, endCoeff, shmem, reg);

  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

  if (threadIdx.x < TargetStride) {
#pragma unroll
    for (std::size_t j = 0; j < Quantities; ++j) {
      target[TargetStride * j + threadIdx.x] = reg[j];
    }
  }
}

template <bool Integral,
          std::size_t SourceQuantities,
          std::size_t TargetQuantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder>
void static taylorSumInternal(std::size_t count,
                              TargetRealT** target,
                              const SourceRealT** source,
                              TargetRealT start,
                              TargetRealT end,
                              void* stream) {
  constexpr std::size_t Quantities = std::min(SourceQuantities, TargetQuantities);
  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

  dim3 threads(Blocksize);
  dim3 blocks(count);

  hipStream_t castedStream = reinterpret_cast<hipStream_t>(stream);

  taylorSumKernel<Integral, Quantities, SourceRealT, TargetRealT, SourceOrder, TargetOrder>
      <<<blocks, threads, 0, castedStream>>>(target, source, start, end);
}
} // namespace

namespace seissol::kernels::time::aux {
void taylorSum(bool integral,
               std::size_t count,
               real** target,
               const real** source,
               real start,
               real end,
               void* stream) {
  if (integral) {
    taylorSumInternal<true,
                      seissol::model::MaterialT::NumQuantities,
                      seissol::model::MaterialT::NumQuantities,
                      real,
                      real,
                      ConvergenceOrder,
                      ConvergenceOrder>(count, target, source, start, end, stream);
  } else {
    taylorSumInternal<false,
                      seissol::model::MaterialT::NumQuantities,
                      seissol::model::MaterialT::NumQuantities,
                      real,
                      real,
                      ConvergenceOrder,
                      ConvergenceOrder>(count, target, source, start, end, stream);
  }
}
} // namespace seissol::kernels::time::aux
#endif

namespace seissol::kernels::local_flux::aux::details {

__global__ void kernelFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
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
  auto stream = reinterpret_cast<hipStream_t>(deviceStream);
  hipLaunchKernelGGL(kernelFreeSurfaceGravity,
                     dim3(grid),
                     dim3(block),
                     0,
                     stream,
                     dofsFaceBoundaryNodalPtrs,
                     displacementDataPtrs,
                     rhos,
                     g,
                     numElements);
}

__global__ void kernelEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
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
      if (tid < iNodalDim0)
        resultTerm[i][tid] = 0.0;
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
  auto stream = reinterpret_cast<hipStream_t>(deviceStream);
  hipLaunchKernelGGL(kernelEasiBoundary,
                     dim3(grid),
                     dim3(block),
                     0,
                     stream,
                     dofsFaceBoundaryNodalPtrs,
                     easiBoundaryMapPtrs,
                     easiBoundaryConstantPtrs,
                     numElements);
}
} // namespace seissol::kernels::local_flux::aux::details

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
  auto stream = reinterpret_cast<hipStream_t>(deviceStream);
  hipLaunchKernelGGL(kernelextractRotationMatrices,
                     dim3(grid),
                     dim3(block),
                     0,
                     stream,
                     displacementToFaceNormalPtrs,
                     displacementToGlobalDataPtrs,
                     TPtrs,
                     TinvPtrs,
                     numElements);
}

__global__ void
    kernelInitializeTaylorSeriesForGravitationalBoundary(real** prevCoefficientsPtrs,
                                                         real** integratedDisplacementNodalPtrs,
                                                         real** rotatedFaceDisplacementPtrs,
                                                         double deltaTInt,
                                                         size_t numElements) {

  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    auto* prevCoefficients = prevCoefficientsPtrs[elementId];
    auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
    const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

    assert(nodal::tensor::nodes2D::Shape[0] <=
           yateto::leadDim<seissol::init::rotatedFaceDisplacement>());

    const int tid = threadIdx.x;
    constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];
    if (tid < num2dNodes) {
      prevCoefficients[tid] = rotatedFaceDisplacement[tid];
      integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
    }
  }
}

void initializeTaylorSeriesForGravitationalBoundary(real** prevCoefficientsPtrs,
                                                    real** integratedDisplacementNodalPtrs,
                                                    real** rotatedFaceDisplacementPtrs,
                                                    double deltaTInt,
                                                    size_t numElements,
                                                    void* deviceStream) {

  dim3 block(yateto::leadDim<seissol::nodal::init::nodes2D>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<hipStream_t>(deviceStream);
  hipLaunchKernelGGL(kernelInitializeTaylorSeriesForGravitationalBoundary,
                     dim3(grid),
                     dim3(block),
                     0,
                     stream,
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

void computeInvAcousticImpedance(
    double* invImpedances, double* rhos, double* lambdas, size_t numElements, void* deviceStream) {
  constexpr size_t blockSize{256};
  dim3 block(blockSize, 1, 1);
  dim3 grid((numElements + blockSize - 1) / blockSize, 1, 1);
  auto stream = reinterpret_cast<hipStream_t>(deviceStream);
  hipLaunchKernelGGL(kernelComputeInvAcousticImpedance,
                     dim3(grid),
                     dim3(block),
                     0,
                     stream,
                     invImpedances,
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

      const double curCoeff =
          uInside - invImpedance * (rho * g * prevCoefficients[tid] + pressureInside);
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

      constexpr auto ldIntegratedFaceDisplacement =
          yateto::leadDim<seissol::init::averageNormalDisplacement>();
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
  auto stream = reinterpret_cast<hipStream_t>(deviceStream);
  hipLaunchKernelGGL(kernelUpdateRotatedFaceDisplacement,
                     dim3(grid),
                     dim3(block),
                     0,
                     stream,
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
