// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include <cstdio>
#include <init.h>
#include <tensor.h>
#include <yateto.h>

#include <Solver/MultipleSimulations.h>

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

#ifdef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
namespace {
constexpr std::size_t Blocksize = 128;

template <typename SourceRealT>
constexpr std::size_t RestFunctions(std::size_t SourceOrder, std::size_t ThisOrder) {
  std::size_t total = 0;
  for (std::size_t j = ThisOrder; j < SourceOrder; ++j) {
    total +=
        seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder - ThisOrder);
  }
  return total;
}

template <std::size_t Quantities,
          std::size_t ThisOrder,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          std::size_t Offset,
          std::size_t SharedOffset,
          typename... Coeffs>
static __device__ __forceinline__ void taylorSumInner(TargetRealT* const __restrict__ target,
                                                      const SourceRealT* const __restrict__ source,
                                                      SourceRealT* const __restrict__ shmem,
                                                      TargetRealT reg[Quantities],
                                                      TargetRealT coeff,
                                                      Coeffs... coeffs) {
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

  static_assert(seissol::tensor::dQ::size(ThisOrder) == SourceMemSize,
                "Tensor size mismatch in explicit kernel.");

  if constexpr (LoadSize > 0) {
    __syncthreads();
    constexpr std::size_t Rounds = LoadSize / Blocksize;
    constexpr std::size_t Rest = LoadSize % Blocksize;
    if constexpr (Rounds > 0) {
#pragma unroll
      for (std::size_t j = 0; j < Rounds; ++j) {
        shmem[j * Blocksize + threadIdx.x] = source[Offset + j * Blocksize + threadIdx.x];
      }
    }
    if constexpr (Rest > 0) {
      if (threadIdx.x < Rest) {
        shmem[Rounds * Blocksize + threadIdx.x] = source[Offset + Rounds * Blocksize + threadIdx.x];
      }
    }
    __syncthreads();
  }

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
  if constexpr (ThisOrder + 1 < std::min(SourceOrder, TargetOrder) && sizeof...(Coeffs) > 0) {
    constexpr std::size_t SharedPosition = UseShared ? SharedOffset + SourceMemSize : 0;
    taylorSumInner<Quantities,
                   ThisOrder + 1,
                   SourceRealT,
                   TargetRealT,
                   SourceOrder,
                   TargetOrder,
                   Offset + LoadSize,
                   SharedPosition>(target, source, shmem, reg, coeffs...);
  }
}

template <std::size_t Quantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          typename... Coeffs>
static __global__ __launch_bounds__(Blocksize) void taylorSumKernel(TargetRealT** targetBatch,
                                                                    const SourceRealT** sourceBatch,
                                                                    Coeffs... coeffs) {
  int batchId = blockIdx.x;

  __shared__ SourceRealT
      shmem[seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder) *
            Quantities];
  TargetRealT reg[Quantities] = {0};

  const SourceRealT* const __restrict__ source =
      const_cast<const SourceRealT*>(sourceBatch[batchId]);
  TargetRealT* const __restrict__ target = targetBatch[batchId];

  taylorSumInner<Quantities, 0, SourceRealT, TargetRealT, SourceOrder, TargetOrder, 0, 0>(
      target, source, shmem, reg, coeffs...);

  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

  if (threadIdx.x < TargetStride) {
#pragma unroll
    for (std::size_t j = 0; j < Quantities; ++j) {
      target[TargetStride * j + threadIdx.x] = reg[j];
    }
  }
}

template <std::size_t SourceQuantities,
          std::size_t TargetQuantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          std::size_t... Indices>
void static taylorSumInternal(std::size_t count,
                              TargetRealT** target,
                              const SourceRealT** source,
                              const TargetRealT* coeffs,
                              void* stream,
                              std::index_sequence<Indices...> indices) {
  constexpr std::size_t Quantities = std::min(SourceQuantities, TargetQuantities);
  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

  dim3 threads(Blocksize);
  dim3 blocks(count);

  StreamT castedStream = reinterpret_cast<StreamT>(stream);

  taylorSumKernel<Quantities, SourceRealT, TargetRealT, SourceOrder, TargetOrder>
      <<<blocks, threads, 0, castedStream>>>(target, source, coeffs[Indices]...);
}
} // namespace

namespace seissol::kernels::time::aux {
void taylorSum(
    std::size_t count, real** target, const real** source, const real* coeffs, void* stream) {
  taylorSumInternal<seissol::model::MaterialT::NumQuantities,
                    seissol::model::MaterialT::NumQuantities,
                    real,
                    real,
                    ConvergenceOrder,
                    ConvergenceOrder>(
      count, target, source, coeffs, stream, std::make_index_sequence<ConvergenceOrder>());
}
} // namespace seissol::kernels::time::aux
#endif

namespace {
using namespace seissol::multisim;

template <typename Tensor>
constexpr size_t leadDim() {
  if constexpr (MultisimEnabled) {
    return Tensor::Stop[1] - Tensor::Start[1];
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

template <typename Tensor>
constexpr size_t linearDim() {
  if constexpr (MultisimEnabled) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

constexpr auto getblock(int size) {
  if constexpr (MultisimEnabled) {
    return dim3(NumSimulations, size);
  } else {
    return dim3(size);
  }
}

__forceinline__ __device__ auto linearidx() {
  if constexpr (MultisimEnabled) {
    return threadIdx.y * NumSimulations + threadIdx.x;
  } else {
    return threadIdx.x;
  }
}

__forceinline__ __device__ auto linearsize() {
  if constexpr (MultisimEnabled) {
    return blockDim.y * blockDim.x;
  } else {
    return blockDim.x;
  }
}

__forceinline__ __device__ auto simidx() {
  if constexpr (MultisimEnabled) {
    return threadIdx.x;
  } else {
    return 0;
  }
}

__forceinline__ __device__ auto validx() {
  if constexpr (MultisimEnabled) {
    return threadIdx.y;
  } else {
    return threadIdx.x;
  }
}
} // namespace

namespace seissol::kernels::local_flux::aux::details {

__global__ void kernelFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                                         real** displacementDataPtrs,
                                         double* rhos,
                                         double g,
                                         size_t numElements) {

  const int tid = linearidx();
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    const double rho = rhos[elementId];
    real* elementBoundaryDofs = dofsFaceBoundaryNodalPtrs[elementId];
    real* elementDisplacement = displacementDataPtrs[elementId];

    constexpr auto numNodes = linearDim<seissol::nodal::init::nodes2D>();
    if (tid < numNodes) {
      constexpr auto ldINodal = linearDim<seissol::init::INodal>();

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
  dim3 block = getblock(leadDim<seissol::nodal::init::nodes2D>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelFreeSurfaceGravity<<<grid, block, 0, stream>>>(
      dofsFaceBoundaryNodalPtrs, displacementDataPtrs, rhos, g, numElements);
}

__global__ void kernelEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                                   real** easiBoundaryMapPtrs,
                                   real** easiBoundaryConstantPtrs,
                                   size_t numElements) {

  const int tid = validx();
  const int elementId = blockIdx.x;

  constexpr auto ldINodalDim = linearDim<seissol::init::INodal>();
  constexpr auto iNodalDim0 = seissol::tensor::INodal::Shape[multisim::BasisFunctionDimension + 0];
  constexpr auto iNodalDim1 = seissol::tensor::INodal::Shape[multisim::BasisFunctionDimension + 1];
  __shared__ __align__(8) real resultTerm[iNodalDim1][iNodalDim0][multisim::NumSimulations];

  constexpr auto ldConstantDim = linearDim<seissol::init::easiBoundaryConstant>();
  constexpr auto constantDim0 =
      seissol::tensor::easiBoundaryConstant::Shape[multisim::BasisFunctionDimension + 0];
  constexpr auto constantDim1 =
      seissol::tensor::easiBoundaryConstant::Shape[multisim::BasisFunctionDimension + 1];
  __shared__ __align__(8) real rightTerm[constantDim0][constantDim1][multisim::NumSimulations];

  constexpr auto ldMapDim = leadDim<seissol::init::easiBoundaryMap>();
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

    for (int i = linearidx(); i < (ldConstantDim * constantDim1); i += linearsize()) {
      const auto sim = i % multisim::NumSimulations;
      const auto subsim = i / multisim::NumSimulations;
      const auto quantity = subsim % constantDim0;
      const auto quadpoint = subsim / constantDim0;
      rightTerm[quantity][quadpoint][sim] = easiBoundaryConstant[i];
    }
    __syncthreads();

    for (int i = 0; i < iNodalDim1; ++i) {
      if (tid < iNodalDim0) {
        resultTerm[i][tid][simidx()] = 0.0;
      }
    }
    __syncthreads();

    for (int quantity = 0; quantity < mapDim1; ++quantity) {
      for (int quadpoint = 0; quadpoint < mapDim2; ++quadpoint) {
        if (tid < mapDim0) {
          leftTerm[tid][quadpoint] =
              easiBoundaryMap[tid + ldMapDim * (quantity + quadpoint * mapDim1)];
        }
      }
      __syncthreads();

      if (tid < mapDim2) {
        const real col = dofsFaceBoundaryNodal[linearidx() + quantity * ldINodalDim];
        for (int quantity2 = 0; quantity2 < mapDim0; ++quantity2) {
          resultTerm[quantity2][tid][simidx()] += leftTerm[quantity2][tid] * col;
        }
      }
      __syncthreads();
    }

    if (tid < iNodalDim0) {
      for (int quantity2 = 0; quantity2 < iNodalDim1; ++quantity2) {
        dofsFaceBoundaryNodal[linearidx() + quantity2 * ldINodalDim] =
            resultTerm[quantity2][tid][simidx()] + rightTerm[quantity2][tid][simidx()];
      }
    }
  }
}

void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream) {

  dim3 block = getblock(leadDim<seissol::init::INodal>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelEasiBoundary<<<grid, block, 0, stream>>>(
      dofsFaceBoundaryNodalPtrs, easiBoundaryMapPtrs, easiBoundaryConstantPtrs, numElements);
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
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelextractRotationMatrices<<<grid, block, 0, stream>>>(
      displacementToFaceNormalPtrs, displacementToGlobalDataPtrs, TPtrs, TinvPtrs, numElements);
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

    assert(linearDim<seissol::nodal::init::nodes2D>() <=
           linearDim<seissol::init::rotatedFaceDisplacement>());

    const int tid = linearidx();
    constexpr auto num2dNodes = linearDim<seissol::nodal::init::nodes2D>();
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

  dim3 block = getblock(leadDim<seissol::nodal::init::nodes2D>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
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

void computeInvAcousticImpedance(
    double* invImpedances, double* rhos, double* lambdas, size_t numElements, void* deviceStream) {
  constexpr size_t blockSize{1024};
  dim3 block(blockSize, 1, 1);
  dim3 grid((numElements + blockSize - 1) / blockSize, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelComputeInvAcousticImpedance<<<grid, block, 0, stream>>>(
      invImpedances, rhos, lambdas, numElements);
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
    constexpr int uIdx = model::MaterialT::TractionQuantities;
    constexpr auto num2dNodes = linearDim<seissol::nodal::init::nodes2D>();

    const int tid = linearidx();
    if (tid < num2dNodes) {

      real* dofsFaceNodal = dofsFaceNodalPtrs[elementId];
      constexpr auto ldINodal = linearDim<seissol::init::INodal>();

      const auto uInside = dofsFaceNodal[tid + (uIdx + 0) * ldINodal];
      const auto vInside = dofsFaceNodal[tid + (uIdx + 1) * ldINodal];
      const auto wInside = dofsFaceNodal[tid + (uIdx + 2) * ldINodal];
      const auto pressureInside = dofsFaceNodal[tid + pIdx * ldINodal];

      real* prevCoefficients = prevCoefficientsPtrs[elementId];
      const auto rho = rhos[elementId];
      const auto invImpedance = invImpedances[elementId];

      const double curCoeff =
          uInside - invImpedance * (rho * g * prevCoefficients[tid] + pressureInside);
      prevCoefficients[tid] = curCoeff;

      constexpr auto ldFaceDisplacement = linearDim<seissol::init::faceDisplacement>();
      static_assert(num2dNodes <= ldFaceDisplacement, "");

      real* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];
      rotatedFaceDisplacement[tid + 0 * ldFaceDisplacement] += factorEvaluated * curCoeff;
      rotatedFaceDisplacement[tid + 1 * ldFaceDisplacement] += factorEvaluated * vInside;
      rotatedFaceDisplacement[tid + 2 * ldFaceDisplacement] += factorEvaluated * wInside;

      constexpr auto ldIntegratedFaceDisplacement =
          linearDim<seissol::init::averageNormalDisplacement>();
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
  dim3 block = getblock(leadDim<seissol::nodal::init::nodes2D>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelUpdateRotatedFaceDisplacement<<<grid, block, 0, stream>>>(rotatedFaceDisplacementPtrs,
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
