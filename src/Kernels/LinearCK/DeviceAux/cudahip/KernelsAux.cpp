// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/LinearCK/DeviceAux/KernelsAux.h"

#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Solver/MultipleSimulations.h"

#include <cstdio>
#include <yateto.h>

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

  static_assert(seissol::tensor::dQ<Cfg>::size(ThisOrder) == SourceMemSize,
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
void taylorSum(std::size_t count,
               Real<Cfg>** target,
               const Real<Cfg>** source,
               const Real<Cfg>* coeffs,
               void* stream) {
  taylorSumInternal<seissol::model::MaterialT::NumQuantities,
                    seissol::model::MaterialT::NumQuantities,
                    Real<Cfg>,
                    Real<Cfg>,
                    Cfg::ConvergenceOrder,
                    Cfg::ConvergenceOrder>(
      count, target, source, coeffs, stream, std::make_index_sequence<Cfg::ConvergenceOrder>());
}
} // namespace seissol::kernels::time::aux
#endif

namespace {
using namespace seissol::multisim;

template <typename Cfg, typename Tensor>
constexpr size_t leadDim() {
  if constexpr (MultisimEnabled<Cfg>) {
    return Tensor::Stop[1] - Tensor::Start[1];
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

template <typename Cfg, typename Tensor>
constexpr size_t linearDim() {
  if constexpr (MultisimEnabled<Cfg>) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

template <typename Cfg>
constexpr auto getblock(int size) {
  if constexpr (MultisimEnabled<Cfg>) {
    return dim3(NumSimulations<Cfg>, size);
  } else {
    return dim3(size);
  }
}

template <typename Cfg>
__forceinline__ __device__ auto linearidx() {
  if constexpr (MultisimEnabled<Cfg>) {
    return threadIdx.y * NumSimulations<Cfg> + threadIdx.x;
  } else {
    return threadIdx.x;
  }
}

template <typename Cfg>
__forceinline__ __device__ auto linearsize() {
  if constexpr (MultisimEnabled<Cfg>) {
    return blockDim.y * blockDim.x;
  } else {
    return blockDim.x;
  }
}

template <typename Cfg>
__forceinline__ __device__ auto simidx() {
  if constexpr (MultisimEnabled<Cfg>) {
    return threadIdx.x;
  } else {
    return 0;
  }
}

template <typename Cfg>
__forceinline__ __device__ auto validx() {
  if constexpr (MultisimEnabled<Cfg>) {
    return threadIdx.y;
  } else {
    return threadIdx.x;
  }
}
} // namespace

namespace seissol::kernels::local_flux::aux {

template <typename Cfg>
__global__ void kernelFreeSurfaceGravity(Real<Cfg>** dofsFaceBoundaryNodalPtrs,
                                         Real<Cfg>** displacementDataPtrs,
                                         double* rhos,
                                         double g,
                                         size_t numElements) {

  const int tid = linearidx<Cfg>();
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    const double rho = rhos[elementId];
    Real<Cfg>* elementBoundaryDofs = dofsFaceBoundaryNodalPtrs[elementId];
    Real<Cfg>* elementDisplacement = displacementDataPtrs[elementId];

    constexpr auto numNodes = linearDim<Cfg, seissol::nodal::init::nodes2D<Cfg>>();
    if (tid < numNodes) {
      constexpr auto ldINodal = linearDim<Cfg, seissol::init::INodal<Cfg>>();

      const auto pressureAtBnd = static_cast<Real<Cfg>>(-1.0) * rho * g * elementDisplacement[tid];

#pragma unroll
      for (int component{0}; component < 3; ++component) {
        elementBoundaryDofs[tid + component * ldINodal] =
            2.0 * pressureAtBnd - elementBoundaryDofs[tid + component * ldINodal];
      }
    }
  }
}

template <typename Cfg>
void FreeSurfaceGravity<Cfg>::dispatch(Real<Cfg>** dofsFaceBoundaryNodalPtrs,
                                       size_t numElements,
                                       void* deviceStream) {
  assert(displacementDataPtrs != nullptr);
  assert(rhos != nullptr);
  dim3 block = getblock<Cfg>(leadDim<Cfg, seissol::nodal::init::nodes2D<Cfg>>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelFreeSurfaceGravity<Cfg><<<grid, block, 0, stream>>>(
      dofsFaceBoundaryNodalPtrs, displacementDataPtrs, rhos, g, numElements);
}

template <typename Cfg>
__global__ void kernelEasiBoundary(Real<Cfg>** dofsFaceBoundaryNodalPtrs,
                                   Real<Cfg>** easiBoundaryMapPtrs,
                                   Real<Cfg>** easiBoundaryConstantPtrs,
                                   size_t numElements) {

  const int tid = validx<Cfg>();
  const int elementId = blockIdx.x;

  constexpr auto ldINodalDim = linearDim<Cfg, seissol::init::INodal<Cfg>>();
  constexpr auto iNodalDim0 = seissol::tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg> + 0];
  constexpr auto iNodalDim1 = seissol::tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  __shared__ __align__(8) Real<Cfg> resultTerm[iNodalDim1][iNodalDim0]
                                              [multisim::NumSimulations<Cfg>];

  constexpr auto ldConstantDim = linearDim<Cfg, seissol::init::easiBoundaryConstant<Cfg>>();
  constexpr auto constantDim0 =
      seissol::tensor::easiBoundaryConstant<Cfg>::Shape[multisim::BasisDim<Cfg> + 0];
  constexpr auto constantDim1 =
      seissol::tensor::easiBoundaryConstant<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  __shared__ __align__(8) Real<Cfg> rightTerm[constantDim0][constantDim1]
                                             [multisim::NumSimulations<Cfg>];

  constexpr auto ldMapDim = leadDim<Cfg, seissol::init::easiBoundaryMap<Cfg>>();
  constexpr auto mapDim0 = seissol::tensor::easiBoundaryMap<Cfg>::Shape[0];
  constexpr auto mapDim1 = seissol::tensor::easiBoundaryMap<Cfg>::Shape[1];
  constexpr auto mapDim2 = seissol::tensor::easiBoundaryMap<Cfg>::Shape[2];
  __shared__ __align__(8) Real<Cfg> leftTerm[mapDim0][mapDim2];

  static_assert(iNodalDim1 == constantDim0, "supposed to be equal");
  static_assert(iNodalDim1 == mapDim0, "supposed to be equal");

  if (elementId < numElements) {
    Real<Cfg>* dofsFaceBoundaryNodal = dofsFaceBoundaryNodalPtrs[elementId];
    Real<Cfg>* easiBoundaryMap = easiBoundaryMapPtrs[elementId];
    auto easiBoundaryConstant = easiBoundaryConstantPtrs[elementId];

    for (int i = linearidx<Cfg>(); i < (ldConstantDim * constantDim1); i += linearsize<Cfg>()) {
      const auto sim = i % multisim::NumSimulations<Cfg>;
      const auto subsim = i / multisim::NumSimulations<Cfg>;
      const auto quantity = subsim % constantDim0;
      const auto quadpoint = subsim / constantDim0;
      rightTerm[quantity][quadpoint][sim] = easiBoundaryConstant[i];
    }
    __syncthreads();

    for (int i = 0; i < iNodalDim1; ++i) {
      if (tid < iNodalDim0) {
        resultTerm[i][tid][simidx<Cfg>()] = 0.0;
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
        const Real<Cfg> col = dofsFaceBoundaryNodal[linearidx<Cfg>() + quantity * ldINodalDim];
        for (int quantity2 = 0; quantity2 < mapDim0; ++quantity2) {
          resultTerm[quantity2][tid][simidx<Cfg>()] += leftTerm[quantity2][tid] * col;
        }
      }
      __syncthreads();
    }

    if (tid < iNodalDim0) {
      for (int quantity2 = 0; quantity2 < iNodalDim1; ++quantity2) {
        dofsFaceBoundaryNodal[linearidx<Cfg>() + quantity2 * ldINodalDim] =
            resultTerm[quantity2][tid][simidx<Cfg>()] + rightTerm[quantity2][tid][simidx<Cfg>()];
      }
    }
  }
}

template <typename Cfg>
void EasiBoundary<Cfg>::dispatch(Real<Cfg>** dofsFaceBoundaryNodalPtrs,
                                 size_t numElements,
                                 void* deviceStream) {
  assert(easiBoundaryMapPtrs != nullptr);
  assert(easiBoundaryConstantPtrs != nullptr);

  dim3 block = getblock<Cfg>(leadDim<Cfg, seissol::init::INodal<Cfg>>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelEasiBoundary<Cfg><<<grid, block, 0, stream>>>(
      dofsFaceBoundaryNodalPtrs, easiBoundaryMapPtrs, easiBoundaryConstantPtrs, numElements);
}

#define SEISSOL_CONFIGITER(cfg) template class FreeSurfaceGravity<cfg>;
#include "ConfigIncludeLinearCK.h"

#define SEISSOL_CONFIGITER(cfg) template class FreeSurfaceGravity<cfg>;
#include "ConfigIncludeSTP.h"

#define SEISSOL_CONFIGITER(cfg) template class EasiBoundary<cfg>;
#include "ConfigIncludeLinearCK.h"

#define SEISSOL_CONFIGITER(cfg) template class EasiBoundary<cfg>;
#include "ConfigIncludeSTP.h"

} // namespace seissol::kernels::local_flux::aux

namespace seissol::kernels::time::aux {
template <typename Cfg>
__global__ void kernelextractRotationMatrices(Real<Cfg>** displacementToFaceNormalPtrs,
                                              Real<Cfg>** displacementToGlobalDataPtrs,
                                              Real<Cfg>** TPtrs,
                                              Real<Cfg>** TinvPtrs,
                                              size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    Real<Cfg>* displacementToFaceNormal = displacementToFaceNormalPtrs[elementId];
    Real<Cfg>* displacementToGlobalData = displacementToGlobalDataPtrs[elementId];
    auto* T = TPtrs[elementId];
    auto* Tinv = TinvPtrs[elementId];

    constexpr auto ldTinv = yateto::leadDim<seissol::init::Tinv<Cfg>>();
    constexpr auto ldT = yateto::leadDim<seissol::init::T<Cfg>>();
    constexpr auto ldDisplacement =
        yateto::leadDim<seissol::init::displacementRotationMatrix<Cfg>>();

    const int i = threadIdx.x;
    const int j = threadIdx.y;

    displacementToFaceNormal[i + j * ldDisplacement] = Tinv[(i + 6) + (j + 6) * ldTinv];
    displacementToGlobalData[i + j * ldDisplacement] = T[(i + 6) + (j + 6) * ldT];
  }
}

template <typename Cfg>
void TimeAux<Cfg>::extractRotationMatrices(Real<Cfg>** displacementToFaceNormalPtrs,
                                           Real<Cfg>** displacementToGlobalDataPtrs,
                                           Real<Cfg>** TPtrs,
                                           Real<Cfg>** TinvPtrs,
                                           size_t numElements,
                                           void* deviceStream) {
  dim3 block(3, 3, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelextractRotationMatrices<Cfg><<<grid, block, 0, stream>>>(
      displacementToFaceNormalPtrs, displacementToGlobalDataPtrs, TPtrs, TinvPtrs, numElements);
}

template <typename Cfg>
__global__ void kernelInitializeTaylorSeriesForGravitationalBoundary(
    Real<Cfg>** prevCoefficientsPtrs,
    Real<Cfg>** integratedDisplacementNodalPtrs,
    Real<Cfg>** rotatedFaceDisplacementPtrs,
    double deltaTInt,
    size_t numElements) {

  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    auto* prevCoefficients = prevCoefficientsPtrs[elementId];
    auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
    const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

    assert((linearDim<Cfg, seissol::nodal::init::nodes2D<Cfg>>()) <=
           (linearDim<Cfg, seissol::init::rotatedFaceDisplacement<Cfg>>()));

    const int tid = linearidx<Cfg>();
    constexpr auto num2dNodes = linearDim<Cfg, seissol::nodal::init::nodes2D<Cfg>>();
    if (tid < num2dNodes) {
      prevCoefficients[tid] = rotatedFaceDisplacement[tid];
      integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
    }
  }
}

template <typename Cfg>
void TimeAux<Cfg>::initializeTaylorSeriesForGravitationalBoundary(
    Real<Cfg>** prevCoefficientsPtrs,
    Real<Cfg>** integratedDisplacementNodalPtrs,
    Real<Cfg>** rotatedFaceDisplacementPtrs,
    double deltaTInt,
    size_t numElements,
    void* deviceStream) {

  dim3 block = getblock<Cfg>(leadDim<Cfg, seissol::nodal::init::nodes2D<Cfg>>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelInitializeTaylorSeriesForGravitationalBoundary<Cfg>
      <<<grid, block, 0, stream>>>(prevCoefficientsPtrs,
                                   integratedDisplacementNodalPtrs,
                                   rotatedFaceDisplacementPtrs,
                                   deltaTInt,
                                   numElements);
}

template <typename Cfg>
__global__ void kernelComputeInvAcousticImpedance(double* invImpedances,
                                                  double* rhos,
                                                  double* lambdas,
                                                  size_t numElements) {

  size_t index = threadIdx.x + blockIdx.x * blockDim.x;
  if (index < numElements) {
    invImpedances[index] = 1.0 / std::sqrt(lambdas[index] * rhos[index]);
  }
}

template <typename Cfg>
void TimeAux<Cfg>::computeInvAcousticImpedance(
    double* invImpedances, double* rhos, double* lambdas, size_t numElements, void* deviceStream) {
  constexpr size_t blockSize{1024};
  dim3 block(blockSize, 1, 1);
  dim3 grid((numElements + blockSize - 1) / blockSize, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelComputeInvAcousticImpedance<Cfg>
      <<<grid, block, 0, stream>>>(invImpedances, rhos, lambdas, numElements);
}

template <typename Cfg>
__global__ void kernelUpdateRotatedFaceDisplacement(Real<Cfg>** rotatedFaceDisplacementPtrs,
                                                    Real<Cfg>** prevCoefficientsPtrs,
                                                    Real<Cfg>** integratedDisplacementNodalPtrs,
                                                    Real<Cfg>** dofsFaceNodalPtrs,
                                                    double* invImpedances,
                                                    double* rhos,
                                                    double g,
                                                    double factorEvaluated,
                                                    double factorInt,
                                                    size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    constexpr int pIdx = 0;
    constexpr int uIdx = model::MaterialTT<Cfg>::TractionQuantities;
    constexpr auto num2dNodes = linearDim<Cfg, seissol::nodal::init::nodes2D<Cfg>>();

    const int tid = linearidx<Cfg>();
    if (tid < num2dNodes) {

      Real<Cfg>* dofsFaceNodal = dofsFaceNodalPtrs[elementId];
      constexpr auto ldINodal = linearDim<Cfg, seissol::init::INodal<Cfg>>();

      const auto uInside = dofsFaceNodal[tid + (uIdx + 0) * ldINodal];
      const auto vInside = dofsFaceNodal[tid + (uIdx + 1) * ldINodal];
      const auto wInside = dofsFaceNodal[tid + (uIdx + 2) * ldINodal];
      const auto pressureInside = dofsFaceNodal[tid + pIdx * ldINodal];

      Real<Cfg>* prevCoefficients = prevCoefficientsPtrs[elementId];
      const auto rho = rhos[elementId];
      const auto invImpedance = invImpedances[elementId];

      const double curCoeff =
          uInside - invImpedance * (rho * g * prevCoefficients[tid] + pressureInside);
      prevCoefficients[tid] = curCoeff;

      constexpr auto ldFaceDisplacement = linearDim<Cfg, seissol::init::faceDisplacement<Cfg>>();
      static_assert(num2dNodes <= ldFaceDisplacement, "");

      Real<Cfg>* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];
      rotatedFaceDisplacement[tid + 0 * ldFaceDisplacement] += factorEvaluated * curCoeff;
      rotatedFaceDisplacement[tid + 1 * ldFaceDisplacement] += factorEvaluated * vInside;
      rotatedFaceDisplacement[tid + 2 * ldFaceDisplacement] += factorEvaluated * wInside;

      constexpr auto ldIntegratedFaceDisplacement =
          linearDim<Cfg, seissol::init::averageNormalDisplacement<Cfg>>();
      static_assert(num2dNodes <= ldIntegratedFaceDisplacement, "");

      Real<Cfg>* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
      integratedDisplacementNodal[tid] += factorInt * curCoeff;
    }
  }
}

template <typename Cfg>
void TimeAux<Cfg>::updateRotatedFaceDisplacement(Real<Cfg>** rotatedFaceDisplacementPtrs,
                                                 Real<Cfg>** prevCoefficientsPtrs,
                                                 Real<Cfg>** integratedDisplacementNodalPtrs,
                                                 Real<Cfg>** dofsFaceNodalPtrs,
                                                 double* invImpedances,
                                                 double* rhos,
                                                 double g,
                                                 double factorEvaluated,
                                                 double factorInt,
                                                 size_t numElements,
                                                 void* deviceStream) {
  dim3 block = getblock<Cfg>(leadDim<Cfg, seissol::nodal::init::nodes2D<Cfg>>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelUpdateRotatedFaceDisplacement<Cfg>
      <<<grid, block, 0, stream>>>(rotatedFaceDisplacementPtrs,
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

#define SEISSOL_CONFIGITER(cfg) template class TimeAux<cfg>;
#include "ConfigIncludeLinearCK.h"

#define SEISSOL_CONFIGITER(cfg) template class TimeAux<cfg>;
#include "ConfigIncludeSTP.h"

} // namespace seissol::kernels::time::aux
