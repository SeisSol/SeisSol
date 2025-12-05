// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/PlasticityAux.h"

#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Model/Plasticity.h"
#include "Solver/MultipleSimulations.h"

#include <cmath>
#include <cstddef>
#include <type_traits>

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

// NOTE: using c++14 because of cuda@10
namespace seissol {
namespace kernels {
namespace device {
namespace aux {

namespace {

template <typename Cfg>
constexpr int NumStressComponents = model::MaterialTT<Cfg>::TractionQuantities;

template <typename Cfg, typename Tensor>
__forceinline__ __device__ constexpr size_t leadDim() {
  if constexpr (multisim::MultisimEnabled<Cfg>) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

template <typename Cfg>
constexpr auto getblock(int size) {
  if constexpr (multisim::MultisimEnabled<Cfg>) {
    return dim3(multisim::NumSimulations<Cfg>, size);
  } else {
    return dim3(size);
  }
}

template <typename Cfg>
__forceinline__ __device__ auto linearidx() {
  if constexpr (multisim::MultisimEnabled<Cfg>) {
    return threadIdx.y * multisim::NumSimulations<Cfg> + threadIdx.x;
  } else {
    return threadIdx.x;
  }
}

template <typename Cfg>
__forceinline__ __device__ auto simidx() {
  if constexpr (multisim::MultisimEnabled<Cfg>) {
    return threadIdx.x;
  } else {
    return 0;
  }
}

template <typename Cfg>
__forceinline__ __device__ auto validx() {
  if constexpr (multisim::MultisimEnabled<Cfg>) {
    return threadIdx.y;
  } else {
    return threadIdx.x;
  }
}

} // namespace

//--------------------------------------------------------------------------------------------------
template <typename Cfg>
__global__ void
    kernel_adjustDeviatoricTensors(Real<Cfg>** nodalStressTensors,
                                   unsigned* isAdjustableVector,
                                   const seissol::model::PlasticityData<Cfg>* plasticity,
                                   const double oneMinusIntegratingFactor) {
  Real<Cfg>* elementTensors = nodalStressTensors[blockIdx.x];
  Real<Cfg> localStresses[NumStressComponents<Cfg>];

  constexpr auto ElementTensorsColumn = leadDim<init::QStressNodal<Cfg>>();
#pragma unroll
  for (int i = 0; i < NumStressComponents<Cfg>; ++i) {
    localStresses[i] = elementTensors[linearidx<Cfg>() + ElementTensorsColumn * i];
  }

  // 1. Compute the mean stress for each node
  const Real<Cfg> meanStress = (localStresses[0] + localStresses[1] + localStresses[2]) / 3.0f;

// 2. Compute deviatoric stress tensor
#pragma unroll
  for (int i = 0; i < 3; ++i) {
    localStresses[i] -= meanStress;
  }

  // 3. Compute the second invariant for each node
  Real<Cfg> tau = 0.5 * (localStresses[0] * localStresses[0] + localStresses[1] * localStresses[1] +
                         localStresses[2] * localStresses[2]);
  tau += (localStresses[3] * localStresses[3] + localStresses[4] * localStresses[4] +
          localStresses[5] * localStresses[5]);
  tau = std::sqrt(tau);

  // 4. Compute the plasticity criteria
  const Real<Cfg> cohesionTimesCosAngularFriction =
      plasticity[blockIdx.x].cohesionTimesCosAngularFriction[simidx<Cfg>()];
  const Real<Cfg> sinAngularFriction = plasticity[blockIdx.x].sinAngularFriction[simidx<Cfg>()];
  Real<Cfg> taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
  taulim = std::max(static_cast<Real<Cfg>>(0.0), taulim);

  __shared__ unsigned isAdjusted;
  if (linearidx<Cfg>() == 0) {
    isAdjusted = static_cast<unsigned>(false);
  }
  __syncthreads();

  // 5. Compute the yield factor
  Real<Cfg> factor = 0.0;
  if (tau > taulim) {
    isAdjusted = static_cast<unsigned>(true);
    factor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
  }

  // 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
  __syncthreads();
  if (isAdjusted) {
#pragma unroll
    for (int i = 0; i < NumStressComponents<Cfg>; ++i) {
      elementTensors[linearidx<Cfg>() + ElementTensorsColumn * i] = localStresses[i] * factor;
    }
  }
  if (linearidx<Cfg>() == 0) {
    isAdjustableVector[blockIdx.x] = isAdjusted;
  }
}

template <typename Cfg>
void Plasticity<Cfg>::adjustDeviatoricTensors(Real<Cfg>** nodalStressTensors,
                                              unsigned* isAdjustableVector,
                                              const seissol::model::PlasticityData<Cfg>* plasticity,
                                              const double oneMinusIntegratingFactor,
                                              const size_t numElements,
                                              void* streamPtr) {
  constexpr unsigned NumNodes = tensor::QStressNodal<Cfg>::Shape[multisim::BasisDim<Cfg>];
  const auto block = getblock<Cfg>(NumNodes);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_adjustDeviatoricTensors<Cfg><<<grid, block, 0, stream>>>(
      nodalStressTensors, isAdjustableVector, plasticity, oneMinusIntegratingFactor);
}

//--------------------------------------------------------------------------------------------------
template <typename Cfg>
__global__ void kernel_computePstrains(Real<Cfg>** pstrains,
                                       const seissol::model::PlasticityData<Cfg>* plasticityData,
                                       Real<Cfg>** dofs,
                                       Real<Cfg>** prevDofs,
                                       Real<Cfg>** dUdTpstrain,
                                       double tV,
                                       double oneMinusIntegratingFactor,
                                       double timeStepWidth,
                                       const unsigned* isAdjustableVector) {
  if (isAdjustableVector[blockIdx.x]) {
    Real<Cfg>* localDofs = dofs[blockIdx.x];
    Real<Cfg>* localPrevDofs = prevDofs[blockIdx.x];
    const seissol::model::PlasticityData<Cfg>* localData = &plasticityData[blockIdx.x];
    Real<Cfg>* localPstrain = pstrains[blockIdx.x];
    Real<Cfg>* localDuDtPstrain = dUdTpstrain[blockIdx.x];

#pragma unroll
    for (int i = 0; i < NumStressComponents<Cfg>; ++i) {
      const int q = linearidx<Cfg>() + i * leadDim<init::Q<Cfg>>();
      const Real<Cfg> factor = localData->mufactor / (tV * oneMinusIntegratingFactor);
      const Real<Cfg> nodeDuDtPstrain = factor * (localPrevDofs[q] - localDofs[q]);

      static_assert(leadDim<init::QStress<Cfg>>() == leadDim<init::Q<Cfg>>(), "");
      localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
      localDuDtPstrain[q] = nodeDuDtPstrain;
    }
  }
}

template <typename Cfg>
void Plasticity<Cfg>::computePstrains(Real<Cfg>** pstrains,
                                      const seissol::model::PlasticityData<Cfg>* plasticityData,
                                      Real<Cfg>** dofs,
                                      Real<Cfg>** prevDofs,
                                      Real<Cfg>** dUdTpstrain,
                                      double tV,
                                      double oneMinusIntegratingFactor,
                                      double timeStepWidth,
                                      unsigned* isAdjustableVector,
                                      size_t numElements,
                                      void* streamPtr) {
  constexpr unsigned NumNodes = tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg>];
  const dim3 block = getblock<Cfg>(NumNodes);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_computePstrains<Cfg><<<grid, block, 0, stream>>>(pstrains,
                                                          plasticityData,
                                                          dofs,
                                                          prevDofs,
                                                          dUdTpstrain,
                                                          tV,
                                                          oneMinusIntegratingFactor,
                                                          timeStepWidth,
                                                          isAdjustableVector);
}

//--------------------------------------------------------------------------------------------------
template <typename Cfg>
__global__ void kernel_updateQEtaNodal(Real<Cfg>** qEtaNodalPtrs,
                                       Real<Cfg>** qStressNodalPtrs,
                                       double timeStepWidth,
                                       const unsigned* isAdjustableVector) {
  if (isAdjustableVector[blockIdx.x]) {
    const size_t tid = linearidx<Cfg>();
    Real<Cfg>* localQEtaNodal = qEtaNodalPtrs[blockIdx.x];
    Real<Cfg>* localQStressNodal = qStressNodalPtrs[blockIdx.x];
    Real<Cfg> factor{0.0};

    constexpr auto Ld = leadDim<init::QStressNodal<Cfg>>();
#pragma unroll
    for (int i = 0; i < NumStressComponents<Cfg>; ++i) {
      factor += localQStressNodal[tid + i * Ld] * localQStressNodal[tid + i * Ld];
    }

    localQEtaNodal[tid] = std::max(static_cast<Real<Cfg>>(0.0), localQEtaNodal[tid]) +
                          timeStepWidth * std::sqrt(static_cast<Real<Cfg>>(0.5) * factor);
  }
}

template <typename Cfg>
void Plasticity<Cfg>::updateQEtaNodal(Real<Cfg>** qEtaNodalPtrs,
                                      Real<Cfg>** qStressNodalPtrs,
                                      double timeStepWidth,
                                      unsigned* isAdjustableVector,
                                      size_t numElements,
                                      void* streamPtr) {
  const dim3 block = getblock<Cfg>(tensor::QStressNodal<Cfg>::Shape[multisim::BasisDim<Cfg>]);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_updateQEtaNodal<Cfg><<<grid, block, 0, stream>>>(
      qEtaNodalPtrs, qStressNodalPtrs, timeStepWidth, isAdjustableVector);
}

} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol
