// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/PlasticityAux.h"
#include "Kernels/Precision.h"
#include "Model/Plasticity.h"
#include "tensor.h"
#include <cmath>
#include <cstddef>
#include <driver_types.h>
#include <init.h>
#include <type_traits>
#include <vector_types.h>

#include <Solver/MultipleSimulations.h>

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
namespace plasticity {

template <typename Tensor>
__forceinline__ __device__ constexpr size_t leadDim() {
  if constexpr (multisim::MultisimEnabled) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

constexpr auto getblock(int size) {
  if constexpr (multisim::MultisimEnabled) {
    return dim3(multisim::NumSimulations, size);
  } else {
    return dim3(size);
  }
}

__forceinline__ __device__ auto linearidx() {
  if constexpr (multisim::MultisimEnabled) {
    return threadIdx.y * multisim::NumSimulations + threadIdx.x;
  } else {
    return threadIdx.x;
  }
}

__forceinline__ __device__ auto simidx() {
  if constexpr (multisim::MultisimEnabled) {
    return threadIdx.x;
  } else {
    return 0;
  }
}

__forceinline__ __device__ auto validx() {
  if constexpr (multisim::MultisimEnabled) {
    return threadIdx.y;
  } else {
    return threadIdx.x;
  }
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_adjustDeviatoricTensors(real** nodalStressTensors,
                                               unsigned* isAdjustableVector,
                                               const seissol::model::PlasticityData* plasticity,
                                               const double oneMinusIntegratingFactor) {
  real* elementTensors = nodalStressTensors[blockIdx.x];
  real localStresses[NUM_STRESS_COMPONENTS];

  constexpr auto ElementTensorsColumn = leadDim<init::QStressNodal>();
#pragma unroll
  for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
    localStresses[i] = elementTensors[linearidx() + ElementTensorsColumn * i];
  }

  // 1. Compute the mean stress for each node
  const real meanStress = (localStresses[0] + localStresses[1] + localStresses[2]) / 3.0f;

// 2. Compute deviatoric stress tensor
#pragma unroll
  for (int i = 0; i < 3; ++i) {
    localStresses[i] -= meanStress;
  }

  // 3. Compute the second invariant for each node
  real tau = 0.5 * (localStresses[0] * localStresses[0] + localStresses[1] * localStresses[1] +
                    localStresses[2] * localStresses[2]);
  tau += (localStresses[3] * localStresses[3] + localStresses[4] * localStresses[4] +
          localStresses[5] * localStresses[5]);
  tau = std::sqrt(tau);

  // 4. Compute the plasticity criteria
  const real cohesionTimesCosAngularFriction =
      plasticity[blockIdx.x].cohesionTimesCosAngularFriction;
  const real sinAngularFriction = plasticity[blockIdx.x].sinAngularFriction;
  real taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
  taulim = std::max(static_cast<real>(0.0), taulim);

  __shared__ unsigned isAdjusted;
  if (validx() == 0 && simidx() == 0) {
    isAdjusted = static_cast<unsigned>(false);
  }
  __syncthreads();

  // 5. Compute the yield factor
  real factor = 0.0;
  if (tau > taulim) {
    isAdjusted = static_cast<unsigned>(true);
    factor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
  }

  // 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
  __syncthreads();
  if (isAdjusted) {
#pragma unroll
    for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
      elementTensors[linearidx() + ElementTensorsColumn * i] = localStresses[i] * factor;
    }
  }
  if (validx() == 0 && simidx() == 0) {
    isAdjustableVector[blockIdx.x] = isAdjusted;
  }
}

void adjustDeviatoricTensors(real** nodalStressTensors,
                             unsigned* isAdjustableVector,
                             const seissol::model::PlasticityData* plasticity,
                             const double oneMinusIntegratingFactor,
                             const size_t numElements,
                             void* streamPtr) {
  constexpr unsigned NumNodes = tensor::QStressNodal::Shape[multisim::BasisFunctionDimension];
  const auto block = getblock(NumNodes);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_adjustDeviatoricTensors<<<grid, block, 0, stream>>>(
      nodalStressTensors, isAdjustableVector, plasticity, oneMinusIntegratingFactor);
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_adjustPointers(real* qEtaNodal,
                                      real** qEtaNodalPtrs,
                                      real* qEtaModal,
                                      real** qEtaModalPtrs,
                                      real* dUdTpstrain,
                                      real** dUdTpstrainPtrs,
                                      size_t numElements) {

  const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < numElements) {
    qEtaNodalPtrs[tid] = &qEtaNodal[tensor::QEtaNodal::Size * tid];
    qEtaModalPtrs[tid] = &qEtaModal[tensor::QEtaModal::Size * tid];
    dUdTpstrainPtrs[tid] = &dUdTpstrain[tensor::QStressNodal::Size * tid];
  }
}

void adjustPointers(real* qEtaNodal,
                    real** qEtaNodalPtrs,
                    real* qEtaModal,
                    real** qEtaModalPtrs,
                    real* dUdTpstrain,
                    real** dUdTpstrainPtrs,
                    size_t numElements,
                    void* streamPtr) {
  const dim3 block(1024, 1, 1);
  const size_t numBlocks = (numElements + block.x - 1) / block.x;
  const dim3 grid(numBlocks, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_adjustPointers<<<grid, block, 0, stream>>>(qEtaNodal,
                                                    qEtaNodalPtrs,
                                                    qEtaModal,
                                                    qEtaModalPtrs,
                                                    dUdTpstrain,
                                                    dUdTpstrainPtrs,
                                                    numElements);
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_computePstrains(real** pstrains,
                                       const seissol::model::PlasticityData* plasticityData,
                                       real** dofs,
                                       real* prevDofs,
                                       real** dUdTpstrain,
                                       double tV,
                                       double oneMinusIntegratingFactor,
                                       double timeStepWidth,
                                       const unsigned* isAdjustableVector) {
  if (isAdjustableVector[blockIdx.x]) {
    real* localDofs = dofs[blockIdx.x];
    real* localPrevDofs = &prevDofs[tensor::Q::Size * blockIdx.x];
    const seissol::model::PlasticityData* localData = &plasticityData[blockIdx.x];
    real* localPstrain = pstrains[blockIdx.x];
    real* localDuDtPstrain = dUdTpstrain[blockIdx.x];

#pragma unroll
    for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
      const int q = linearidx() + i * leadDim<init::Q>();
      const real factor = localData->mufactor / (tV * oneMinusIntegratingFactor);
      const real nodeDuDtPstrain = factor * (localPrevDofs[q] - localDofs[q]);

      static_assert(leadDim<init::QStress>() == leadDim<init::Q>(), "");
      localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
      localDuDtPstrain[q] = nodeDuDtPstrain;
    }
  }
}

void computePstrains(real** pstrains,
                     const seissol::model::PlasticityData* plasticityData,
                     real** dofs,
                     real* prevDofs,
                     real** dUdTpstrain,
                     double tV,
                     double oneMinusIntegratingFactor,
                     double timeStepWidth,
                     unsigned* isAdjustableVector,
                     size_t numElements,
                     void* streamPtr) {
  constexpr unsigned NumNodes = tensor::Q::Shape[multisim::BasisFunctionDimension];
  const dim3 block = getblock(NumNodes);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_computePstrains<<<grid, block, 0, stream>>>(pstrains,
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
__global__ void kernel_pstrainToQEtaModal(real** pstrains,
                                          real** qEtaModalPtrs,
                                          const unsigned* isAdjustableVector) {
  static_assert(tensor::QEtaModal::Size == leadDim<init::QStressNodal>(), "");

  if (isAdjustableVector[blockIdx.x]) {
    real* localQEtaModal = qEtaModalPtrs[blockIdx.x];
    real* localPstrain = pstrains[blockIdx.x];
    for (int i = threadIdx.x; i < tensor::QEtaModal::Size; i += 1024) {
      localQEtaModal[i] = localPstrain[NUM_STRESS_COMPONENTS * leadDim<init::QStressNodal>() + i];
    }
  }
}

void pstrainToQEtaModal(real** pstrains,
                        real** qEtaModalPtrs,
                        unsigned* isAdjustableVector,
                        size_t numElements,
                        void* streamPtr) {
  const dim3 block(1024, 1, 1);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_pstrainToQEtaModal<<<grid, block, 0, stream>>>(
      pstrains, qEtaModalPtrs, isAdjustableVector);
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_qEtaModalToPstrain(real** qEtaModalPtrs,
                                          real** pstrains,
                                          const unsigned* isAdjustableVector) {
  static_assert(tensor::QEtaModal::Size == leadDim<init::QStressNodal>(), "");

  if (isAdjustableVector[blockIdx.x]) {
    real* localQEtaModal = qEtaModalPtrs[blockIdx.x];
    real* localPstrain = pstrains[blockIdx.x];
    for (int i = threadIdx.x; i < tensor::QEtaModal::Size; i += 1024) {
      localPstrain[NUM_STRESS_COMPONENTS * leadDim<init::QStressNodal>() + i] = localQEtaModal[i];
    }
  }
}

void qEtaModalToPstrain(real** qEtaModalPtrs,
                        real** pstrains,
                        unsigned* isAdjustableVector,
                        size_t numElements,
                        void* streamPtr) {
  const dim3 block(1024, 1, 1);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_qEtaModalToPstrain<<<grid, block, 0, stream>>>(
      qEtaModalPtrs, pstrains, isAdjustableVector);
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_updateQEtaNodal(real** qEtaNodalPtrs,
                                       real** qStressNodalPtrs,
                                       double timeStepWidth,
                                       const unsigned* isAdjustableVector) {
  if (isAdjustableVector[blockIdx.x]) {
    const size_t tid = linearidx();
    real* localQEtaNodal = qEtaNodalPtrs[blockIdx.x];
    real* localQStressNodal = qStressNodalPtrs[blockIdx.x];
    real factor{0.0};

    constexpr auto Ld = leadDim<init::QStressNodal>();
#pragma unroll
    for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
      factor += localQStressNodal[tid + i * Ld] * localQStressNodal[tid + i * Ld];
    }

    localQEtaNodal[tid] = std::max(static_cast<real>(0.0), localQEtaNodal[tid]) +
                          timeStepWidth * std::sqrt(static_cast<real>(0.5) * factor);
  }
}

void updateQEtaNodal(real** qEtaNodalPtrs,
                     real** qStressNodalPtrs,
                     double timeStepWidth,
                     unsigned* isAdjustableVector,
                     size_t numElements,
                     void* streamPtr) {
  const dim3 block = getblock(tensor::QStressNodal::Shape[multisim::BasisFunctionDimension]);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_updateQEtaNodal<<<grid, block, 0, stream>>>(
      qEtaNodalPtrs, qStressNodalPtrs, timeStepWidth, isAdjustableVector);
}

} // namespace plasticity
} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol
