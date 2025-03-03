// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
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


// NOTE: using c++14 because of cuda@10
namespace seissol {
namespace kernels {
namespace device {
namespace aux {
namespace plasticity {

template<typename T>
__forceinline__ __device__ typename std::enable_if<std::is_floating_point<T>::value, T>::type
squareRoot(T x) {
  return std::is_same<T, double>::value ? sqrt(x) : sqrtf(x);
}

template<typename T>
__forceinline__ __device__ typename std::enable_if<std::is_floating_point<T>::value, T>::type
maxValue(T x, T y) {
  return std::is_same<T, double>::value ? fmax(x, y) : fmaxf(x, y);
}

template<typename Tensor>
__forceinline__  __device__
constexpr size_t leadDim() {
  return Tensor::Stop[0] - Tensor::Start[0];
}


//--------------------------------------------------------------------------------------------------
__global__ void kernel_adjustDeviatoricTensors(real **nodalStressTensors,
                                               unsigned *isAdjustableVector,
                                               const seissol::model::PlasticityData *plasticity,
                                               const double oneMinusIntegratingFactor) {
  real *elementTensors = nodalStressTensors[blockIdx.x];
  real localStresses[NUM_STRESS_COMPONENTS];


  constexpr auto ElementTensorsColumn = leadDim<init::QStressNodal>();
  #pragma unroll
  for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
    localStresses[i] = elementTensors[threadIdx.x + ElementTensorsColumn * i];
  }

  // 1. Compute the mean stress for each node
  real const meanStress = (localStresses[0] + localStresses[1] + localStresses[2]) / 3.0f;

  // 2. Compute deviatoric stress tensor
  #pragma unroll
  for (int i = 0; i < 3; ++i) {
    localStresses[i] -= meanStress;
  }

  // 3. Compute the second invariant for each node
  real tau = 0.5 * (localStresses[0] * localStresses[0] +
                    localStresses[1] * localStresses[1] +
                    localStresses[2] * localStresses[2]);
  tau += (localStresses[3] * localStresses[3] +
          localStresses[4] * localStresses[4] +
          localStresses[5] * localStresses[5]);
  tau = squareRoot(tau);

  // 4. Compute the plasticity criteria
  const real cohesionTimesCosAngularFriction = plasticity[blockIdx.x].cohesionTimesCosAngularFriction;
  const real sinAngularFriction = plasticity[blockIdx.x].sinAngularFriction;
  real taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
  taulim = maxValue(static_cast<real>(0.0), taulim);

  __shared__ unsigned isAdjusted;
  if (threadIdx.x == 0) { isAdjusted = static_cast<unsigned>(false); }
  __syncthreads();

  // 5. Compute the yield factor
  real factor = 0.0;
  if (tau > taulim) {
    isAdjusted = static_cast<unsigned >(true);
    factor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
  }

  // 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
  __syncthreads();
  if (isAdjusted) {
    #pragma unroll
    for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
      elementTensors[threadIdx.x + ElementTensorsColumn * i] = localStresses[i] * factor;
    }
  }
  if (threadIdx.x == 0) {
    isAdjustableVector[blockIdx.x] = isAdjusted;
  }
}

void adjustDeviatoricTensors(real **nodalStressTensors,
                             unsigned *isAdjustableVector,
                             const seissol::model::PlasticityData *plasticity,
                             const double oneMinusIntegratingFactor,
                             const size_t numElements,
                             void *streamPtr) {
  constexpr unsigned NumNodes = tensor::QStressNodal::Shape[0];
  dim3 const block(NumNodes, 1, 1);
  dim3 const grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_adjustDeviatoricTensors<<<grid, block, 0, stream>>>(nodalStressTensors,
                                                             isAdjustableVector,
                                                             plasticity,
                                                             oneMinusIntegratingFactor);
}


//--------------------------------------------------------------------------------------------------
__global__ void kernel_adjustPointers(real *qEtaNodal,
                                      real **qEtaNodalPtrs,
                                      real *qEtaModal,
                                      real **qEtaModalPtrs,
                                      real *dUdTpstrain,
                                      real **dUdTpstrainPtrs,
                                      size_t numElements) {

  size_t const tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < numElements) {
    qEtaNodalPtrs[tid] = &qEtaNodal[tensor::QEtaNodal::Size * tid];
    qEtaModalPtrs[tid] = &qEtaModal[tensor::QEtaModal::Size * tid];
    dUdTpstrainPtrs[tid] = &dUdTpstrain[tensor::QStressNodal::Size * tid];
  }
}

void adjustPointers(real *qEtaNodal,
                    real **qEtaNodalPtrs,
                    real *qEtaModal,
                    real **qEtaModalPtrs,
                    real *dUdTpstrain,
                    real **dUdTpstrainPtrs,
                    size_t numElements,
                    void *streamPtr) {
  dim3 const block(128, 1, 1);
  size_t const numBlocks = (numElements + block.x - 1) / block.x;
  dim3 const grid(numBlocks, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_adjustPointers<<<grid, block, 0, stream>>>(qEtaNodal,
                                                    qEtaNodalPtrs,
                                                    qEtaModal,
                                                    qEtaModalPtrs,
                                                    dUdTpstrain,
                                                    dUdTpstrainPtrs,
                                                    numElements);
}


//--------------------------------------------------------------------------------------------------
__global__ void kernel_computePstrains(real **pstrains,
                                       const seissol::model::PlasticityData *plasticityData,
                                       real **dofs,
                                       real *prevDofs,
                                       real **dUdTpstrain,
                                       double tV,
                                       double oneMinusIntegratingFactor,
                                       double timeStepWidth,
                                       const unsigned *isAdjustableVector) {
  if (isAdjustableVector[blockIdx.x]) {
    real *localDofs = dofs[blockIdx.x];
    real *localPrevDofs = &prevDofs[tensor::Q::Size * blockIdx.x];
    const seissol::model::PlasticityData *localData = &plasticityData[blockIdx.x];
    real *localPstrain = pstrains[blockIdx.x];
    real *localDuDtPstrain = dUdTpstrain[blockIdx.x];

    #pragma unroll
    for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
      int const q = threadIdx.x + i * leadDim<init::Q>();
      real const factor = localData->mufactor / (tV * oneMinusIntegratingFactor);
      real const nodeDuDtPstrain = factor * (localPrevDofs[q] - localDofs[q]);

      static_assert(leadDim<init::QStress>() == leadDim<init::Q>(), "");
      localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
      localDuDtPstrain[q] = nodeDuDtPstrain;

    }
  }
}

void computePstrains(real **pstrains,
                     const seissol::model::PlasticityData *plasticityData,
                     real **dofs,
                     real *prevDofs,
                     real **dUdTpstrain,
                     double tV,
                     double oneMinusIntegratingFactor,
                     double timeStepWidth,
                     unsigned *isAdjustableVector,
                     size_t numElements,
                     void *streamPtr) {
  constexpr unsigned NumNodes = tensor::Q::Shape[0];
  dim3 const block(NumNodes, 1, 1);
  dim3 const grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
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
__global__ void kernel_pstrainToQEtaModal(real **pstrains,
                                          real **qEtaModalPtrs,
                                          const unsigned *isAdjustableVector) {
    static_assert(tensor::QEtaModal::Size == leadDim<init::QStressNodal>(), "");

    if (isAdjustableVector[blockIdx.x]) {
    real *localQEtaModal = qEtaModalPtrs[blockIdx.x];
    real *localPstrain = pstrains[blockIdx.x];
    localQEtaModal[threadIdx.x] = localPstrain[NUM_STRESS_COMPONENTS * leadDim<init::QStressNodal>() + threadIdx.x];
  }
}

void pstrainToQEtaModal(real **pstrains,
                       real **qEtaModalPtrs,
                       unsigned *isAdjustableVector,
                       size_t numElements,
                       void *streamPtr) {
  dim3 const block(tensor::QEtaModal::Size, 1, 1);
  dim3 const grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_pstrainToQEtaModal<<<grid, block, 0, stream>>>(pstrains, qEtaModalPtrs, isAdjustableVector);
}


//--------------------------------------------------------------------------------------------------
__global__ void kernel_qEtaModalToPstrain(real **qEtaModalPtrs,
                                          real **pstrains,
                                          const unsigned *isAdjustableVector) {
    static_assert(tensor::QEtaModal::Size == leadDim<init::QStressNodal>(), "");

    if (isAdjustableVector[blockIdx.x]) {
    real *localQEtaModal = qEtaModalPtrs[blockIdx.x];
    real *localPstrain = pstrains[blockIdx.x];
    localPstrain[NUM_STRESS_COMPONENTS * leadDim<init::QStressNodal>() + threadIdx.x] = localQEtaModal[threadIdx.x];
  }
}

void qEtaModalToPstrain(real **qEtaModalPtrs,
                        real **pstrains,
                        unsigned *isAdjustableVector,
                        size_t numElements,
                        void *streamPtr) {
  dim3 const block(tensor::QEtaModal::Size, 1, 1);
  dim3 const grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_qEtaModalToPstrain<<<grid, block, 0, stream>>>(qEtaModalPtrs, pstrains, isAdjustableVector);
}


//--------------------------------------------------------------------------------------------------
__global__ void kernel_updateQEtaNodal(real **qEtaNodalPtrs,
                                       real **qStressNodalPtrs,
                                       double timeStepWidth,
                                       const unsigned *isAdjustableVector) {
  if (isAdjustableVector[blockIdx.x]) {
    size_t const tid = threadIdx.x;
    real *localQEtaNodal = qEtaNodalPtrs[blockIdx.x];
    real *localQStressNodal = qStressNodalPtrs[blockIdx.x];
    real factor{0.0};

    constexpr auto Ld = leadDim<init::QStressNodal>();
    #pragma unroll
    for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
      factor += localQStressNodal[tid + i * Ld] * localQStressNodal[tid + i * Ld];
    }

    localQEtaNodal[tid] = maxValue(static_cast<real>(0.0), localQEtaNodal[tid])
        + timeStepWidth * squareRoot(static_cast<real>(0.5) * factor);
  }
}

void updateQEtaNodal(real **qEtaNodalPtrs,
                     real **qStressNodalPtrs,
                     double timeStepWidth,
                     unsigned *isAdjustableVector,
                     size_t numElements,
                     void *streamPtr) {
  dim3 const block(tensor::QStressNodal::Shape[0], 1, 1);
  dim3 const grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_updateQEtaNodal<<<grid, block, 0, stream>>>(qEtaNodalPtrs, qStressNodalPtrs, timeStepWidth, isAdjustableVector);
}

} // namespace plasticity
} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol

