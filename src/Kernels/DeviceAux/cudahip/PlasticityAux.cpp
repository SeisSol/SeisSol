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
namespace plasticity {

template <typename Tensor>
__forceinline__ __device__ __host__ constexpr size_t leadDim() {
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

__global__ void
    kernel_plasticityNonlinear(real** __restrict nodalStressTensors,
                               real** __restrict prevNodal,
                               real** __restrict pstrainPtr,
                               unsigned* __restrict isAdjustableVector,
                               std::size_t* __restrict yieldCounter,
                               const seissol::model::PlasticityData* __restrict plasticity,
                               double oneMinusIntegratingFactor,
                               double tV,
                               double timeStepWidth) {
  real* __restrict qStressNodal = nodalStressTensors[blockIdx.x];
  real localStresses[NumStressComponents];

  constexpr auto ElementTensorsColumn = leadDim<init::QStressNodal>();
#pragma unroll
  for (int i = 0; i < NumStressComponents; ++i) {
    localStresses[i] = qStressNodal[linearidx() + ElementTensorsColumn * i];
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
      plasticity[blockIdx.x].cohesionTimesCosAngularFriction[simidx()];
  const real sinAngularFriction = plasticity[blockIdx.x].sinAngularFriction[simidx()];
  const real taulim = std::max(static_cast<real>(0.0),
                               cohesionTimesCosAngularFriction - meanStress * sinAngularFriction);

  __shared__ bool isAdjusted;
  if (linearidx() == 0) {
    isAdjusted = false;
  }
  __syncthreads();

  // 5. Compute the yield factor
  real yieldfactor = 0.0;
  if (tau > taulim) {
    isAdjusted = true;
    yieldfactor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
  }

  // 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
  __syncthreads();
  if (isAdjusted) {
    const real factor = plasticity[blockIdx.x].mufactor / (tV * oneMinusIntegratingFactor);

    const real* __restrict localPrevNodal = prevNodal[blockIdx.x];
    real* __restrict eta = pstrainPtr[blockIdx.x] + tensor::QStressNodal::size();
    real* __restrict localPstrain = pstrainPtr[blockIdx.x];

    real dudtUpdate = 0;

#pragma unroll
    for (int i = 0; i < NumStressComponents; ++i) {
      const int q = linearidx() + ElementTensorsColumn * i;

      const auto updatedStressNodal = localStresses[i] * yieldfactor;

      const real nodeDuDtPstrain = factor * (localPrevNodal[q] - updatedStressNodal);

      localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
      qStressNodal[q] = updatedStressNodal;

      dudtUpdate += nodeDuDtPstrain * nodeDuDtPstrain;
    }

    eta[linearidx()] += timeStepWidth * std::sqrt(static_cast<real>(0.5) * dudtUpdate);

    // update the FLOPs that we've been here
    atomicAdd(reinterpret_cast<unsigned long long*>(yieldCounter), 1);
  }
  if (linearidx() == 0) {
    isAdjustableVector[blockIdx.x] = isAdjusted;
  }
}

void plasticityNonlinear(real** __restrict nodalStressTensors,
                         real** __restrict prevNodal,
                         real** __restrict pstrainPtr,
                         unsigned* __restrict isAdjustableVector,
                         std::size_t* __restrict yieldCounter,
                         const seissol::model::PlasticityData* __restrict plasticity,
                         double oneMinusIntegratingFactor,
                         double tV,
                         double timeStepWidth,
                         size_t numElements,
                         void* streamPtr) {
  // use Stop/Start to include padding (and possibly avoid masked warps/wavefronts)
  constexpr unsigned NumNodes = init::QStressNodal::Stop[multisim::BasisFunctionDimension] -
                                init::QStressNodal::Start[multisim::BasisFunctionDimension];
  const auto block = getblock(NumNodes);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_plasticityNonlinear<<<grid, block, 0, stream>>>(nodalStressTensors,
                                                         prevNodal,
                                                         pstrainPtr,
                                                         isAdjustableVector,
                                                         yieldCounter,
                                                         plasticity,
                                                         oneMinusIntegratingFactor,
                                                         tV,
                                                         timeStepWidth);
}

} // namespace plasticity
} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol
