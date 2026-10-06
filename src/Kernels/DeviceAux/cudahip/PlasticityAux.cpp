// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/PlasticityAux.h"

#include "Config.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Model/Plasticity.h"
#include "Solver/MultipleSimulations.h"

#include <cmath>
#include <cstddef>

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

namespace seissol::kernels::device::aux::plasticity {

template <typename Cfg>
constexpr auto getblock(int size) {
  if constexpr (Cfg::NumSimulations > 1) {
    return dim3(Cfg::NumSimulations, size);
  } else {
    return dim3(size);
  }
}

template <typename Cfg>
__forceinline__ __device__ auto linearidx() {
  if constexpr (Cfg::NumSimulations > 1) {
    return threadIdx.y * Cfg::NumSimulations + threadIdx.x;
  } else {
    return threadIdx.x;
  }
}

template <typename Cfg>
__global__ void
    kernel_plasticityNonlinear(Real<Cfg>** __restrict nodalStressTensors,
                               Real<Cfg>** __restrict pstrainPtr,
                               unsigned* __restrict isAdjustableVector,
                               std::size_t* __restrict yieldCounter,
                               const seissol::model::PlasticityData<Cfg>* __restrict plasticity,
                               Real<Cfg> oneMinusIntegratingFactor,
                               Real<Cfg> tV,
                               Real<Cfg> timeStepWidth) {

  using real = Real<Cfg>;

  real* __restrict qStressNodal = nodalStressTensors[blockIdx.x];
  real localStresses[NumStressComponents<Cfg>];

  constexpr auto ElementTensorsColumn = multisim::linearDim<Cfg, init::QStressNodal<Cfg>>();
#pragma unroll
  for (int i = 0; i < NumStressComponents<Cfg>; ++i) {
    localStresses[i] = qStressNodal[linearidx<Cfg>() + ElementTensorsColumn * i];
  }

  // 1. Compute the mean stress for each node
  const real meanStress =
      (localStresses[0] + localStresses[1] + localStresses[2]) / static_cast<real>(3);

// 2. Compute deviatoric stress tensor
#pragma unroll
  for (int i = 0; i < 3; ++i) {
    localStresses[i] -= meanStress;
  }

  // 3. Compute the second invariant for each node
  real tau = static_cast<real>(0.5) *
             (localStresses[0] * localStresses[0] + localStresses[1] * localStresses[1] +
              localStresses[2] * localStresses[2]);
  tau += (localStresses[3] * localStresses[3] + localStresses[4] * localStresses[4] +
          localStresses[5] * localStresses[5]);
  tau = std::sqrt(tau);

  // 4. Compute the plasticity criteria
  const real cohesionTimesCosAngularFriction =
      plasticity[blockIdx.x].cohesionTimesCosAngularFriction[linearidx<Cfg>()];
  const real sinAngularFriction = plasticity[blockIdx.x].sinAngularFriction[linearidx<Cfg>()];
  const real taulim = std::max(static_cast<real>(0.0),
                               cohesionTimesCosAngularFriction - meanStress * sinAngularFriction);

  __shared__ bool isAdjusted;
  if (linearidx<Cfg>() == 0) {
    isAdjusted = false;
  }
  __syncthreads();

  // 5. Compute the yield factor
  real yieldfactor{};
  if (tau > taulim) {
    isAdjusted = true;
    yieldfactor = ((taulim / tau) - static_cast<real>(1.0)) * oneMinusIntegratingFactor;
  }

  // 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
  __syncthreads();
  if (isAdjusted) {
    const real factor = plasticity[blockIdx.x].mufactor / (tV * oneMinusIntegratingFactor);

    real* __restrict eta = pstrainPtr[blockIdx.x] + tensor::QStressNodal<Cfg>::size();
    real* __restrict localPstrain = pstrainPtr[blockIdx.x];

    real dudtUpdate = 0;

#pragma unroll
    for (int i = 0; i < NumStressComponents<Cfg>; ++i) {
      const int q = linearidx<Cfg>() + ElementTensorsColumn * i;

      const auto updatedStressNodal = localStresses[i] * yieldfactor;

      const real nodeDuDtPstrain = -factor * updatedStressNodal;

      localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
      qStressNodal[q] = updatedStressNodal;

      dudtUpdate += nodeDuDtPstrain * nodeDuDtPstrain;
    }

    eta[linearidx<Cfg>()] += timeStepWidth * std::sqrt(static_cast<real>(0.5) * dudtUpdate);

    // update the FLOPs that we've been here
    // (there's no atomicAdd for unsigned long / sometimes size_t, so take one of the other ones)
    static_assert(sizeof(unsigned long long) == sizeof(std::size_t));
    atomicAdd(reinterpret_cast<unsigned long long*>(yieldCounter), 1);
  }
  if (linearidx<Cfg>() == 0) {
    isAdjustableVector[blockIdx.x] = isAdjusted;
  }
}

template <typename Cfg>
void plasticityNonlinear(Real<Cfg>** __restrict nodalStressTensors,
                         Real<Cfg>** __restrict pstrainPtr,
                         unsigned* __restrict isAdjustableVector,
                         std::size_t* __restrict yieldCounter,
                         const seissol::model::PlasticityData<Cfg>* __restrict plasticity,
                         Real<Cfg> oneMinusIntegratingFactor,
                         Real<Cfg> tV,
                         Real<Cfg> timeStepWidth,
                         size_t numElements,
                         void* streamPtr) {
  // use Stop/Start to include padding (and possibly avoid masked warps/wavefronts)
  constexpr unsigned NumNodes = init::QStressNodal<Cfg>::Stop[multisim::BasisDim<Cfg>] -
                                init::QStressNodal<Cfg>::Start[multisim::BasisDim<Cfg>];
  const auto block = getblock<Cfg>(NumNodes);
  const dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(streamPtr);
  kernel_plasticityNonlinear<<<grid, block, 0, stream>>>(nodalStressTensors,
                                                         pstrainPtr,
                                                         isAdjustableVector,
                                                         yieldCounter,
                                                         plasticity,
                                                         oneMinusIntegratingFactor,
                                                         tV,
                                                         timeStepWidth);
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template void plasticityNonlinear<Cfg>(                                                          \
      Real<Cfg> * * __restrict nodalStressTensors,                                                 \
      Real<Cfg> * * __restrict pstrainPtr,                                                         \
      unsigned* __restrict isAdjustableVector,                                                     \
      std::size_t* __restrict yieldCounter,                                                        \
      const seissol::model::PlasticityData<Cfg>* __restrict plasticity,                            \
      Real<Cfg> oneMinusIntegratingFactor,                                                         \
      Real<Cfg> tV,                                                                                \
      Real<Cfg> timeStepWidth,                                                                     \
      size_t numElements,                                                                          \
      void* streamPtr);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::device::aux::plasticity
