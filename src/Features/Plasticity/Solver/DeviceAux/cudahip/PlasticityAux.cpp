// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Features/Plasticity/Solver/DeviceAux/PlasticityAux.h"

#include "Config.h"
#include "Core/Common/MultipleSimulations.h"
#include "Features/Plasticity/Physics/DruckerPrager.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Physics/Datastructures.h"

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
                               const unsigned* __restrict cellIndices,
                               Real<Cfg> oneMinusIntegratingFactor) {

  using real = Real<Cfg>;

  real* __restrict qStressNodal = nodalStressTensors[blockIdx.x];
  const seissol::model::PlasticityData<Cfg>& cellPlasticity = plasticity[cellIndices[blockIdx.x]];

  // the yield criterion at the node of this thread (see model::plasticity::yieldOfNode)
  constexpr auto ElementTensorsColumn = multisim::linearDim<Cfg, init::QStressNodal<Cfg>>();
  real localStresses[model::plasticity::StressComponents];
#pragma unroll
  for (int i = 0; i < model::plasticity::StressComponents; ++i) {
    localStresses[i] = qStressNodal[linearidx<Cfg>() + ElementTensorsColumn * i];
  }
  const auto node = model::plasticity::yieldOfNode(
      localStresses,
      cellPlasticity.cohesionTimesCosAngularFriction[linearidx<Cfg>()],
      cellPlasticity.sinAngularFriction[linearidx<Cfg>()],
      oneMinusIntegratingFactor);

  // the cell is adjusted where one of its nodes yields
  __shared__ bool isAdjusted;
  if (linearidx<Cfg>() == 0) {
    isAdjusted = false;
  }
  __syncthreads();
  if (node.yields) {
    isAdjusted = true;
  }
  __syncthreads();

  // the relaxation of the node (see model::plasticity::relaxNode); the nodal stresses hold the
  // update of the stress afterwards
  if (isAdjusted) {
    real* __restrict eta = pstrainPtr[blockIdx.x] + tensor::QStressNodal<Cfg>::size();
    real* __restrict localPstrain = pstrainPtr[blockIdx.x];

    real stressUpdate[model::plasticity::StressComponents];
    real pstrainUpdate[model::plasticity::StressComponents];
    const real etaGrowth =
        model::plasticity::relaxNode(node, cellPlasticity.mufactor, stressUpdate, pstrainUpdate);

#pragma unroll
    for (int i = 0; i < model::plasticity::StressComponents; ++i) {
      const int q = linearidx<Cfg>() + ElementTensorsColumn * i;
      localPstrain[q] += pstrainUpdate[i];
      qStressNodal[q] = stressUpdate[i];
    }

    eta[linearidx<Cfg>()] += etaGrowth;
  }
  if (linearidx<Cfg>() == 0) {
    isAdjustableVector[blockIdx.x] = isAdjusted;
    if (isAdjusted) {
      // count the cell for the FLOPs of the yield, once, as the CPU does
      // (there's no atomicAdd for unsigned long / sometimes size_t, so take one of the other ones)
      static_assert(sizeof(unsigned long long) == sizeof(std::size_t));
      atomicAdd(reinterpret_cast<unsigned long long*>(yieldCounter), 1);
    }
  }
}

template <typename Cfg>
void plasticityNonlinear(Real<Cfg>** __restrict nodalStressTensors,
                         Real<Cfg>** __restrict pstrainPtr,
                         unsigned* __restrict isAdjustableVector,
                         std::size_t* __restrict yieldCounter,
                         const seissol::model::PlasticityData<Cfg>* __restrict plasticity,
                         const unsigned* __restrict cellIndices,
                         Real<Cfg> oneMinusIntegratingFactor,
                         size_t numElements,
                         void* streamPtr) {
  // a material without the stresses of a solid has nothing that yields; the kernel, which reads the
  // six stresses, is not even compiled for it
  if constexpr (model::MaterialOf<Cfg>::SupportsPlasticity) {
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
                                                           cellIndices,
                                                           oneMinusIntegratingFactor);
  }
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template void plasticityNonlinear<Cfg>(                                                          \
      Real<Cfg> * * __restrict nodalStressTensors,                                                 \
      Real<Cfg> * * __restrict pstrainPtr,                                                         \
      unsigned* __restrict isAdjustableVector,                                                     \
      std::size_t* __restrict yieldCounter,                                                        \
      const seissol::model::PlasticityData<Cfg>* __restrict plasticity,                            \
      const unsigned* __restrict cellIndices,                                                      \
      Real<Cfg> oneMinusIntegratingFactor,                                                         \
      size_t numElements,                                                                          \
      void* streamPtr);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::device::aux::plasticity
