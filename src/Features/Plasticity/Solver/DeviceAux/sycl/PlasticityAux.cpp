// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Features/Plasticity/Solver/DeviceAux/PlasticityAux.h"

#include "Config.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Common/Real.h"
#include "Features/Plasticity/Physics/DruckerPrager.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <cstddef>
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::plasticity {

namespace {
template <typename Cfg>
auto getrange(std::size_t size, std::size_t numElements) {
  if constexpr (Cfg::NumSimulations > 1) {
    return sycl::nd_range<1>({numElements * Cfg::NumSimulations * size},
                             {Cfg::NumSimulations * size});
  } else {
    return sycl::nd_range<1>({numElements * size}, {size});
  }
}

template <typename Cfg>
void plasticityNonlinearOfSolid(Real<Cfg>** __restrict nodalStressTensors,
                                Real<Cfg>** __restrict pstrainPtr,
                                unsigned* __restrict isAdjustableVector,
                                std::size_t* __restrict yieldCounter,
                                const seissol::model::PlasticityData<Cfg>* __restrict plasticity,
                                const unsigned* __restrict cellIndices,
                                Real<Cfg> oneMinusIntegratingFactor,
                                const size_t numElements,
                                void* streamPtr) {

  using real = Real<Cfg>;

  constexpr unsigned NumNodes = init::QStressNodal<Cfg>::Stop[multisim::BasisDim<Cfg>] -
                                init::QStressNodal<Cfg>::Start[multisim::BasisDim<Cfg>];

  auto& queue = *reinterpret_cast<sycl::queue*>(streamPtr);
  auto rng = getrange<Cfg>(NumNodes, numElements);

  queue.submit([&](sycl::handler& cgh) {
    const sycl::local_accessor<int> isAdjusted(1, cgh);

    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      const auto wid = item.get_group().get_group_id(0);
      const auto tid = item.get_local_id(0);

      real* qStressNodal = nodalStressTensors[wid];
      const seissol::model::PlasticityData<Cfg>& cellPlasticity = plasticity[cellIndices[wid]];

      // the yield criterion at the node of this work item (see model::plasticity::yieldOfNode)
      constexpr auto ElementTensorsColumn = multisim::linearDim<Cfg, init::QStressNodal<Cfg>>();
      real localStresses[model::plasticity::StressComponents];
#pragma unroll
      for (int i = 0; i < model::plasticity::StressComponents; ++i) {
        localStresses[i] = qStressNodal[tid + ElementTensorsColumn * i];
      }
      const auto node =
          model::plasticity::yieldOfNode(localStresses,
                                         cellPlasticity.cohesionTimesCosAngularFriction[tid],
                                         cellPlasticity.sinAngularFriction[tid],
                                         oneMinusIntegratingFactor);

      // the cell is adjusted where one of its nodes yields
      if (tid == 0) {
        isAdjusted[0] = static_cast<unsigned>(false);
      }
      item.barrier();
      if (node.yields) {
        isAdjusted[0] = static_cast<unsigned>(true);
      }
      item.barrier();

      // the relaxation of the node (see model::plasticity::relaxNode); the nodal stresses hold
      // the update of the stress afterwards
      if (isAdjusted[0]) {
        real* __restrict eta = pstrainPtr[wid] + tensor::QStressNodal<Cfg>::size();
        real* __restrict localPstrain = pstrainPtr[wid];

        real stressUpdate[model::plasticity::StressComponents];
        real pstrainUpdate[model::plasticity::StressComponents];
        const real etaGrowth = model::plasticity::relaxNode(
            node, cellPlasticity.mufactor, stressUpdate, pstrainUpdate);

#pragma unroll
        for (int i = 0; i < model::plasticity::StressComponents; ++i) {
          const int q = tid + ElementTensorsColumn * i;
          localPstrain[q] += pstrainUpdate[i];
          qStressNodal[q] = stressUpdate[i];
        }

        eta[tid] += etaGrowth;
      }
      if (tid == 0) {
        isAdjustableVector[wid] = isAdjusted[0];
        if (isAdjusted[0]) {
          const auto yieldCounterAR =
              sycl::atomic_ref<std::size_t,
                               sycl::memory_order::relaxed,
                               sycl::memory_scope::device,
                               sycl::access::address_space::global_space>(*yieldCounter);

          // count the cell for the FLOPs of the yield, once, as the CPU does
          yieldCounterAR.fetch_add(std::size_t{1}, sycl::memory_order::relaxed);
        }
      }
    });
  });
}
} // namespace

template <typename Cfg>
void plasticityNonlinear(Real<Cfg>** __restrict nodalStressTensors,
                         Real<Cfg>** __restrict pstrainPtr,
                         unsigned* __restrict isAdjustableVector,
                         std::size_t* __restrict yieldCounter,
                         const seissol::model::PlasticityData<Cfg>* __restrict plasticity,
                         const unsigned* __restrict cellIndices,
                         Real<Cfg> oneMinusIntegratingFactor,
                         const size_t numElements,
                         void* streamPtr) {
  // a material without the stresses of a solid has nothing that yields; the kernel, which reads the
  // six stresses, is not even compiled for it
  if constexpr (model::MaterialOf<Cfg>::SupportsPlasticity) {
    plasticityNonlinearOfSolid<Cfg>(nodalStressTensors,
                                    pstrainPtr,
                                    isAdjustableVector,
                                    yieldCounter,
                                    plasticity,
                                    cellIndices,
                                    oneMinusIntegratingFactor,
                                    numElements,
                                    streamPtr);
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
