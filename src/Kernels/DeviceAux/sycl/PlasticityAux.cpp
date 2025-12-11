// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/PlasticityAux.h"

#include "GeneratedCode/init.h"
#include "Solver/MultipleSimulations.h"

#include <cmath>
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::plasticity {

template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type squareRoot(T x) {
  return std::is_same<T, double>::value ? sqrt(x) : sqrtf(x);
}

template <typename Tensor>
constexpr size_t leadDim() {
  if constexpr (multisim::MultisimEnabled) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

auto getrange(std::size_t size, std::size_t numElements) {
  if constexpr (multisim::MultisimEnabled) {
    return sycl::nd_range<1>({numElements * multisim::NumSimulations * size},
                             {multisim::NumSimulations * size});
  } else {
    return sycl::nd_range<1>({numElements * size}, {size});
  }
}

void adjustDeviatoricTensors(real** __restrict nodalStressTensors,
                             real** __restrict prevNodal,
                             real** __restrict pstrainPtr,
                             unsigned* __restrict isAdjustableVector,
                             std::size_t* __restrict yieldCounter,
                             const seissol::model::PlasticityData* __restrict plasticity,
                             double oneMinusIntegratingFactor,
                             double tV,
                             double timeStepWidth,
                             const size_t numElements,
                             void* queuePtr) {
  constexpr unsigned NumNodes = tensor::QStressNodal::Shape[multisim::BasisFunctionDimension];
  auto queue = reinterpret_cast<sycl::queue*>(queuePtr);
  auto rng = getrange(NumNodes, numElements);

  queue->submit([&](sycl::handler& cgh) {
    sycl::local_accessor<int> isAdjusted(1, cgh);

    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      auto wid = item.get_group().get_group_id(0);
      auto tid = item.get_local_id(0);

      real* qStressNodal = nodalStressTensors[wid];
      real localStresses[NumStressComponents];

      constexpr auto ElementTensorsColumn = leadDim<init::QStressNodal>();
#pragma unroll
      for (int i = 0; i < NumStressComponents; ++i) {
        localStresses[i] = qStressNodal[tid + ElementTensorsColumn * i];
      }

      // 1. Compute the mean stress for each node
      real meanStress = (localStresses[0] + localStresses[1] + localStresses[2]) / 3.0f;

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
          plasticity[wid].cohesionTimesCosAngularFriction[tid % seissol::multisim::NumSimulations];
      const real sinAngularFriction =
          plasticity[wid].sinAngularFriction[tid % seissol::multisim::NumSimulations];
      real taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
      taulim = std::max(static_cast<real>(0.0), taulim);

      if (tid == 0) {
        isAdjusted[0] = static_cast<unsigned>(false);
      }
      item.barrier();

      // 5. Compute the yield factor
      real yieldfactor = 0.0;
      if (tau > taulim) {
        isAdjusted[0] = static_cast<unsigned>(true);
        yieldfactor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
      }

      // 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
      item.barrier();
      if (isAdjusted[0]) {
        const real factor = plasticity[wid].mufactor / (tV * oneMinusIntegratingFactor);

        const real* __restrict localPrevNodal = prevNodal[wid];
        real* __restrict eta = pstrainPtr[wid] + tensor::QStressNodal::size();
        real* __restrict localPstrain = pstrainPtr[wid];

        real dudtUpdate = 0;

#pragma unroll
        for (int i = 0; i < NumStressComponents; ++i) {
          const int q = tid + ElementTensorsColumn * i;

          const auto updatedStressNodal = localStresses[i] * yieldfactor;

          const real nodeDuDtPstrain = factor * (localPrevNodal[q] - updatedStressNodal);

          localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
          qStressNodal[q] = updatedStressNodal;

          dudtUpdate += nodeDuDtPstrain * nodeDuDtPstrain;
        }

        eta[tid] += timeStepWidth * std::sqrt(static_cast<real>(0.5) * dudtUpdate);

        auto yieldCounterAR =
            sycl::atomic_ref<std::size_t,
                             sycl::memory_order::relaxed,
                             sycl::memory_scope::device,
                             sycl::access::address_space::global_space>(yieldCounter);

        // update the FLOPs that we've been here
        yieldCounterAR.fetch_add(1, sycl::memory_order::relaxed);
      }
      if (tid == 0) {
        isAdjustableVector[wid] = isAdjusted[0];
      }
    });
  });
}

void computePstrains(real** pstrains,
                     const seissol::model::PlasticityData* plasticityData,
                     real** dofs,
                     real** prevDofs,
                     real** dUdTpstrain,
                     double T_v,
                     double oneMinusIntegratingFactor,
                     double timeStepWidth,
                     unsigned* isAdjustableVector,
                     size_t numElements,
                     void* queuePtr) {
  constexpr unsigned numNodes = tensor::Q::Shape[multisim::BasisFunctionDimension];
  auto queue = reinterpret_cast<sycl::queue*>(queuePtr);
  auto rng = getrange(numNodes, numElements);

  queue->submit([&](sycl::handler& cgh) {
    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      auto wid = item.get_group().get_group_id(0);
      auto lid = item.get_local_id(0);

      if (isAdjustableVector[wid]) {
        real* localDofs = dofs[wid];
        real* localPrevDofs = prevDofs[wid];
        const seissol::model::PlasticityData* localData = &plasticityData[wid];
        real* localPstrain = pstrains[wid];
        real* localDuDtPstrain = dUdTpstrain[wid];

#pragma unroll
        for (int i = 0; i < NumStressComponents; ++i) {
          int q = lid + i * leadDim<init::Q>();
          real factor = localData->mufactor / (T_v * oneMinusIntegratingFactor);
          real nodeDuDtPstrain = factor * (localPrevDofs[q] - localDofs[q]);

          static_assert(leadDim<init::QStress>() == leadDim<init::Q>());
          localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
          localDuDtPstrain[q] = nodeDuDtPstrain;
        }
      }
    });
  });
}

void updateQEtaNodal(real** QEtaNodalPtrs,
                     real** QStressNodalPtrs,
                     double timeStepWidth,
                     unsigned* isAdjustableVector,
                     size_t numElements,
                     void* queuePtr) {
  auto queue = reinterpret_cast<sycl::queue*>(queuePtr);
  auto rng = getrange(tensor::QStressNodal::Shape[multisim::BasisFunctionDimension], numElements);

  queue->submit([&](sycl::handler& cgh) {
    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      auto wid = item.get_group().get_group_id(0);
      auto lid = item.get_local_id(0);

      if (isAdjustableVector[wid]) {
        real* __restrict localQEtaNodal = QEtaNodalPtrs[wid] + tensor::QStress::size();
        real* __restrict localQStressNodal = QStressNodalPtrs[wid];
        real factor{0.0};

        constexpr auto ld = leadDim<init::QStressNodal>();
#pragma unroll
        for (int i = 0; i < NumStressComponents; ++i) {
          factor += localQStressNodal[lid + i * ld] * localQStressNodal[lid + i * ld];
        }

        localQEtaNodal[lid] += timeStepWidth * std::sqrt(static_cast<real>(0.5) * factor);
      }
    });
  });
}

} // namespace seissol::kernels::device::aux::plasticity
