// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Stephanie Wollherr

#include "Plasticity.h"

#include "Alignment.h"
#include "Common/Marker.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Model/Plasticity.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>

#ifdef ACL_DEVICE
#include "DeviceAux/PlasticityAux.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
using namespace device;
#endif

#ifndef ACL_DEVICE
#include <utils/logger.h>
#endif

#ifndef NDEBUG
#include <cstdint>
#endif

namespace seissol::kernels {
template <typename Cfg>
std::size_t
    Plasticity<Cfg>::computePlasticity(real oneMinusIntegratingFactor,
                                       real timeStepWidth,
                                       real tV,
                                       const GlobalData<Cfg>* global,
                                       const seissol::model::PlasticityData<Cfg>* plasticityData,
                                       real degreesOfFreedom[tensor::Q<Cfg>::size()],
                                       real* pstrain) {
  // a material without the stresses of a solid has nothing that yields
  if constexpr (!model::MaterialOf<Cfg>::SupportsPlasticity) {
    return 0;
  }

  assert(reinterpret_cast<uintptr_t>(degreesOfFreedom) % Vectorsize == 0);

  alignas(Alignment) real qStressNodal[tensor::QStressNodal<Cfg>::size()]{};

  alignas(Alignment) real meanStress[tensor::meanStress<Cfg>::size()]{};
  alignas(Alignment) real secondInvariant[tensor::secondInvariant<Cfg>::size()]{};
  alignas(Alignment) real tau[tensor::secondInvariant<Cfg>::size()]{};
  alignas(Alignment) real taulim[tensor::meanStress<Cfg>::size()]{};
  alignas(Alignment) real yieldFactor[tensor::yieldFactor<Cfg>::size()]{};

  static_assert(tensor::secondInvariant<Cfg>::size() == tensor::meanStress<Cfg>::size(),
                "Second invariant tensor and mean stress tensor must be of the same size().");
  static_assert(tensor::yieldFactor<Cfg>::size() <= tensor::meanStress<Cfg>::size(),
                "Yield factor tensor must be smaller than mean stress tensor.");

  /* Convert modal to nodal and add sigma0.
   * Stores s_{ij} := sigma_{ij} + sigma0_{ij} for every node.
   * sigma0 is constant
   * also stores the previous DOFs before adding the new initial loading
   */

  kernel::plConvertToNodal<Cfg> m2nKrnl;
  m2nKrnl.bindGlobals(*global);
  m2nKrnl.QStress = degreesOfFreedom;
  m2nKrnl.QStressNodal = qStressNodal;
  m2nKrnl.initialLoading = plasticityData->initialLoading;
  m2nKrnl.execute();

  // Computes m = s_{ii} / 3.0 for every node
  kernel::plComputeMean<Cfg> cmKrnl;
  cmKrnl.bindGlobals(*global);
  cmKrnl.meanStress = meanStress;
  cmKrnl.QStressNodal = qStressNodal;
  cmKrnl.execute();

  /* Compute s_{ij} := s_{ij} - m delta_{ij},
   * where delta_{ij} = 1 if i == j else 0.
   * Thus, s_{ij} contains the deviatoric stresses. */
  kernel::plSubtractMean<Cfg> smKrnl;
  smKrnl.bindGlobals(*global);
  smKrnl.meanStress = meanStress;
  smKrnl.QStressNodal = qStressNodal;
  smKrnl.execute();

  // Compute I_2 = 0.5 s_{ij} s_ji for every node
  kernel::plComputeSecondInvariant<Cfg> siKrnl;
  siKrnl.bindGlobals(*global);
  siKrnl.secondInvariant = secondInvariant;
  siKrnl.QStressNodal = qStressNodal;
  siKrnl.execute();

// tau := sqrt(I_2) for every node
#pragma omp simd
  for (std::size_t ip = 0; ip < tensor::secondInvariant<Cfg>::size(); ++ip) {
    tau[ip] = std::sqrt(secondInvariant[ip]);
  }

// Compute tau_c for every node
#pragma omp simd
  for (std::size_t ip = 0; ip < tensor::meanStress<Cfg>::size(); ++ip) {
    taulim[ip] = std::max(static_cast<real>(0.0),
                          plasticityData->cohesionTimesCosAngularFriction[ip] -
                              meanStress[ip] * plasticityData->sinAngularFriction[ip]);
  }

  int32_t adjust = 0;

#pragma omp simd reduction(max : adjust)
  for (std::size_t ip = 0; ip < tensor::yieldFactor<Cfg>::size(); ++ip) {
    // Compute yield := (t_c / tau - 1) r for every node,
    // where r = 1 - exp(-timeStepWidth / tV)
    const auto doesYield = tau[ip] > taulim[ip];
    // Combine with the previous value: under the reduction, each SIMD lane keeps only its own
    // last value, so a plain assignment would discard all but the last chunk of nodes.
    adjust = std::max(adjust, static_cast<int32_t>(doesYield ? 1 : 0));
    const auto ifYield =
        (taulim[ip] / tau[ip] - static_cast<real>(1.0)) * oneMinusIntegratingFactor;
    yieldFactor[ip] = doesYield ? ifYield : 0;
  }

  if (adjust != 0) {
    const real factor = plasticityData->mufactor / (tV * oneMinusIntegratingFactor);

    // calculate plastic strain
    constexpr std::size_t NumNodes = init::QStressNodal<Cfg>::Stop[multisim::BasisDim<Cfg>] -
                                     init::QStressNodal<Cfg>::Start[multisim::BasisDim<Cfg>];

    real* __restrict qEtaNodal = &pstrain[tensor::QStressNodal<Cfg>::size()];

    /**
     * Compute sigma_{ij} := sigma_{ij} + yield s_{ij} for every node
     * and store as modal basis.
     *
     * Remark: According to Wollherr et al., the update formula (13) should be
     *
     * sigmaNew_{ij} := f^* s_{ij} + m delta_{ij} - sigma0_{ij}
     *
     * where f^* = r tau_c / tau + (1 - r) = 1 + yield. Adding 0 to (13) gives
     *
     * sigmaNew_{ij} := f^* s_{ij} + m delta_{ij} - sigma0_{ij}
     *                  + sigma_{ij} + sigma0_{ij} - sigma_{ij} - sigma0_{ij}
     *                = f^* s_{ij} + sigma_{ij} - s_{ij}
     *                = sigma_{ij} + (f^* - 1) s_{ij}
     *                = sigma_{ij} + yield s_{ij}
     */

    constexpr auto NumTotalPoints = NumNodes * Cfg::NumSimulations;

#pragma omp simd
    for (std::size_t qp = 0; qp < NumTotalPoints; ++qp) {

      real dudtPstrainSqAcc = 0;

#pragma unroll
      for (std::size_t x = 0; x < 6; ++x) {
        const auto q = qp + NumTotalPoints * x;
        /**
         * Equation (10) from Wollherr et al.:
         *
         * d/dt strain_{ij} = (sigma_{ij} + sigma0_{ij} - P_{ij}(sigma)) / (2mu tV)
         *
         * where (11)
         *
         * P_{ij}(sigma) = { tau_c/tau s_{ij} + m delta_{ij}         if     tau >= taulim
         *                 { sigma_{ij} + sigma0_{ij}                else
         *
         * Thus,
         *
         * d/dt strain_{ij} = { (1 - tau_c/tau) / (2mu tV) s_{ij}   if     tau >= taulim
         *                    { 0                                    else
         *
         * Consider tau >= taulim first. We have (1 - tau_c/tau) = -yield / r. Therefore,
         *
         * d/dt strain_{ij} = -1 / (2mu tV r) yield s_{ij}
         *                  = -1 / (2mu tV r) (sigmaNew_{ij} - sigma_{ij})
         *                  = (sigma_{ij} - sigmaNew_{ij}) / (2mu tV r)
         *                  = -yield s_{ij} / (2mu tV r)
         *
         * If tau < taulim, then sigma_{ij} - sigmaNew_{ij} = 0.
         */
        const auto qStressNodalUpdate = qStressNodal[q] * yieldFactor[qp];
        const auto dudtPstrain = -factor * qStressNodalUpdate;

        // Integrate with explicit Euler
        pstrain[q] += timeStepWidth * dudtPstrain;

        // now contains the update for qStressNodal (cf. below)
        qStressNodal[q] = qStressNodalUpdate;

        dudtPstrainSqAcc += dudtPstrain * dudtPstrain;
      }

      // eta := int_0^t sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt) dt
      // Approximate with eta += timeStepWidth * sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt)

      qEtaNodal[qp] += timeStepWidth * std::sqrt(static_cast<real>(0.5) * dudtPstrainSqAcc);
    }

    kernel::plConvertToModal<Cfg> adjKrnl;
    adjKrnl.QStress = degreesOfFreedom;
    adjKrnl.bindGlobals(*global);
    adjKrnl.QStressNodal = qStressNodal;
    adjKrnl.execute();

    return 1;
  }
  return 0;
}

template <typename Cfg>
void Plasticity<Cfg>::computePlasticityBatched(
    SEISSOL_GPU_PARAM real timeStepWidth,
    SEISSOL_GPU_PARAM real tV,
    SEISSOL_GPU_PARAM const GlobalData<Cfg>* global,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::model::PlasticityData<Cfg>* plasticityData,
    SEISSOL_GPU_PARAM std::size_t* yieldCounter,
    SEISSOL_GPU_PARAM unsigned* isAdjustableVector,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // a material without the stresses of a solid has nothing that yields
  if constexpr (!model::MaterialOf<Cfg>::SupportsPlasticity) {
    return;
  }

  using namespace seissol::recording;

  static_assert(tensor::Q<Cfg>::Shape[0] == tensor::QStressNodal<Cfg>::Shape[0],
                "modal and nodal dofs must have the same leading dimensions");
  static_assert(tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg>] == tensor::v<Cfg>::Shape[0],
                "modal dofs and vandermonde matrix must have the same leading dimensions");

  const ConditionalKey key(*KernelNames::Plasticity);
  auto* defaultStream = runtime.stream();

  if (table.find(key) != table.end()) {
    const auto oneMinusIntegratingFactor = computeRelaxTime(tV, timeStepWidth);

    auto& entry = table[key];
    const size_t numElements = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize();

    // Convert modal to nodal
    real** modalStressTensors = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    real** nodalStressTensors =
        (entry.get<real*>(inner_keys::Wp::Id::NodalStressTensor))->getDeviceDataPtr();

    static_assert(kernel::gpu_plConvertToNodal<Cfg>::TmpMaxMemRequiredInBytes == 0);
    real** initLoad = (entry.get<real*>(inner_keys::Wp::Id::InitialLoad))->getDeviceDataPtr();
    kernel::gpu_plConvertToNodal<Cfg> m2nKrnl;
    m2nKrnl.bindGlobals(*global);
    m2nKrnl.QStress = const_cast<const real**>(modalStressTensors);
    m2nKrnl.QStressNodal = nodalStressTensors;
    m2nKrnl.initialLoading = const_cast<const real**>(initLoad);
    m2nKrnl.streamPtr = defaultStream;
    m2nKrnl.numElements = numElements;
    m2nKrnl.execute();

    real** pstrains = entry.get<real*>(inner_keys::Wp::Id::Pstrains)->getDeviceDataPtr();

    device::aux::plasticity::plasticityNonlinear(nodalStressTensors,
                                                 pstrains,
                                                 isAdjustableVector,
                                                 yieldCounter,
                                                 plasticityData,
                                                 oneMinusIntegratingFactor,
                                                 tV,
                                                 timeStepWidth,
                                                 numElements,
                                                 defaultStream);

    kernel::gpu_plConvertToModal<Cfg> n2mKrnl;
    n2mKrnl.bindGlobals(*global);
    n2mKrnl.QStressNodal = const_cast<const real**>(nodalStressTensors);
    n2mKrnl.QStress = modalStressTensors;
    n2mKrnl.streamPtr = defaultStream;
    n2mKrnl.flags = isAdjustableVector;
    n2mKrnl.numElements = numElements;
    n2mKrnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif // ACL_DEVICE
}

template <typename Cfg>
std::pair<PerformanceEstimate, PerformanceEstimate> Plasticity<Cfg>::metrics() {
  // reset flops
  PerformanceEstimate check;
  PerformanceEstimate yield;

  // flops from checking, i.e. outside if (adjust) {}
  check += PerformanceEstimate::fromKernel<kernel::plConvertToNodal<Cfg>>();

  // compute mean stress
  check += PerformanceEstimate::fromKernel<kernel::plComputeMean<Cfg>>();

  // subtract mean stress
  check += PerformanceEstimate::fromKernel<kernel::plSubtractMean<Cfg>>();

  // compute second invariant
  check += PerformanceEstimate::fromKernel<kernel::plComputeSecondInvariant<Cfg>>();

  // compute taulim (1 add, 1 mul, max NOT counted)
  check.nonzeroFlop += static_cast<std::uint64_t>(2 * tensor::meanStress<Cfg>::size());
  check.hardwareFlop += static_cast<std::uint64_t>(2 * tensor::meanStress<Cfg>::size());

  // check for yield (NOT counted, as it would require counting the number of yielding points)

  // flops from plastic yielding, i.e. inside if (adjust) {}
  yield += PerformanceEstimate::fromKernel<kernel::plConvertToModal<Cfg>>();

  // manually counted
  yield.nonzeroFlop += static_cast<std::uint64_t>(tensor::QStressNodal<Cfg>::size() * 6);
  yield.hardwareFlop += static_cast<std::uint64_t>(tensor::QStressNodal<Cfg>::size() * 6);
  yield.nonzeroFlop += static_cast<std::uint64_t>(tensor::QEtaNodal<Cfg>::size() * 3);
  yield.hardwareFlop += static_cast<std::uint64_t>(tensor::QEtaNodal<Cfg>::size() * 3);

  return {check, yield};
}
#define SEISSOL_INSTANTIATE(Cfg) template class Plasticity<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
