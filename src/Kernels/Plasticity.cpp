// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Stephanie Wollherr

#include "Plasticity.h"

#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include <Alignment.h>
#include <DataTypes/ConditionalTable.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Model/Plasticity.h>
#include <Parallel/Runtime/Stream.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>

#include "utils/logger.h"

#ifdef ACL_DEVICE
#include "DeviceAux/PlasticityAux.h"
#include "device.h"
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#include <Solver/MultipleSimulations.h>
using namespace device;
#endif

#ifndef NDEBUG
#include <cstdint>
#endif

namespace seissol::kernels {
template<typename Cfg>
std::size_t Plasticity<Cfg>::computePlasticity(double oneMinusIntegratingFactor,
                                          double timeStepWidth,
                                          double tV,
                                          const GlobalData* global,
                                          const seissol::model::PlasticityData* plasticityData,
                                          real degreesOfFreedom[tensor::Q<Cfg>::size()],
                                          real* pstrain) {
  if constexpr (multisim::MultisimEnabled) {
    // TODO: really still the case?
    logError() << "Plasticity does not work with multiple simulations";
  }
  assert(reinterpret_cast<uintptr_t>(degreesOfFreedom) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(global->vandermondeMatrix) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(global->vandermondeMatrixInverse) % Alignment == 0);

  alignas(Alignment) real qStressNodal[tensor::QStressNodal<Cfg>::size()];
  alignas(Alignment) real qEtaNodal[tensor::QEtaNodal<Cfg>::size()];
  alignas(Alignment) real qEtaModal[tensor::QEtaModal<Cfg>::size()];
  alignas(Alignment) real meanStress[tensor::meanStress<Cfg>::size()];
  alignas(Alignment) real secondInvariant[tensor::secondInvariant<Cfg>::size()];
  alignas(Alignment) real tau[tensor::secondInvariant<Cfg>::size()];
  alignas(Alignment) real taulim[tensor::meanStress<Cfg>::size()];
  alignas(Alignment) real yieldFactor[tensor::yieldFactor<Cfg>::size()];
  alignas(Alignment) real dudtPstrain[tensor::QStress<Cfg>::size()];

  static_assert(tensor::secondInvariant<Cfg>::size() == tensor::meanStress<Cfg>::size(),
                "Second invariant tensor and mean stress tensor must be of the same size().");
  static_assert(tensor::yieldFactor<Cfg>::size() <= tensor::meanStress<Cfg>::size(),
                "Yield factor tensor must be smaller than mean stress tensor.");

  // copy dofs for later comparison, only first dof of stresses required
  //  @todo multiple sims

  real prevDegreesOfFreedom[tensor::QStress<Cfg>::size()];
  for (std::size_t q = 0; q < tensor::QStress<Cfg>::size(); ++q) {
    prevDegreesOfFreedom[q] = degreesOfFreedom[q];
  }

  /* Convert modal to nodal and add sigma0.
   * Stores s_{ij} := sigma_{ij} + sigma0_{ij} for every node.
   * sigma0 is constant */
  kernel::plConvertToNodal<Cfg> m2nKrnl;
  m2nKrnl.v = global->vandermondeMatrix;
  m2nKrnl.QStress = degreesOfFreedom;
  m2nKrnl.QStressNodal = qStressNodal;
  m2nKrnl.replicateInitialLoading = init::replicateInitialLoading<Cfg>::Values;
  m2nKrnl.initialLoading = plasticityData->initialLoading;
  m2nKrnl.execute();

  // Computes m = s_{ii} / 3.0 for every node
  kernel::plComputeMean<Cfg> cmKrnl;
  cmKrnl.meanStress = meanStress;
  cmKrnl.QStressNodal = qStressNodal;
  cmKrnl.selectBulkAverage = init::selectBulkAverage<Cfg>::Values;
  cmKrnl.execute();

  /* Compute s_{ij} := s_{ij} - m delta_{ij},
   * where delta_{ij} = 1 if i == j else 0.
   * Thus, s_{ij} contains the deviatoric stresses. */
  kernel::plSubtractMean<Cfg> smKrnl;
  smKrnl.meanStress = meanStress;
  smKrnl.QStressNodal = qStressNodal;
  smKrnl.selectBulkNegative = init::selectBulkNegative<Cfg>::Values;
  smKrnl.execute();

  // Compute I_2 = 0.5 s_{ij} s_ji for every node
  kernel::plComputeSecondInvariant<Cfg> siKrnl;
  siKrnl.secondInvariant = secondInvariant;
  siKrnl.QStressNodal = qStressNodal;
  siKrnl.weightSecondInvariant = init::weightSecondInvariant<Cfg>::Values;
  siKrnl.execute();

  // tau := sqrt(I_2) for every node
  for (std::size_t ip = 0; ip < tensor::secondInvariant<Cfg>::size(); ++ip) {
    tau[ip] = sqrt(secondInvariant[ip]);
  }

  // Compute tau_c for every node
  for (std::size_t ip = 0; ip < tensor::meanStress<Cfg>::size(); ++ip) {
    taulim[ip] = std::max(static_cast<real>(0.0),
                          plasticityData->cohesionTimesCosAngularFriction -
                              meanStress[ip] * plasticityData->sinAngularFriction);
  }

  bool adjust = false;
  for (std::size_t ip = 0; ip < tensor::yieldFactor<Cfg>::size(); ++ip) {
    // Compute yield := (t_c / tau - 1) r for every node,
    // where r = 1 - exp(-timeStepWidth / tV)
    if (tau[ip] > taulim[ip]) {
      adjust = true;
      yieldFactor[ip] = (taulim[ip] / tau[ip] - 1.0) * oneMinusIntegratingFactor;
    } else {
      yieldFactor[ip] = 0.0;
    }
  }

  if (adjust) {
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
    kernel::plAdjustStresses<Cfg> adjKrnl;
    adjKrnl.QStress = degreesOfFreedom;
    adjKrnl.vInv = global->vandermondeMatrixInverse;
    adjKrnl.QStressNodal = qStressNodal;
    adjKrnl.yieldFactor = yieldFactor;
    adjKrnl.execute();

    // calculate plastic strain
    for (unsigned q = 0; q < tensor::QStress<Cfg>::size(); ++q) {
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
       *
       * If tau < taulim, then sigma_{ij} - sigmaNew_{ij} = 0.
       */
      const real factor = plasticityData->mufactor / (tV * oneMinusIntegratingFactor);
      dudtPstrain[q] = factor * (prevDegreesOfFreedom[q] - degreesOfFreedom[q]);
      // Integrate with explicit Euler
      pstrain[q] += timeStepWidth * dudtPstrain[q];
    }
    /* Convert modal to nodal */
    kernel::plConvertToNodalNoLoading<Cfg> m2nKrnlDudtPstrain;
    m2nKrnlDudtPstrain.v = global->vandermondeMatrix;
    m2nKrnlDudtPstrain.QStress = dudtPstrain;
    m2nKrnlDudtPstrain.QStressNodal = qStressNodal;
    m2nKrnlDudtPstrain.execute();

    // Sizes:
    for (unsigned q = 0; q < tensor::QEtaModal<Cfg>::size(); ++q) {
      qEtaModal[q] = pstrain[tensor::QStress<Cfg>::size() + q];
    }

    /* Convert modal to nodal */
    kernel::plConvertEtaModal2Nodal<Cfg> m2nEtaKrnl;
    m2nEtaKrnl.v = global->vandermondeMatrix;
    m2nEtaKrnl.QEtaModal = qEtaModal;
    m2nEtaKrnl.QEtaNodal = qEtaNodal;
    m2nEtaKrnl.execute();

    auto qStressNodalView = init::QStressNodal<Cfg>::view::create(qStressNodal);
    const auto numNodes = qStressNodalView.shape(multisim::BasisFunctionDimension);
    for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
      for (std::size_t i = 0; i < numNodes; ++i) {
        // eta := int_0^t sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt) dt
        // Approximate with eta += timeStepWidth * sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt)
        qEtaNodal[i * multisim::NumSimulations + s] =
            std::max(static_cast<real>(0.0), qEtaNodal[i * multisim::NumSimulations + s]) +
            timeStepWidth * sqrt(0.5 * (multisim::multisimWrap(qStressNodalView, s, i, 0) *
                                            multisim::multisimWrap(qStressNodalView, s, i, 0) +
                                        multisim::multisimWrap(qStressNodalView, s, i, 1) *
                                            multisim::multisimWrap(qStressNodalView, s, i, 1) +
                                        multisim::multisimWrap(qStressNodalView, s, i, 2) *
                                            multisim::multisimWrap(qStressNodalView, s, i, 2) +
                                        multisim::multisimWrap(qStressNodalView, s, i, 3) *
                                            multisim::multisimWrap(qStressNodalView, s, i, 3) +
                                        multisim::multisimWrap(qStressNodalView, s, i, 4) *
                                            multisim::multisimWrap(qStressNodalView, s, i, 4) +
                                        multisim::multisimWrap(qStressNodalView, s, i, 5) *
                                            multisim::multisimWrap(qStressNodalView, s, i, 5)));
      }
    }

    /* Convert nodal to modal */
    kernel::plConvertEtaNodal2Modal<Cfg> n2mEtaKrnl;
    n2mEtaKrnl.vInv = global->vandermondeMatrixInverse;
    n2mEtaKrnl.QEtaNodal = qEtaNodal;
    n2mEtaKrnl.QEtaModal = qEtaModal;
    n2mEtaKrnl.execute();
    for (std::size_t q = 0; q < tensor::QEtaModal<Cfg>::size(); ++q) {
      pstrain[tensor::QStress<Cfg>::size() + q] = qEtaModal[q];
    }
    return 1;
  }
  return 0;
}

template<typename Cfg>
void Plasticity<Cfg>::computePlasticityBatched(
    double timeStepWidth,
    double tV,
    const GlobalData* global,
    initializer::recording::ConditionalPointersToRealsTable& table,
    seissol::model::PlasticityData* plasticityData,
    std::size_t* yieldCounter,
    unsigned* isAdjustableVector,
    seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  static_assert(tensor::Q<Cfg>::Shape[0] == tensor::QStressNodal<Cfg>::Shape[0],
                "modal and nodal dofs must have the same leading dimensions");
  static_assert(tensor::Q<Cfg>::Shape[multisim::BasisFunctionDimension] == tensor::v<Cfg>::Shape[0],
                "modal dofs and vandermonde matrix must hage the same leading dimensions");

  DeviceInstance& device = DeviceInstance::getInstance();
  ConditionalKey key(*KernelNames::Plasticity);
  auto defaultStream = runtime.stream();

  if (table.find(key) != table.end()) {
    const auto oneMinusIntegratingFactor = computeRelaxTime(tV, timeStepWidth);

    auto& entry = table[key];
    const size_t numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();

    // copy dofs for later comparison, only first dof of stresses required
    constexpr auto DofsSize = tensor::Q<Cfg>::Size;

    real** prevDofs = (entry.get(inner_keys::Wp::Id::PrevDofs))->getDeviceDataPtr();
    real** dofsPtrs = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    device.algorithms.streamBatchedData(
        const_cast<const real**>(dofsPtrs), prevDofs, DofsSize, numElements, defaultStream);

    // Convert modal to nodal
    real** modalStressTensors = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    real** nodalStressTensors =
        (entry.get(inner_keys::Wp::Id::NodalStressTensor))->getDeviceDataPtr();

    assert(global->replicateStresses != nullptr && "replicateStresses has not been initialized");
    static_assert(kernel::gpu_plConvertToNodal<Cfg>::TmpMaxMemRequiredInBytes == 0);
    real** initLoad = (entry.get(inner_keys::Wp::Id::InitialLoad))->getDeviceDataPtr();
    kernel::gpu_plConvertToNodal<Cfg> m2nKrnl;
    m2nKrnl.v = global->vandermondeMatrix;
    m2nKrnl.QStress = const_cast<const real**>(modalStressTensors);
    m2nKrnl.QStressNodal = nodalStressTensors;
    m2nKrnl.replicateInitialLoadingM = global->replicateStresses;
    m2nKrnl.initialLoadingM = const_cast<const real**>(initLoad);
    m2nKrnl.streamPtr = defaultStream;
    m2nKrnl.numElements = numElements;
    m2nKrnl.execute();

    device::aux::plasticity::adjustDeviatoricTensors(nodalStressTensors,
                                                     isAdjustableVector,
                                                     plasticityData,
                                                     oneMinusIntegratingFactor,
                                                     numElements,
                                                     defaultStream);

    // count how many elements needs to be adjusted
    device.algorithms.reduceVector(yieldCounter,
                                   isAdjustableVector,
                                   true,
                                   numElements,
                                   ::device::ReductionType::Add,
                                   defaultStream);

    // convert back to modal (taking into account the adjustment)
    static_assert(kernel::gpu_plConvertToModal<Cfg>::TmpMaxMemRequiredInBytes == 0);
    kernel::gpu_plConvertToModal<Cfg> n2mKrnl;
    n2mKrnl.vInv = global->vandermondeMatrixInverse;
    n2mKrnl.QStressNodal = const_cast<const real**>(nodalStressTensors);
    n2mKrnl.QStress = modalStressTensors;
    n2mKrnl.streamPtr = defaultStream;
    n2mKrnl.flags = isAdjustableVector;
    n2mKrnl.numElements = numElements;
    n2mKrnl.execute();

    // prepare memory
    real** qEtaNodalPtrs = (entry.get(inner_keys::Wp::Id::QEtaNodal))->getDeviceDataPtr();
    real** dUdTpstrainPtrs = (entry.get(inner_keys::Wp::Id::DuDtStrain))->getDeviceDataPtr();

    static_assert(tensor::QStress<Cfg>::Size == tensor::QStressNodal<Cfg>::Size);

    // ------------------------------------------------------------------------------
    real** pstrains = entry.get(inner_keys::Wp::Id::Pstrains)->getDeviceDataPtr();
    real** dofs = modalStressTensors;
    device::aux::plasticity::computePstrains(pstrains,
                                             plasticityData,
                                             dofs,
                                             prevDofs,
                                             dUdTpstrainPtrs,
                                             tV,
                                             oneMinusIntegratingFactor,
                                             timeStepWidth,
                                             isAdjustableVector,
                                             numElements,
                                             defaultStream);

    // Convert modal to nodal
    static_assert(kernel::gpu_plConvertToNodalNoLoading<Cfg>::TmpMaxMemRequiredInBytes == 0);
    kernel::gpu_plConvertToNodalNoLoading<Cfg> m2nKrnlDudtPstrain;
    m2nKrnlDudtPstrain.v = global->vandermondeMatrix;
    m2nKrnlDudtPstrain.QStress = const_cast<const real**>(dUdTpstrainPtrs);
    m2nKrnlDudtPstrain.QStressNodal = nodalStressTensors;
    m2nKrnlDudtPstrain.streamPtr = defaultStream;
    m2nKrnlDudtPstrain.flags = isAdjustableVector;
    m2nKrnlDudtPstrain.numElements = numElements;
    m2nKrnlDudtPstrain.execute();

    // Convert modal to nodal
    static_assert(kernel::gpu_plConvertEtaModal2Nodal<Cfg>::TmpMaxMemRequiredInBytes == 0);
    kernel::gpu_plConvertEtaModal2Nodal<Cfg> m2nEtaKrnl;
    m2nEtaKrnl.v = global->vandermondeMatrix;
    m2nEtaKrnl.QEtaModal = const_cast<const real**>(pstrains);
    m2nEtaKrnl.extraOffset_QEtaModal = tensor::QStress<Cfg>::size();
    m2nEtaKrnl.QEtaNodal = qEtaNodalPtrs;
    m2nEtaKrnl.streamPtr = defaultStream;
    m2nEtaKrnl.flags = isAdjustableVector;
    m2nEtaKrnl.numElements = numElements;
    m2nEtaKrnl.execute();

    // adjust: QEtaNodal
    device::aux::plasticity::updateQEtaNodal(qEtaNodalPtrs,
                                             nodalStressTensors,
                                             timeStepWidth,
                                             isAdjustableVector,
                                             numElements,
                                             defaultStream);

    // Convert nodal to modal
    static_assert(kernel::gpu_plConvertEtaNodal2Modal<Cfg>::TmpMaxMemRequiredInBytes == 0);
    kernel::gpu_plConvertEtaNodal2Modal<Cfg> n2mEtaKrnl;
    n2mEtaKrnl.vInv = global->vandermondeMatrixInverse;
    n2mEtaKrnl.QEtaNodal = const_cast<const real**>(qEtaNodalPtrs);
    n2mEtaKrnl.QEtaModal = pstrains;
    n2mEtaKrnl.extraOffset_QEtaModal = tensor::QStress<Cfg>::size();
    n2mEtaKrnl.streamPtr = defaultStream;
    n2mEtaKrnl.flags = isAdjustableVector;
    n2mEtaKrnl.numElements = numElements;
    n2mEtaKrnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif // ACL_DEVICE
}

template<typename Cfg>
void Plasticity<Cfg>::flopsPlasticity(std::uint64_t& nonZeroFlopsCheck,
                                 std::uint64_t& hardwareFlopsCheck,
                                 std::uint64_t& nonZeroFlopsYield,
                                 std::uint64_t& hardwareFlopsYield) {
  // reset flops
  nonZeroFlopsCheck = 0;
  hardwareFlopsCheck = 0;
  nonZeroFlopsYield = 0;
  hardwareFlopsYield = 0;

  // flops from checking, i.e. outside if (adjust) {}
  nonZeroFlopsCheck += kernel::plConvertToNodal<Cfg>::NonZeroFlops;
  hardwareFlopsCheck += kernel::plConvertToNodal<Cfg>::HardwareFlops;

  // compute mean stress
  nonZeroFlopsCheck += kernel::plComputeMean<Cfg>::NonZeroFlops;
  hardwareFlopsCheck += kernel::plComputeMean<Cfg>::HardwareFlops;

  // subtract mean stress
  nonZeroFlopsCheck += kernel::plSubtractMean<Cfg>::NonZeroFlops;
  hardwareFlopsCheck += kernel::plSubtractMean<Cfg>::HardwareFlops;

  // compute second invariant
  nonZeroFlopsCheck += kernel::plComputeSecondInvariant<Cfg>::NonZeroFlops;
  hardwareFlopsCheck += kernel::plComputeSecondInvariant<Cfg>::HardwareFlops;

  // compute taulim (1 add, 1 mul, max NOT counted)
  nonZeroFlopsCheck += static_cast<std::uint64_t>(2 * tensor::meanStress<Cfg>::size());
  hardwareFlopsCheck += static_cast<std::uint64_t>(2 * tensor::meanStress<Cfg>::size());

  // check for yield (NOT counted, as it would require counting the number of yielding points)

  // flops from plastic yielding, i.e. inside if (adjust) {}
  nonZeroFlopsYield += kernel::plAdjustStresses<Cfg>::NonZeroFlops;
  hardwareFlopsYield += kernel::plAdjustStresses<Cfg>::HardwareFlops;
}

#define _H_(cfg) template class Plasticity<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels
