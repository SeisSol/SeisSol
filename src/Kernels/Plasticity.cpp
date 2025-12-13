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
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Model/Plasticity.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "DeviceAux/PlasticityAux.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Solver/MultipleSimulations.h"

#include <Device/device.h>
using namespace device;
#endif

#ifndef NDEBUG
#include <cstdint>
#endif

namespace seissol::kernels {
std::size_t Plasticity::computePlasticity(double oneMinusIntegratingFactor,
                                          double timeStepWidth,
                                          double tV,
                                          const GlobalData* global,
                                          const seissol::model::PlasticityData* plasticityData,
                                          real degreesOfFreedom[tensor::Q::size()],
                                          real* pstrain) {

  assert(reinterpret_cast<uintptr_t>(degreesOfFreedom) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(global->vandermondeMatrix) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(global->vandermondeMatrixInverse) % Alignment == 0);

  alignas(Alignment) real qStressNodal[tensor::QStressNodal::size()]{};

  alignas(Alignment) real meanStress[tensor::meanStress::size()]{};
  alignas(Alignment) real secondInvariant[tensor::secondInvariant::size()]{};
  alignas(Alignment) real tau[tensor::secondInvariant::size()]{};
  alignas(Alignment) real taulim[tensor::meanStress::size()]{};
  alignas(Alignment) real yieldFactor[tensor::yieldFactor::size()]{};

  static_assert(tensor::secondInvariant::size() == tensor::meanStress::size(),
                "Second invariant tensor and mean stress tensor must be of the same size().");
  static_assert(tensor::yieldFactor::size() <= tensor::meanStress::size(),
                "Yield factor tensor must be smaller than mean stress tensor.");

  /* Convert modal to nodal and add sigma0.
   * Stores s_{ij} := sigma_{ij} + sigma0_{ij} for every node.
   * sigma0 is constant
   * also stores the previous DOFs before adding the new initial loading
   */

  kernel::plConvertToNodal m2nKrnl;
  m2nKrnl.v = global->vandermondeMatrix;
  m2nKrnl.QStress = degreesOfFreedom;
  m2nKrnl.QStressNodal = qStressNodal;
  m2nKrnl.replicateInitialLoading = init::replicateInitialLoading::Values;
  m2nKrnl.initialLoading = plasticityData->initialLoading;
  m2nKrnl.execute();

  // Computes m = s_{ii} / 3.0 for every node
  kernel::plComputeMean cmKrnl;
  cmKrnl.meanStress = meanStress;
  cmKrnl.QStressNodal = qStressNodal;
  cmKrnl.selectBulkAverage = init::selectBulkAverage::Values;
  cmKrnl.execute();

  /* Compute s_{ij} := s_{ij} - m delta_{ij},
   * where delta_{ij} = 1 if i == j else 0.
   * Thus, s_{ij} contains the deviatoric stresses. */
  kernel::plSubtractMean smKrnl;
  smKrnl.meanStress = meanStress;
  smKrnl.QStressNodal = qStressNodal;
  smKrnl.selectBulkNegative = init::selectBulkNegative::Values;
  smKrnl.execute();

  // Compute I_2 = 0.5 s_{ij} s_ji for every node
  kernel::plComputeSecondInvariant siKrnl;
  siKrnl.secondInvariant = secondInvariant;
  siKrnl.QStressNodal = qStressNodal;
  siKrnl.weightSecondInvariant = init::weightSecondInvariant::Values;
  siKrnl.execute();

// tau := sqrt(I_2) for every node
#pragma omp simd
  for (std::size_t ip = 0; ip < tensor::secondInvariant::size(); ++ip) {
    tau[ip] = std::sqrt(secondInvariant[ip]);
  }

// Compute tau_c for every node
#pragma omp simd
  for (std::size_t ip = 0; ip < tensor::meanStress::size(); ++ip) {
    taulim[ip] = std::max(
        static_cast<real>(0.0),
        plasticityData->cohesionTimesCosAngularFriction[ip % seissol::multisim::NumSimulations] -
            meanStress[ip] *
                plasticityData->sinAngularFriction[ip % seissol::multisim::NumSimulations]);
  }

  int adjust = 0;

#pragma omp simd reduction(max : adjust)
  for (std::size_t ip = 0; ip < tensor::yieldFactor::size(); ++ip) {
    // Compute yield := (t_c / tau - 1) r for every node,
    // where r = 1 - exp(-timeStepWidth / tV)
    if (tau[ip] > taulim[ip]) {
      adjust = 1;
      yieldFactor[ip] = (taulim[ip] / tau[ip] - 1.0) * oneMinusIntegratingFactor;
    } else {
      yieldFactor[ip] = 0.0;
    }
  }

  if (adjust != 0) {
    const real factor = plasticityData->mufactor / (tV * oneMinusIntegratingFactor);

    // calculate plastic strain
    constexpr std::size_t NumNodes = init::QStressNodal::Stop[multisim::BasisFunctionDimension] -
                                     init::QStressNodal::Start[multisim::BasisFunctionDimension];

    real* __restrict qEtaNodal = &pstrain[tensor::QStressNodal::size()];

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

#pragma omp simd collapse(2)
    for (std::size_t i = 0; i < NumNodes; ++i) {
      for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
        const auto qp = s + multisim::NumSimulations * i;

        real dudtPstrainSqAcc = 0;

        for (std::size_t x = 0; x < 6; ++x) {
          const auto q = qp + multisim::NumSimulations * NumNodes * x;
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

        qEtaNodal[qp] += timeStepWidth * std::sqrt(0.5 * dudtPstrainSqAcc);
      }
    }

    kernel::plConvertToModal adjKrnl;
    adjKrnl.QStress = degreesOfFreedom;
    adjKrnl.vInv = global->vandermondeMatrixInverse;
    adjKrnl.QStressNodal = qStressNodal;
    adjKrnl.execute();

    return 1;
  }
  return 0;
}

void Plasticity::computePlasticityBatched(
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM double tV,
    SEISSOL_GPU_PARAM const GlobalData* global,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::model::PlasticityData* plasticityData,
    SEISSOL_GPU_PARAM std::size_t* yieldCounter,
    SEISSOL_GPU_PARAM unsigned* isAdjustableVector,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  static_assert(tensor::Q::Shape[0] == tensor::QStressNodal::Shape[0],
                "modal and nodal dofs must have the same leading dimensions");
  static_assert(tensor::Q::Shape[multisim::BasisFunctionDimension] == tensor::v::Shape[0],
                "modal dofs and vandermonde matrix must have the same leading dimensions");

  ConditionalKey key(*KernelNames::Plasticity);
  auto defaultStream = runtime.stream();

  if (table.find(key) != table.end()) {
    const auto oneMinusIntegratingFactor = computeRelaxTime(tV, timeStepWidth);

    auto& entry = table[key];
    const size_t numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();

    // Convert modal to nodal
    real** modalStressTensors = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    real** nodalStressTensors =
        (entry.get(inner_keys::Wp::Id::NodalStressTensor))->getDeviceDataPtr();

    assert(global->replicateStresses != nullptr && "replicateStresses has not been initialized");
    static_assert(kernel::gpu_plConvertToNodal::TmpMaxMemRequiredInBytes == 0);
    real** initLoad = (entry.get(inner_keys::Wp::Id::InitialLoad))->getDeviceDataPtr();
    kernel::gpu_plConvertToNodal m2nKrnl;
    m2nKrnl.v = global->vandermondeMatrix;
    m2nKrnl.QStress = const_cast<const real**>(modalStressTensors);
    m2nKrnl.QStressNodal = nodalStressTensors;
    m2nKrnl.replicateInitialLoadingM = global->replicateStresses;
    m2nKrnl.initialLoadingM = const_cast<const real**>(initLoad);
    m2nKrnl.streamPtr = defaultStream;
    m2nKrnl.numElements = numElements;
    m2nKrnl.execute();

    real** pstrains = entry.get(inner_keys::Wp::Id::Pstrains)->getDeviceDataPtr();

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

    kernel::gpu_plConvertToModal n2mKrnl;
    n2mKrnl.vInv = global->vandermondeMatrixInverse;
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

void Plasticity::flopsPlasticity(std::uint64_t& nonZeroFlopsCheck,
                                 std::uint64_t& hardwareFlopsCheck,
                                 std::uint64_t& nonZeroFlopsYield,
                                 std::uint64_t& hardwareFlopsYield) {
  // reset flops
  nonZeroFlopsCheck = 0;
  hardwareFlopsCheck = 0;
  nonZeroFlopsYield = 0;
  hardwareFlopsYield = 0;

  // flops from checking, i.e. outside if (adjust) {}
  nonZeroFlopsCheck += kernel::plConvertToNodal::NonZeroFlops;
  hardwareFlopsCheck += kernel::plConvertToNodal::HardwareFlops;

  // compute mean stress
  nonZeroFlopsCheck += kernel::plComputeMean::NonZeroFlops;
  hardwareFlopsCheck += kernel::plComputeMean::HardwareFlops;

  // subtract mean stress
  nonZeroFlopsCheck += kernel::plSubtractMean::NonZeroFlops;
  hardwareFlopsCheck += kernel::plSubtractMean::HardwareFlops;

  // compute second invariant
  nonZeroFlopsCheck += kernel::plComputeSecondInvariant::NonZeroFlops;
  hardwareFlopsCheck += kernel::plComputeSecondInvariant::HardwareFlops;

  // compute taulim (1 add, 1 mul, max NOT counted)
  nonZeroFlopsCheck += static_cast<std::uint64_t>(2 * tensor::meanStress::size());
  hardwareFlopsCheck += static_cast<std::uint64_t>(2 * tensor::meanStress::size());

  // check for yield (NOT counted, as it would require counting the number of yielding points)

  // flops from plastic yielding, i.e. inside if (adjust) {}
  nonZeroFlopsYield += kernel::plConvertToModal::NonZeroFlops;
  hardwareFlopsYield += kernel::plConvertToModal::HardwareFlops;

  // manually counted
  nonZeroFlopsYield += tensor::QStressNodal::size() * 6;
  hardwareFlopsYield += tensor::QStressNodal::size() * 6;
  nonZeroFlopsYield += tensor::QEtaNodal::size() * 3;
  hardwareFlopsYield += tensor::QEtaNodal::size() * 3;
}
} // namespace seissol::kernels
