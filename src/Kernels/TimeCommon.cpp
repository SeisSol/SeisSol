// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "TimeCommon.h"
#include <Common/ConfigHelper.h>
#include <Common/Constants.h>
#include <Config.h>
#include <DataTypes/ConditionalTable.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/LtsSetup.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/MultipleSimulations.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <stdint.h>

#include "utils/logger.h"

#ifdef ACL_DEVICE
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#endif

#ifndef NDEBUG
#include "Alignment.h"
#include <cstdint>
#endif

namespace seissol::kernels {

template <typename Cfg>
void TimeCommon<Cfg>::computeIntegrals(Time<Cfg>& time,
                                       const CellLocalInformation& cellInfo,
                                       const real* timeCoeffs,
                                       const real* subtimeCoeffs,
                                       void* const timeDofs[4],
                                       real integrationBuffer[4][tensor::I<Cfg>::size()],
                                       real* timeIntegrated[4]) {
  // call the more general assembly
  /*
   * assert valid input.
   */

  const auto& faceTypes = cellInfo.faceTypes;
  const auto& ltsSetup = cellInfo.ltsSetup;

#ifndef NDEBUG
  // alignment of the time derivatives/integrated dofs and the buffer
  for (std::size_t dofneighbor = 0; dofneighbor < Cell::NumFaces; dofneighbor++) {
    assert(reinterpret_cast<uintptr_t>(timeDofs[dofneighbor]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(integrationBuffer[dofneighbor]) % Alignment == 0);
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for (std::size_t dofneighbor = 0; dofneighbor < Cell::NumFaces; ++dofneighbor) {
    // collect information only in the case that neighboring element contributions are required
    if (faceTypes[dofneighbor] != FaceType::Outflow &&
        faceTypes[dofneighbor] != FaceType::DynamicRupture) {
      // check if the time integration is already done (-> copy pointer)
      if (!ltsSetup.neighborHasDerivatives(dofneighbor)) {
        const auto neighborConfig = cellInfo.neighborConfigIds[dofneighbor];

        if (ConfigVariant(Cfg()).index() == neighborConfig) {
          timeIntegrated[dofneighbor] = static_cast<real*>(timeDofs[dofneighbor]);
        } else {
          // we might need to convert at least precision-wise
          std::visit(
              [&](auto cfg) {
                using CfgNeighbor = decltype(cfg);

                if constexpr (!(tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg>] ==
                                    tensor::I<
                                        CfgNeighbor>::Shape[multisim::BasisDim<CfgNeighbor>] &&
                                tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg> + 1] ==
                                    tensor::I<CfgNeighbor>::Shape[multisim::BasisDim<CfgNeighbor> +
                                                                  1] &&
                                multisim::NumSimulations<Cfg> ==
                                    multisim::NumSimulations<CfgNeighbor>)) {
                  logError() << "Fatal error: wanted to compare differently-sized buffers.";
                }

                using RealNeighbor = Real<CfgNeighbor>;
                auto* neighborPtr = static_cast<RealNeighbor*>(timeDofs[dofneighbor]);

                if constexpr (std::is_same_v<real, RealNeighbor>) {
                  // same precision; just assign the pointer
                  timeIntegrated[dofneighbor] = neighborPtr;
                } else {
                  // convert precision

                  auto ownView = init::I<Cfg>::view::create(integrationBuffer[dofneighbor]);
                  auto neighborView = init::I<CfgNeighbor>::view::create(neighborPtr);

#pragma omp simd collapse(3)
                  for (std::size_t i = 0; i < multisim::NumSimulations<Cfg>; ++i) {
                    for (std::size_t j = 0; j < ownView.shape(multisim::BasisDim<Cfg>); ++j) {
                      for (std::size_t k = 0; k < ownView.shape(multisim::BasisDim<Cfg> + 1); ++k) {
                        multisim::multisimWrap<Cfg>(ownView, i, j, k) = static_cast<real>(
                            multisim::multisimWrap<CfgNeighbor>(neighborView, i, j, k));
                      }
                    }
                  }
                  timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
                }
              },
              ConfigVariantList[neighborConfig]);
        }
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        const auto* coeffs = ltsSetup.neighborGTS(dofneighbor) ? timeCoeffs : subtimeCoeffs;
        time.evaluate(
            coeffs, static_cast<real*>(timeDofs[dofneighbor]), integrationBuffer[dofneighbor]);

        timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
      }
    }
  }
}

template <typename Cfg>
void TimeCommon<Cfg>::computeBatchedIntegrals(Time<Cfg>& time,
                                              const real* timeCoeffs,
                                              const real* subtimeCoeffs,
                                              ConditionalPointersToRealsTable& table,
                                              seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // Compute time integrated dofs using neighbors derivatives using the GTS relation,
  // i.e. the expansion point is around 'timeStepStart'
  ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.evaluateBatched(
        timeCoeffs,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }

  // Compute time integrated dofs using neighbors derivatives using the LTS relation,
  // i.e. the expansion point is around '0'
  key = ConditionalKey(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.evaluateBatched(
        subtimeCoeffs,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

#define _H_(cfg) template class TimeCommon<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels
