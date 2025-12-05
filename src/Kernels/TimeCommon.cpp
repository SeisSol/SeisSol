// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "TimeCommon.h"

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Parallel/Runtime/Stream.h"

#include <Common/ConfigHelper.h>
#include <Common/Constants.h>
#include <GeneratedCode/init.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/CellLocalInformation.h>
#include <Initializer/LtsSetup.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/MultipleSimulations.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <stdint.h>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

#ifndef NDEBUG
#include "Alignment.h"

#include <cstdint>
#endif

namespace seissol::kernels {

namespace {
template <typename Cfg, typename CfgNeighbor>
void checkCompatible() {
  if constexpr (!(tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg>] ==
                      tensor::I<CfgNeighbor>::Shape[multisim::BasisDim<CfgNeighbor>] &&
                  tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg> + 1] ==
                      tensor::I<CfgNeighbor>::Shape[multisim::BasisDim<CfgNeighbor> + 1] &&
                  multisim::NumSimulations<Cfg> == multisim::NumSimulations<CfgNeighbor>)) {
    logError() << "Fatal error: wanted to compare differently-sized buffers.";
  }
}

template <typename Cfg, typename CfgNeighbor>
void copyView(Real<Cfg>* ownPtr, Real<CfgNeighbor>* neighborPtr) {
  // convert precision / copy and zero relevant data

  auto ownView = init::I<Cfg>::view::create(ownPtr);
  auto neighborView = init::I<CfgNeighbor>::view::create(neighborPtr);

  using real = Real<Cfg>;

  const auto dataRange = std::min(ownView.shape(multisim::BasisDim<Cfg>),
                                  neighborView.shape(multisim::BasisDim<CfgNeighbor>));

#pragma omp simd collapse(3)
  for (std::size_t j = 0; j < dataRange; ++j) {
    for (std::size_t k = 0; k < ownView.shape(multisim::BasisDim<Cfg> + 1); ++k) {
      for (std::size_t i = 0; i < multisim::NumSimulations<Cfg>; ++i) {
        multisim::multisimWrap<Cfg>(ownView, i, j, k) =
            static_cast<real>(multisim::multisimWrap<CfgNeighbor>(neighborView, i, j, k));
      }
    }
  }
#pragma omp simd collapse(3)
  for (std::size_t j = dataRange; j < ownView.shape(multisim::BasisDim<Cfg>); ++j) {
    for (std::size_t k = 0; k < ownView.shape(multisim::BasisDim<Cfg> + 1); ++k) {

      for (std::size_t i = 0; i < multisim::NumSimulations<Cfg>; ++i) {
        multisim::multisimWrap<Cfg>(ownView, i, j, k) = 0;
      }
    }
  }
}
} // namespace

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
      const auto neighborConfig = cellInfo.neighborConfigIds[dofneighbor];

      // check if the time integration is already done (-> copy pointer)
      if (!ltsSetup.neighborHasDerivatives(dofneighbor)) {
        if (ConfigVariant(Cfg()).index() == neighborConfig) {
          timeIntegrated[dofneighbor] = static_cast<real*>(timeDofs[dofneighbor]);
        } else {
          // we might need to convert at least precision-wise
          std::visit(
              [&](auto cfg) {
                using CfgNeighbor = decltype(cfg);

                checkCompatible<Cfg, CfgNeighbor>();

                using RealNeighbor = Real<CfgNeighbor>;
                auto* neighborPtr = static_cast<RealNeighbor*>(timeDofs[dofneighbor]);

                if constexpr (std::is_same_v<real, RealNeighbor> &&
                              Cfg::ConvergenceOrder <= CfgNeighbor::ConvergenceOrder) {
                  // same precision; just assign the pointer (also works for lower orders)
                  timeIntegrated[dofneighbor] = neighborPtr;
                } else {
                  // convert precision / copy and zero relevant data

                  copyView<Cfg, CfgNeighbor>(integrationBuffer[dofneighbor], neighborPtr);
                  timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
                }
              },
              ConfigVariantList[neighborConfig]);
        }
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        // select the time coefficients next: dependent on if we have GTS as a neighbor, use GTS
        // coefficients (timeCoeffs); otherwise use LTS coefficients (subtimeCoeffs) for a large
        // time cluster neighbor. IMPORTANT: make sure to not land in this code path if the neighbor
        // time cluster is less than the local time cluster. It shouldn't happen with the current
        // setup; but just be aware of it when changing things. In that case, enforce the "GTS
        // relation" instead; then everything will work again.

        const auto* coeffs = ltsSetup.neighborGTSRelation(dofneighbor) ? timeCoeffs : subtimeCoeffs;

        if (ConfigVariant(Cfg()).index() == neighborConfig) {
          time.evaluate(
              coeffs, static_cast<real*>(timeDofs[dofneighbor]), integrationBuffer[dofneighbor]);

          timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
        } else {
          std::visit(
              [&](auto cfg) {
                using CfgNeighbor = decltype(cfg);

                checkCompatible<Cfg, CfgNeighbor>();

                using RealNeighbor = Real<CfgNeighbor>;

                alignas(Alignment) RealNeighbor bufferN[tensor::I<CfgNeighbor>::size()];
                RealNeighbor coeffsN[CfgNeighbor::ConvergenceOrder]{};

                for (std::size_t i = 0;
                     i < std::min(CfgNeighbor::ConvergenceOrder, Cfg::ConvergenceOrder);
                     ++i) {
                  coeffsN[i] = coeffs[i];
                }

                Time<CfgNeighbor> timeNeighbor;

                RealNeighbor* buffer = nullptr;

                if constexpr (std::is_same_v<real, RealNeighbor> &&
                              Cfg::ConvergenceOrder == CfgNeighbor::ConvergenceOrder) {
                  buffer = integrationBuffer[dofneighbor];
                } else {
                  // use temporary buffer
                  buffer = bufferN;
                }

                timeNeighbor.evaluate(
                    coeffsN, static_cast<RealNeighbor*>(timeDofs[dofneighbor]), buffer);

                if constexpr (!(std::is_same_v<real, RealNeighbor> &&
                                Cfg::ConvergenceOrder == CfgNeighbor::ConvergenceOrder)) {
                  copyView<Cfg, CfgNeighbor>(static_cast<real*>(integrationBuffer[dofneighbor]),
                                             buffer);
                }

                timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
              },
              ConfigVariantList[neighborConfig]);
        }
      }
    }
  }
}

template <typename Cfg>
void TimeCommon<Cfg>::computeBatchedIntegrals(
    SEISSOL_GPU_PARAM Time<Cfg>& time,
    SEISSOL_GPU_PARAM const real* timeCoeffs,
    SEISSOL_GPU_PARAM const real* subtimeCoeffs,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  // Compute time integrated dofs using neighbors derivatives using the GTS relation,
  // i.e. the expansion point is around 'timeStepStart'
  ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.evaluateBatched(
        timeCoeffs,
        const_cast<const real**>(
            (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtrAs<real*>()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtrAs<real*>(),
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
        const_cast<const real**>(
            (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtrAs<real*>()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtrAs<real*>(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

#define SEISSOL_CONFIGITER(cfg) template class TimeCommon<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels
