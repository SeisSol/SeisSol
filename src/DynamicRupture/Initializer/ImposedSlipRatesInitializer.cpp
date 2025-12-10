// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ImposedSlipRatesInitializer.h"

#include "DynamicRupture/Misc.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/ParameterDB.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <Geometry/MeshReader.h>
#include <Solver/MultipleSimulations.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::initializer {

namespace {
/**
 * Rotate slip from strike/dip cooordinate system to the fault aligned coordinate system.
 * @param it
 * @param strikeSlip: Slip in strike direction
 * @param dipSlip: Slip in dip direction
 * @param imposedSlipDirection1: Slip in fault aligned direction 1
 * @param imposedSlipDirection2: Slip in fault aligned direction 2
 */
template <typename Cfg>
void rotateSlipToFaultCS(
    DynamicRupture::Layer& layer,
    const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& strikeSlip,
    const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& dipSlip,
    Real<Cfg> (*imposedSlipDirection1)[misc::NumPaddedPoints<Cfg>],
    Real<Cfg> (*imposedSlipDirection2)[misc::NumPaddedPoints<Cfg>],
    const geometry::MeshReader& meshReader) {
  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
    const auto meshFace = drFaceInformation[ltsFace].meshFace;
    const Fault& fault = meshReader.getFault().at(meshFace);

    VrtxCoords strikeVector{};
    VrtxCoords dipVector{};
    misc::computeStrikeAndDipVectors(fault.normal, strikeVector, dipVector);

    // cos^2 can be greater than 1 because of rounding errors
    const auto cos = std::clamp(MeshTools::dot(strikeVector, fault.tangent1), -1.0, 1.0);
    VrtxCoords crossProduct{};
    MeshTools::cross(strikeVector, fault.tangent1, crossProduct);
    const auto scalarProduct = MeshTools::dot(crossProduct, fault.normal);
    const auto sin = std::sqrt(1 - cos * cos) * std::copysign(1.0, scalarProduct);
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      imposedSlipDirection1[ltsFace][pointIndex] =
          cos * strikeSlip[ltsFace][pointIndex] + sin * dipSlip[ltsFace][pointIndex];
      imposedSlipDirection2[ltsFace][pointIndex] =
          -sin * strikeSlip[ltsFace][pointIndex] + cos * dipSlip[ltsFace][pointIndex];
    }
  }
}
} // namespace

void ImposedSlipRatesInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  logInfo() << "Initializing Fault...";
  for (auto& layer : drStorage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>;
      // parameters to be read from fault parameters yaml file
      std::unordered_map<std::string, void*> parameterToStorageMap;

      auto* imposedSlipDirection1 = layer.var<LTSImposedSlipRates::ImposedSlipDirection1>(cfg);
      auto* imposedSlipDirection2 = layer.var<LTSImposedSlipRates::ImposedSlipDirection2>(cfg);
      auto* onsetTime = layer.var<LTSImposedSlipRates::OnsetTime>(cfg);

      // First read slip in strike/dip direction. Later we will rotate this to the face aligned
      // coordinate system.
      using VectorOfArraysT = std::vector<std::array<real, misc::NumPaddedPoints<Cfg>>>;
      VectorOfArraysT strikeSlip(layer.size());
      VectorOfArraysT dipSlip(layer.size());
      parameterToStorageMap.insert({"strike_slip", strikeSlip.data()->data()});
      parameterToStorageMap.insert({"dip_slip", dipSlip.data()->data()});
      parameterToStorageMap.insert({"rupture_onset", reinterpret_cast<real*>(onsetTime)});

      // get additional parameters (for derived friction laws)
      addAdditionalParameters(parameterToStorageMap, layer);

      for (std::size_t i = 0; i < multisim::NumSimulations<Cfg>; ++i) {
        seissol::initializer::FaultParameterDB<Real<Cfg>> faultParameterDB(
            i, multisim::NumSimulations<Cfg>);
        // read parameters from yaml file
        for (const auto& parameterStoragePair : parameterToStorageMap) {
          faultParameterDB.addParameter(parameterStoragePair.first,
                                        reinterpret_cast<real*>(parameterStoragePair.second));
        }
        const auto faceIDs = getFaceIDsInIterator(layer);
        queryModel(faultParameterDB, faceIDs, i, layer.getIdentifier().config.index());
      }

      rotateSlipToFaultCS<Cfg>(layer,
                               strikeSlip,
                               dipSlip,
                               imposedSlipDirection1,
                               imposedSlipDirection2,
                               seissolInstance.meshReader());

      auto* initialStressInFaultCS = layer.var<DynamicRupture::InitialStressInFaultCS>(cfg);
      auto* initialPressure = layer.var<DynamicRupture::InitialPressure>(cfg);
      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
          for (unsigned int dim = 0; dim < 6; ++dim) {
            initialStressInFaultCS[ltsFace][dim][pointIndex] = 0;
          }
          initialPressure[ltsFace][pointIndex] = 0;
        }
      }

      for (unsigned i = 0; i < drParameters->nucleationCount; ++i) {
        auto* nucleationStressInFaultCS = layer.var<DynamicRupture::NucleationStressInFaultCS>(cfg);
        auto* nucleationPressure = layer.var<DynamicRupture::NucleationPressure>(cfg);
        for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
          for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>;
               ++pointIndex) {
            for (unsigned int dim = 0; dim < 6; ++dim) {
              nucleationStressInFaultCS[ltsFace * drParameters->nucleationCount + i][dim]
                                       [pointIndex] = 0;
            }
            nucleationPressure[ltsFace * drParameters->nucleationCount + i][pointIndex] = 0;
          }
        }
      }

      // Set initial and nucleation stress to zero, these are not needed for this FL

      fixInterpolatedSTFParameters(layer);

      initializeOtherVariables(layer);
    });
  }
}

void ImposedSlipRatesInitializer::fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) {
  // do nothing
}

void ImposedSlipRatesYoffeInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;
    real(*tauS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauS>(cfg);
    real(*tauR)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauR>(cfg);
    parameterToStorageMap.insert({"tau_S", reinterpret_cast<real*>(tauS)});
    parameterToStorageMap.insert({"tau_R", reinterpret_cast<real*>(tauR)});
  });
}

void ImposedSlipRatesYoffeInitializer::fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) {
  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;
    real(*tauS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauS>(cfg);
    real(*tauR)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauR>(cfg);
    // ensure that tauR is larger than tauS and that tauS and tauR are greater than 0 (the contrary
    // can happen due to ASAGI interpolation)
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
        tauS[ltsFace][pointIndex] = std::max(static_cast<real>(0.0), tauS[ltsFace][pointIndex]);
        tauR[ltsFace][pointIndex] = std::max(tauR[ltsFace][pointIndex], tauS[ltsFace][pointIndex]);
      }
    }
  });
}

void ImposedSlipRatesGaussianInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;
    real(*riseTime)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSImposedSlipRatesGaussian::RiseTime>(cfg);
    parameterToStorageMap.insert({"rupture_rise_time", reinterpret_cast<real*>(riseTime)});
  });
}

void ImposedSlipRatesDeltaInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {}
} // namespace seissol::dr::initializer
