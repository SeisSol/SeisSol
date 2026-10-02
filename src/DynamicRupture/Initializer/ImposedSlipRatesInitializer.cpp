// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ImposedSlipRatesInitializer.h"

#include "Common/ConfigDispatch.h"
#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/ParameterDB.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::dr::initializer {
void ImposedSlipRatesInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  logQuadratureRules(drStorage);
  for (auto& layer : drStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
      // parameters to be read from fault parameters yaml file
      std::unordered_map<std::string, void*> parameterToStorageMap;

      auto* imposedSlipDirection1 = layer.var<LTSImposedSlipRates::ImposedSlipDirection1>(Cfg());
      auto* imposedSlipDirection2 = layer.var<LTSImposedSlipRates::ImposedSlipDirection2>(Cfg());
      auto* onsetTime = layer.var<LTSImposedSlipRates::OnsetTime>(Cfg());

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

      for (std::size_t i = 0; i < Cfg::NumSimulations; ++i) {
        seissol::initializer::FaultParameterDB<real> faultParameterDB(i, Cfg::NumSimulations);
        // read parameters from yaml file
        for (const auto& parameterStoragePair : parameterToStorageMap) {
          faultParameterDB.addParameter(parameterStoragePair.first,
                                        static_cast<real*>(parameterStoragePair.second));
        }
        const auto faceIDs = getFaceIDsInIterator(layer);
        queryModel<Cfg>(faultParameterDB, faceIDs, i);
      }

      rotateSlipToFaultCS<Cfg>(
          layer, strikeSlip, dipSlip, imposedSlipDirection1, imposedSlipDirection2);

      const auto sourceCount = stressSourceCount(*drParameters_);
      auto* stressInFaultCS = layer.var<DynamicRupture::StressSourceInFaultCS>(Cfg());
      auto* pressure = layer.var<DynamicRupture::StressSourcePressure>(Cfg());
      auto* onset = layer.var<DynamicRupture::StressSourceOnset>(Cfg());
      auto* riseTime = layer.var<DynamicRupture::StressSourceRiseTime>(Cfg());
      for (std::uint32_t source = 0; source < sourceCount; ++source) {
        for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
          for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>;
               ++pointIndex) {
            for (std::uint32_t dim = 0; dim < 6; ++dim) {
              stressInFaultCS[ltsFace * sourceCount + source][dim][pointIndex] = 0;
            }
            pressure[ltsFace * sourceCount + source][pointIndex] = 0;
            onset[ltsFace * sourceCount + source][pointIndex] = 0;
            riseTime[ltsFace * sourceCount + source][pointIndex] = 0;
          }
        }
      }

      // Set initial and nucleation stress to zero, these are not needed for this FL

      fixInterpolatedSTFParameters(layer);

      initializeOtherVariables(layer);
    });
  }
}

template <typename Cfg>
void ImposedSlipRatesInitializer::rotateSlipToFaultCS(
    DynamicRupture::Layer& layer,
    const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& strikeSlip,
    const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& dipSlip,
    Real<Cfg> (*imposedSlipDirection1)[misc::NumPaddedPoints<Cfg>],
    Real<Cfg> (*imposedSlipDirection2)[misc::NumPaddedPoints<Cfg>]) {
  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
    const auto meshFace = drFaceInformation[ltsFace].meshFace;
    const Fault& fault = seissolInstance_.meshReader().getFault().at(meshFace);

    CoordinateT strikeVector{};
    CoordinateT dipVector{};
    misc::computeStrikeAndDipVectors(fault.normal, strikeVector, dipVector);

    // cos^2 can be greater than 1 because of rounding errors
    const double cos = std::clamp(MeshTools::dot(strikeVector, fault.tangent1), -1.0, 1.0);
    CoordinateT crossProduct{};
    MeshTools::cross(strikeVector, fault.tangent1, crossProduct);
    const double scalarProduct = MeshTools::dot(crossProduct, fault.normal);
    const double sin = std::sqrt(1 - cos * cos) * std::copysign(1.0, scalarProduct);
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      imposedSlipDirection1[ltsFace][pointIndex] =
          cos * strikeSlip[ltsFace][pointIndex] + sin * dipSlip[ltsFace][pointIndex];
      imposedSlipDirection2[ltsFace][pointIndex] =
          -sin * strikeSlip[ltsFace][pointIndex] + cos * dipSlip[ltsFace][pointIndex];
    }
  }
}

void ImposedSlipRatesInitializer::fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) {
  // do nothing
}

void ImposedSlipRatesYoffeInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*tauS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauS>(Cfg());
    real(*tauR)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauR>(Cfg());
    parameterToStorageMap.insert({"tau_S", reinterpret_cast<real*>(tauS)});
    parameterToStorageMap.insert({"tau_R", reinterpret_cast<real*>(tauR)});
  });
}

void ImposedSlipRatesYoffeInitializer::fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*tauS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauS>(Cfg());
    real(*tauR)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauR>(Cfg());
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
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*riseTime)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSImposedSlipRatesGaussian::RiseTime>(Cfg());
    parameterToStorageMap.insert({"rupture_rise_time", reinterpret_cast<real*>(riseTime)});
  });
}

void ImposedSlipRatesDeltaInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {}
} // namespace seissol::dr::initializer
