// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
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
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::initializer {
void ImposedSlipRatesInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  logInfo() << "Initializing Fault, using a quadrature rule with " << misc::NumBoundaryGaussPoints
            << " points.";
  seissol::initializer::FaultParameterDB faultParameterDB;

  for (auto& layer : dynRupTree->leaves(Ghost)) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSImposedSlipRates*>(dynRup);
    auto* imposedSlipDirection1 = layer.var(concreteLts->imposedSlipDirection1);
    auto* imposedSlipDirection2 = layer.var(concreteLts->imposedSlipDirection2);
    auto* onsetTime = layer.var(concreteLts->onsetTime);

    // First read slip in strike/dip direction. Later we will rotate this to the face aligned
    // coordinate system.
    using VectorOfArraysT = std::vector<std::array<real, misc::NumPaddedPoints>>;
    VectorOfArraysT strikeSlip(layer.getNumberOfCells());
    VectorOfArraysT dipSlip(layer.getNumberOfCells());
    parameterToStorageMap.insert({"strike_slip", strikeSlip.data()->data()});
    parameterToStorageMap.insert({"dip_slip", dipSlip.data()->data()});
    parameterToStorageMap.insert({"rupture_onset", reinterpret_cast<real*>(onsetTime)});

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, layer);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, layer);
    queryModel(faultParameterDB, faceIDs);

    rotateSlipToFaultCS(
        dynRup, layer, strikeSlip, dipSlip, imposedSlipDirection1, imposedSlipDirection2);

    real(*nucleationStressInFaultCS)[misc::NumPaddedPoints][6] =
        layer.var(dynRup->nucleationStressInFaultCS);
    real(*initialStressInFaultCS)[misc::NumPaddedPoints][6] =
        layer.var(dynRup->initialStressInFaultCS);

    // Set initial and nucleation stress to zero, these are not needed for this FL
    for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        for (unsigned int dim = 0; dim < 6; ++dim) {
          initialStressInFaultCS[ltsFace][pointIndex][dim] = 0;
          nucleationStressInFaultCS[ltsFace][pointIndex][dim] = 0;
        }
      }
    }

    fixInterpolatedSTFParameters(dynRup, layer);

    initializeOtherVariables(dynRup, layer);
  }
}

void ImposedSlipRatesInitializer::rotateSlipToFaultCS(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer,
    const std::vector<std::array<real, misc::NumPaddedPoints>>& strikeSlip,
    const std::vector<std::array<real, misc::NumPaddedPoints>>& dipSlip,
    real (*imposedSlipDirection1)[misc::NumPaddedPoints],
    real (*imposedSlipDirection2)[misc::NumPaddedPoints]) {
  for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = layer.var(dynRup->faceInformation);
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissolInstance.meshReader().getFault().at(meshFace);

    VrtxCoords strikeVector{};
    VrtxCoords dipVector{};
    misc::computeStrikeAndDipVectors(fault.normal, strikeVector, dipVector);

    // cos^2 can be greater than 1 because of rounding errors
    const real cos = std::clamp(MeshTools::dot(strikeVector, fault.tangent1), -1.0, 1.0);
    VrtxCoords crossProduct{};
    MeshTools::cross(strikeVector, fault.tangent1, crossProduct);
    const real scalarProduct = MeshTools::dot(crossProduct, fault.normal);
    const real sin = std::sqrt(1 - cos * cos) * std::copysign(1.0, scalarProduct);
    for (size_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      imposedSlipDirection1[ltsFace][pointIndex] =
          cos * strikeSlip[ltsFace][pointIndex] + sin * dipSlip[ltsFace][pointIndex];
      imposedSlipDirection2[ltsFace][pointIndex] =
          -sin * strikeSlip[ltsFace][pointIndex] + cos * dipSlip[ltsFace][pointIndex];
    }
  }
}

void ImposedSlipRatesInitializer::fixInterpolatedSTFParameters(
    const seissol::initializer::DynamicRupture* const dynRup, seissol::initializer::Layer& layer) {
  // do nothing
}

void ImposedSlipRatesYoffeInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesYoffe*>(dynRup);
  real(*tauS)[misc::NumPaddedPoints] = layer.var(concreteLts->tauS);
  real(*tauR)[misc::NumPaddedPoints] = layer.var(concreteLts->tauR);
  parameterToStorageMap.insert({"tau_S", reinterpret_cast<real*>(tauS)});
  parameterToStorageMap.insert({"tau_R", reinterpret_cast<real*>(tauR)});
}

void ImposedSlipRatesYoffeInitializer::fixInterpolatedSTFParameters(
    const seissol::initializer::DynamicRupture* const dynRup, seissol::initializer::Layer& layer) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesYoffe*>(dynRup);
  real(*tauS)[misc::NumPaddedPoints] = layer.var(concreteLts->tauS);
  real(*tauR)[misc::NumPaddedPoints] = layer.var(concreteLts->tauR);
  // ensure that tauR is larger than tauS and that tauS and tauR are greater than 0 (the contrary
  // can happen due to ASAGI interpolation)
  for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
    for (unsigned int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      tauS[ltsFace][pointIndex] = std::max(static_cast<real>(0.0), tauS[ltsFace][pointIndex]);
      tauR[ltsFace][pointIndex] = std::max(tauR[ltsFace][pointIndex], tauS[ltsFace][pointIndex]);
    }
  }
}

void ImposedSlipRatesGaussianInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesGaussian*>(dynRup);
  real(*riseTime)[misc::NumPaddedPoints] = layer.var(concreteLts->riseTime);
  parameterToStorageMap.insert({"rupture_rise_time", reinterpret_cast<real*>(riseTime)});
}

void ImposedSlipRatesDeltaInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {}
} // namespace seissol::dr::initializer
