#include "ImposedSlipRatesInitializer.h"

#include "Model/common.hpp"
#include "SeisSol.h"
#include <utils/logger.h>

namespace seissol::dr::initializers {
void ImposedSlipRatesInitializer::initializeFault(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  const int rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Initializing Fault, using a quadrature rule with "
                << misc::numberOfBoundaryGaussPoints << " points.";
  seissol::initializers::FaultParameterDB faultParameterDB;

  for (auto it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSImposedSlipRates const* const>(dynRup);
    auto* imposedSlipDirection1 = it->var(concreteLts->imposedSlipDirection1);
    auto* imposedSlipDirection2 = it->var(concreteLts->imposedSlipDirection2);
    auto* onsetTime = it->var(concreteLts->onsetTime);

    // First read slip in strike/dip direction. Later we will rotate this to the face aligned
    // coordinate system.
    using VectorOfArraysT = std::vector<std::array<real, misc::numPaddedPoints>>;
    VectorOfArraysT strikeSlip(it->getNumberOfCells());
    VectorOfArraysT dipSlip(it->getNumberOfCells());
    parameterToStorageMap.insert({"strike_slip", strikeSlip.data()->data()});
    parameterToStorageMap.insert({"dip_slip", dipSlip.data()->data()});
    parameterToStorageMap.insert({"rupture_onset", (real*)onsetTime});

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, it);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, it);
    queryModel(faultParameterDB, faceIDs);

    rotateSlipToFaultCS(
        dynRup, it, strikeSlip, dipSlip, imposedSlipDirection1, imposedSlipDirection2);

    real(*nucleationStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(dynRup->nucleationStressInFaultCS);
    real(*initialStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(dynRup->initialStressInFaultCS);

    // Set initial and nucleation stress to zero, these are not needed for this FL
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        for (unsigned int dim = 0; dim < 6; ++dim) {
          initialStressInFaultCS[ltsFace][pointIndex][dim] = 0;
          nucleationStressInFaultCS[ltsFace][pointIndex][dim] = 0;
        }
      }
    }

    fixInterpolatedSTFParameters(dynRup, it);

    initializeOtherVariables(dynRup, it);
  }
}

void ImposedSlipRatesInitializer::rotateSlipToFaultCS(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    std::vector<std::array<real, misc::numPaddedPoints>> const& strikeSlip,
    std::vector<std::array<real, misc::numPaddedPoints>> const& dipSlip,
    real (*imposedSlipDirection1)[misc::numPaddedPoints],
    real (*imposedSlipDirection2)[misc::numPaddedPoints]) {
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissolInstance.meshReader().getFault().at(meshFace);

    VrtxCoords strikeVector{};
    VrtxCoords dipVector{};
    misc::computeStrikeAndDipVectors(fault.normal, strikeVector, dipVector);

    // cos^2 can be greater than 1 because of rounding errors
    real cos = std::clamp(MeshTools::dot(strikeVector, fault.tangent1), -1.0, 1.0);
    VrtxCoords crossProduct{};
    MeshTools::cross(strikeVector, fault.tangent1, crossProduct);
    real scalarProduct = MeshTools::dot(crossProduct, fault.normal);
    real sin = std::sqrt(1 - cos * cos) * std::copysign(1.0, scalarProduct);
    for (size_t pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      imposedSlipDirection1[ltsFace][pointIndex] =
          cos * strikeSlip[ltsFace][pointIndex] + sin * dipSlip[ltsFace][pointIndex];
      imposedSlipDirection2[ltsFace][pointIndex] =
          -sin * strikeSlip[ltsFace][pointIndex] + cos * dipSlip[ltsFace][pointIndex];
    }
  }
}

void ImposedSlipRatesInitializer::fixInterpolatedSTFParameters(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // do nothing
}

void ImposedSlipRatesYoffeInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesYoffe const* const>(dynRup);
  real(*tauS)[misc::numPaddedPoints] = it->var(concreteLts->tauS);
  real(*tauR)[misc::numPaddedPoints] = it->var(concreteLts->tauR);
  parameterToStorageMap.insert({"tau_S", (real*)tauS});
  parameterToStorageMap.insert({"tau_R", (real*)tauR});
}

void ImposedSlipRatesYoffeInitializer::fixInterpolatedSTFParameters(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesYoffe const* const>(dynRup);
  real(*tauS)[misc::numPaddedPoints] = it->var(concreteLts->tauS);
  real(*tauR)[misc::numPaddedPoints] = it->var(concreteLts->tauR);
  // ensure that tauR is larger than tauS and that tauS and tauR are greater than 0 (the contrary
  // can happen due to ASAGI interpolation)
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      tauS[ltsFace][pointIndex] = std::max(static_cast<real>(0.0), tauS[ltsFace][pointIndex]);
      tauR[ltsFace][pointIndex] = std::max(tauR[ltsFace][pointIndex], tauS[ltsFace][pointIndex]);
    }
  }
}

void ImposedSlipRatesGaussianInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesGaussian const* const>(dynRup);
  real(*riseTime)[misc::numPaddedPoints] = it->var(concreteLts->riseTime);
  parameterToStorageMap.insert({"rupture_rise_time", (real*)riseTime});
}
} // namespace seissol::dr::initializers
