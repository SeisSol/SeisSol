#include "ImposedSlipRatesInitializer.h"

#include <utils/logger.h>
#include "Model/common.hpp"
#include "SeisSol.h"

namespace seissol::dr::initializers {
void ImposedSlipRatesInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                                  seissol::initializers::LTSTree* dynRupTree,
                                                  seissol::Interoperability* eInteroperability) {
  logInfo() << "Initializing Fault, using a quadrature rule with "
            << misc::numberOfBoundaryGaussPoints << " points.";
  seissol::initializers::FaultParameterDB faultParameterDB;

  for (auto it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRates*>(dynRup);
    real(*strikeSlip)[misc::numPaddedPoints] = it->var(concreteLts->strikeSlip);
    real(*dipSlip)[misc::numPaddedPoints] = it->var(concreteLts->dipSlip);
    real(*onsetTime)[misc::numPaddedPoints] = it->var(concreteLts->onsetTime);

    parameterToStorageMap.insert({"strike_slip", (real*)strikeSlip});
    parameterToStorageMap.insert({"dip_slip", (real*)dipSlip});
    parameterToStorageMap.insert({"rupture_onset", (real*)onsetTime});

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, it);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, it);
    queryModel(faultParameterDB, faceIDs);

    rotateSlipToFaultCS(dynRup, it, strikeSlip, dipSlip);

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

    ensureCorrectness(dynRup, it);

    initializeOtherVariables(dynRup, it, eInteroperability);
  }
}
void ImposedSlipRatesInitializer::rotateSlipToFaultCS(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    real (*strikeSlip)[misc::numPaddedPoints],
    real (*dipSlip)[misc::numPaddedPoints]) {
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissol::SeisSol::main.meshReader().getFault().at(meshFace);

    VrtxCoords strikeVector{};
    VrtxCoords dipVector{};
    misc::computeStrikeAndDipVectors(fault.normal, strikeVector, dipVector);

    // cos^2 can be greater than 1 because of rounding errors
    real cos = std::min(1.0, MeshTools::dot(strikeVector, fault.tangent1));
    VrtxCoords crossProduct{};
    MeshTools::cross(strikeVector, fault.tangent1, crossProduct);
    real scalarProduct = MeshTools::dot(crossProduct, fault.normal);
    real sin = std::sqrt(1 - cos * cos) * std::copysign(1.0, scalarProduct);
    for (size_t pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      real tmpStrikeSlip = strikeSlip[ltsFace][pointIndex];
      real tmpDipSlip = dipSlip[ltsFace][pointIndex];
      strikeSlip[ltsFace][pointIndex] = cos * tmpStrikeSlip + sin * tmpDipSlip;
      dipSlip[ltsFace][pointIndex] = -sin * tmpStrikeSlip + cos * tmpDipSlip;
    }
  }
}

void ImposedSlipRatesYoffeInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesYoffe*>(dynRup);
  real(*tauS)[misc::numPaddedPoints] = it->var(concreteLts->tauS);
  real(*tauR)[misc::numPaddedPoints] = it->var(concreteLts->tauR);
  parameterToStorageMap.insert({"tau_S", (real*)tauS});
  parameterToStorageMap.insert({"tau_R", (real*)tauR});
}

void ImposedSlipRatesInitializer::ensureCorrectness(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // do nothing
}

void ImposedSlipRatesYoffeInitializer::ensureCorrectness(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesYoffe*>(dynRup);
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
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesGaussian*>(dynRup);
  real(*riseTime)[misc::numPaddedPoints] = it->var(concreteLts->riseTime);
  parameterToStorageMap.insert({"rupture_rise_time", (real*)riseTime});
}
} // namespace seissol::dr::initializers