#include "RateAndStateInitializer.h"

namespace seissol::dr::initializers {
// void RateAndStateFastVelocityInitializer::initializeFrictionMatrices(
//    seissol::initializers::DynamicRupture* dynRup,
//    seissol::initializers::LTSTree* dynRupTree,
//    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
//    std::unordered_map<std::string, double*> faultParameters,
//    unsigned* ltsFaceToMeshFace,
//    seissol::Interoperability& e_interoperability) {
//  BaseDRInitializer::initializeFrictionMatrices(
//      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
//  auto concreteLts =
//      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);
//
//  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;
//
//  for (seissol::initializers::LTSTree::leaf_iterator it =
//           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
//       it != dynRupTree->endLeaf();
//       ++it) {
//
//    real(*nucleationStressInFaultCS)[numPaddedPoints][6] =
//        it->var(concreteLts->nucleationStressInFaultCS);           // get from fortran
//    real(*RS_sl0)[numPaddedPoints] = it->var(concreteLts->RS_sl0); // get from faultParameters
//    real(*RS_a)[numPaddedPoints] = it->var(concreteLts->RS_a);     // get from faultParameters
//    real(*rs_srW)[numPaddedPoints] = it->var(concreteLts->rs_srW); // get from faultParameters
//    bool(*DS)[numPaddedPoints] = it->var(concreteLts->DS);         // par file
//    real* averaged_Slip = it->var(concreteLts->averagedSlip);      // = 0
//    real(*stateVar)[numPaddedPoints] =
//        it->var(concreteLts->stateVariable); // get from Fortran = EQN%IniStateVar
//    real(*dynStressTime)[numPaddedPoints] = it->var(concreteLts->dynStressTime); // = 0
//
//    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
//      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
//
//      e_interoperability.getDynRupStateVar(ltsFace, meshFace, stateVar);
//      e_interoperability.getDynRupNucStress(ltsFace, meshFace, nucleationStressInFaultCS);
//
//      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints;
//           ++pointIndex) { // loop includes padded elements
//        dynStressTime[ltsFace][pointIndex] = 0.0;
//        DS[ltsFace][pointIndex] = drParameters.isDsOutputOn;
//      }
//      averaged_Slip[ltsFace] = 0.0;
//
//      for (unsigned pointIndex = 0; pointIndex < numberOfPoints; ++pointIndex) {
//        RS_a[ltsFace][pointIndex] =
//            static_cast<real>(faultParameters["rs_a"][(meshFace)*numberOfPoints + pointIndex]);
//        rs_srW[ltsFace][pointIndex] =
//            static_cast<real>(faultParameters["rs_srW"][(meshFace)*numberOfPoints + pointIndex]);
//        RS_sl0[ltsFace][pointIndex] =
//            static_cast<real>(faultParameters["RS_sl0"][(meshFace)*numberOfPoints + pointIndex]);
//      }
//      // initialize padded elements for vectorization
//      for (unsigned pointIndex = numberOfPoints; pointIndex < numPaddedPoints; ++pointIndex) {
//        RS_a[ltsFace][pointIndex] = 0.0;
//        rs_srW[ltsFace][pointIndex] = 0.0;
//        RS_sl0[ltsFace][pointIndex] = 0.0;
//      }
//
//    } // lts-face loop
//    layerLtsFaceToMeshFace += it->getNumberOfCells();
//  } // leaf_iterator loop
//}
//
// void RateAndStateFL103TPInitializer::initializeFrictionMatrices(
//    seissol::initializers::DynamicRupture* dynRup,
//    seissol::initializers::LTSTree* dynRupTree,
//    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
//    std::unordered_map<std::string, double*> faultParameters,
//    unsigned* ltsFaceToMeshFace,
//    seissol::Interoperability& e_interoperability) {
//  // BaseDrInitializer::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters,
//  // ltsFaceToMeshFace, e_interoperability);
//  RateAndStateFL103Initializer::initializeFrictionMatrices(
//      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
//
//  auto concreteLts =
//      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);
//  seissol::dr::friction_law::RateAndStateThermalFL103* SolverFL103 =
//      dynamic_cast<seissol::dr::friction_law::RateAndStateThermalFL103*>(FrictionLaw);
//  SolverFL103->initializeTP(e_interoperability); //
//  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;
//
//  for (seissol::initializers::LTSTree::leaf_iterator it =
//           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
//       it != dynRupTree->endLeaf();
//       ++it) {
//
//    real(*temperature)[numPaddedPoints] = it->var(concreteLts->temperature);
//    real(*pressure)[numPaddedPoints] = it->var(concreteLts->pressure);
//
//    real(*TP_Theta)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_theta);
//    real(*TP_sigma)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_sigma);
//
//    real(*TP_halfWidthShearZone)[numPaddedPoints] =
//        it->var(concreteLts->TP_halfWidthShearZone);
//    real(*alphaHy)[numPaddedPoints] = it->var(concreteLts->alphaHy);
//
//    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
//      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
//
//      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
//        temperature[ltsFace][pointIndex] = drParameters.iniTemp;
//        pressure[ltsFace][pointIndex] = drParameters.iniPressure;
//        TP_halfWidthShearZone[ltsFace][pointIndex] = static_cast<real>(
//            faultParameters["TP_halfWidthShearZone"][(meshFace)*numberOfPoints + pointIndex]);
//        alphaHy[ltsFace][pointIndex] =
//            static_cast<real>(faultParameters["alphaHy"][(meshFace)*numberOfPoints +
//            pointIndex]);
//        for (unsigned iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; ++iTP_grid_nz) {
//          TP_Theta[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
//          TP_sigma[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
//        }
//      }
//    } // lts-face loop
//    layerLtsFaceToMeshFace += it->getNumberOfCells();
//  } // leaf_iterator loop
//}
void RateAndStateInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                              seissol::initializers::LTSTree* dynRupTree,
                                              seissol::Interoperability* e_interoperability) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree, e_interoperability);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    bool(*DS)[numPaddedPoints] = it->var(concreteLts->DS);
    real* averagedSlip = it->var(concreteLts->averagedSlip);
    real(*dynStressTime)[numPaddedPoints] = it->var(concreteLts->dynStressTime);
    real(*slipRateStrike)[numPaddedPoints] = it->var(concreteLts->slipRateStrike);
    real(*slipRateDip)[numPaddedPoints] = it->var(concreteLts->slipRateDip);
    real(*mu)[numPaddedPoints] = it->var(concreteLts->mu);

    real(*stateVariable)[numPaddedPoints] = it->var(concreteLts->stateVariable);
    real(*rs_sl0)[numPaddedPoints] = it->var(concreteLts->rs_sl0);
    real(*rs_a)[numPaddedPoints] = it->var(concreteLts->rs_a);
    real(*initialStressInFaultCS)[numPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    real initialSlipRate = std::sqrt(std::pow(drParameters.rs_initialSlipRate1, 2) +
                                     std::pow(drParameters.rs_initialSlipRate1, 2));

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        DS[ltsFace][pointIndex] = drParameters.isDsOutputOn;
        dynStressTime[ltsFace][pointIndex] = 0.0;
        slipRateStrike[ltsFace][pointIndex] = drParameters.rs_initialSlipRate1;
        slipRateDip[ltsFace][pointIndex] = drParameters.rs_initialSlipRate2;
        // compute initial friction and state
        std::tie(stateVariable[ltsFace][pointIndex], mu[ltsFace][pointIndex]) =
            computeInitialStateAndFriction(initialStressInFaultCS[ltsFace][pointIndex][3],
                                           initialStressInFaultCS[ltsFace][pointIndex][5],
                                           initialStressInFaultCS[ltsFace][pointIndex][0],
                                           rs_a[ltsFace][pointIndex],
                                           drParameters.rs_b,
                                           rs_sl0[ltsFace][pointIndex],
                                           drParameters.rs_sr0,
                                           drParameters.rs_f0,
                                           initialSlipRate);
      }
      averagedSlip[ltsFace] = 0.0;
    }
    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      e_interoperability->copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averagedSlip, dynStressTime, slipRateStrike, slipRateDip, mu);
      e_interoperability->copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVariable);
    }
  }
}

std::pair<real, real>
    RateAndStateInitializer::computeInitialStateAndFriction(real tractionXY,
                                                            real tractionXZ,
                                                            real pressure,
                                                            real rs_a,
                                                            real rs_b,
                                                            real rs_sl0,
                                                            real rs_sr0,
                                                            real rs_f0,
                                                            real initialSlipRate) {
  real absoluteTraction = std::sqrt(std::pow(tractionXY, 2) + std::pow(tractionXZ, 2));
  real tmp = std::abs(absoluteTraction / (rs_a * pressure));
  real stateVariable = rs_a / rs_sr0 *
                       std::exp((rs_a * std::log(std::exp(tmp) - std::exp(-tmp)) - rs_f0 -
                                 rs_a * std::log(initialSlipRate / rs_sr0)) /
                                drParameters.rs_b);
  real tmp2 = initialSlipRate * 0.5 / rs_sr0 *
              std::exp((drParameters.rs_f0 +
                        drParameters.rs_b * std::log(drParameters.rs_f0 * stateVariable / rs_sl0)) /
                       rs_a);
  // asinh(x)=log(x+sqrt(x^2+1))
  real mu = rs_a * std::log(tmp2 + std::exp(tmp2 * tmp2 + 1.0));
  return {stateVariable, mu};
}

void RateAndStateInitializer::addAdditionalParameters(
    std::map<std::string, double*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);
  real(*rs_sl0)[numPaddedPoints] = it->var(concreteLts->rs_sl0);
  real(*rs_a)[numPaddedPoints] = it->var(concreteLts->rs_a);
  parameterToStorageMap.insert({"RS_sl0", (double*)rs_sl0});
  parameterToStorageMap.insert({"RS_a", (double*)rs_a});
}

std::pair<real, real>
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(real tractionXY,
                                                                        real tractionXZ,
                                                                        real pressure,
                                                                        real rs_a,
                                                                        real rs_b,
                                                                        real rs_sl0,
                                                                        real rs_sr0,
                                                                        real rs_f0,
                                                                        real initialSlipRate) {
  real absoluteTraction = std::sqrt(std::pow(tractionXY, 2) + std::pow(tractionXZ, 2));
  real tmp = std::abs(absoluteTraction / (rs_a * pressure));
  real stateVariable =
      rs_a * std::log(2.0 * rs_sr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  real tmp2 = initialSlipRate * 0.5 / rs_sr0 * std::exp(stateVariable / rs_a);
  // asinh(x)=log(x+sqrt(x^2+1))
  real mu = rs_a * std::log(tmp2 + std::sqrt(std::pow(tmp2, 2) + 1.0));
  return {stateVariable, mu};
}

void RateAndStateFastVelocityInitializer::addAdditionalParameters(
    std::map<std::string, double*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);
  real(*rs_srW)[numPaddedPoints] = it->var(concreteLts->rs_srW);
  real(*rs_fW)[numPaddedPoints] = it->var(concreteLts->rs_fW);
  parameterToStorageMap.insert({"RS_srW", (double*)rs_srW});
  parameterToStorageMap.insert({"rs_fW", (double*)rs_fW});
}

void RateAndStateThermalPressurisationInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* e_interoperability) {
  RateAndStateInitializer::initializeFault(dynRup, dynRupTree, e_interoperability);

  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*temperature)[numPaddedPoints] = it->var(concreteLts->temperature);
    real(*pressure)[numPaddedPoints] = it->var(concreteLts->pressure);
    real(*TP_Theta)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_theta);
    real(*TP_sigma)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_sigma);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters.iniTemp;
        pressure[ltsFace][pointIndex] = drParameters.iniPressure;
        for (unsigned iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; ++iTP_grid_nz) {
          TP_Theta[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
          TP_sigma[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
        }
      }
    }
  }
}

void RateAndStateThermalPressurisationInitializer::addAdditionalParameters(
    std::map<std::string, double*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateFastVelocityInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);

  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);

  real(*TP_halfWidthShearZone)[numPaddedPoints] = it->var(concreteLts->TP_halfWidthShearZone);
  real(*alphaHy)[numPaddedPoints] = it->var(concreteLts->alphaHy);
  parameterToStorageMap.insert({"TP_halfWidthShearZone", (double*)TP_halfWidthShearZone});
  parameterToStorageMap.insert({"alphaHy", (double*)alphaHy});
}
} // namespace seissol::dr::initializers