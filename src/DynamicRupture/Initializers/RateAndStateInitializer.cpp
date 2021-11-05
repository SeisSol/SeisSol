#include "RateAndStateInitializer.h"

namespace seissol::dr::initializers {
void RateAndStateFL103Initializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  BaseDRInitializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    real(*nucleationStressInFaultCS)[numPaddedPoints][6] =
        it->var(concreteLts->nucleationStressInFaultCS);           // get from fortran
    real(*RS_sl0)[numPaddedPoints] = it->var(concreteLts->RS_sl0); // get from faultParameters
    real(*RS_a)[numPaddedPoints] = it->var(concreteLts->RS_a);     // get from faultParameters
    real(*RS_srW)[numPaddedPoints] = it->var(concreteLts->RS_srW); // get from faultParameters
    bool(*DS)[numPaddedPoints] = it->var(concreteLts->DS);         // par file
    real* averaged_Slip = it->var(concreteLts->averagedSlip);      // = 0
    real(*stateVar)[numPaddedPoints] =
        it->var(concreteLts->stateVariable); // get from Fortran = EQN%IniStateVar
    real(*dynStressTime)[numPaddedPoints] = it->var(concreteLts->dynStressTime); // = 0

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      e_interoperability.getDynRupStateVar(ltsFace, meshFace, stateVar);
      e_interoperability.getDynRupNucStress(ltsFace, meshFace, nucleationStressInFaultCS);

      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints;
           ++pointIndex) { // loop includes padded elements
        dynStressTime[ltsFace][pointIndex] = 0.0;
        DS[ltsFace][pointIndex] = drParameters.isDsOutputOn;
      }
      averaged_Slip[ltsFace] = 0.0;

      for (unsigned pointIndex = 0; pointIndex < numberOfPoints; ++pointIndex) {
        RS_a[ltsFace][pointIndex] =
            static_cast<real>(faultParameters["rs_a"][(meshFace)*numberOfPoints + pointIndex]);
        RS_srW[ltsFace][pointIndex] =
            static_cast<real>(faultParameters["rs_srW"][(meshFace)*numberOfPoints + pointIndex]);
        RS_sl0[ltsFace][pointIndex] =
            static_cast<real>(faultParameters["RS_sl0"][(meshFace)*numberOfPoints + pointIndex]);
      }
      // initialize padded elements for vectorization
      for (unsigned pointIndex = numberOfPoints; pointIndex < numPaddedPoints; ++pointIndex) {
        RS_a[ltsFace][pointIndex] = 0.0;
        RS_srW[ltsFace][pointIndex] = 0.0;
        RS_sl0[ltsFace][pointIndex] = 0.0;
      }

    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}

void RateAndStateFL103TPInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  // BaseDrInitializer::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters,
  // ltsFaceToMeshFace, e_interoperability);
  RateAndStateFL103Initializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);

  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);
  seissol::dr::friction_law::RateAndStateThermalFL103* SolverFL103 =
      dynamic_cast<seissol::dr::friction_law::RateAndStateThermalFL103*>(FrictionLaw);
  SolverFL103->initializeTP(e_interoperability);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    real(*temperature)[numPaddedPoints] = it->var(concreteLts->temperature);
    real(*pressure)[numPaddedPoints] = it->var(concreteLts->pressure);

    real(*TP_Theta)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_theta);
    real(*TP_sigma)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_sigma);

    real(*TP_half_width_shear_zone)[numPaddedPoints] =
        it->var(concreteLts->TP_half_width_shear_zone);
    real(*alpha_hy)[numPaddedPoints] = it->var(concreteLts->alpha_hy);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters.iniTemp;
        pressure[ltsFace][pointIndex] = drParameters.iniPressure;
        TP_half_width_shear_zone[ltsFace][pointIndex] = static_cast<real>(
            faultParameters["TP_half_width_shear_zone"][(meshFace)*numberOfPoints + pointIndex]);
        alpha_hy[ltsFace][pointIndex] =
            static_cast<real>(faultParameters["alpha_hy"][(meshFace)*numberOfPoints + pointIndex]);
        for (unsigned iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; ++iTP_grid_nz) {
          TP_Theta[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
          TP_sigma[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
        }
      }
    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}
} // namespace seissol::dr::initializers