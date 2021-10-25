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
  seissol::initializers::LTS_RateAndStateFL103* ConcreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFL103*>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    real(*nucleationStressInFaultCS)[numPaddedPoints][6] =
        it->var(ConcreteLts->nucleationStressInFaultCS); // get from fortran
    real(*RS_sl0_array)[numPaddedPoints] =
        it->var(ConcreteLts->RS_sl0_array); // get from faultParameters
    real(*RS_a_array)[numPaddedPoints] =
        it->var(ConcreteLts->RS_a_array); // get from faultParameters
    real(*RS_srW_array)[numPaddedPoints] =
        it->var(ConcreteLts->RS_srW_array);                    // get from faultParameters
    bool(*DS)[numPaddedPoints] = it->var(ConcreteLts->DS);     // par file
    real* averaged_Slip = it->var(ConcreteLts->averaged_Slip); // = 0
    real(*stateVar)[numPaddedPoints] =
        it->var(ConcreteLts->stateVar); // get from Fortran = EQN%IniStateVar
    real(*dynStress_time)[numPaddedPoints] = it->var(ConcreteLts->dynStress_time); // = 0

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      e_interoperability.getDynRupStateVar(ltsFace, meshFace, stateVar);
      e_interoperability.getDynRupNucStress(ltsFace, meshFace, nucleationStressInFaultCS);

      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints;
           ++pointIndex) { // loop includes padded elements
        dynStress_time[ltsFace][pointIndex] = 0.0;
        DS[ltsFace][pointIndex] = m_Params->IsDsOutputOn;
      }
      averaged_Slip[ltsFace] = 0.0;

      for (unsigned pointIndex = 0; pointIndex < numberOfPoints; ++pointIndex) {
        RS_a_array[ltsFace][pointIndex] =
            static_cast<real>(faultParameters["rs_a"][(meshFace)*numberOfPoints + pointIndex]);
        RS_srW_array[ltsFace][pointIndex] =
            static_cast<real>(faultParameters["rs_srW"][(meshFace)*numberOfPoints + pointIndex]);
        RS_sl0_array[ltsFace][pointIndex] =
            static_cast<real>(faultParameters["RS_sl0"][(meshFace)*numberOfPoints + pointIndex]);
      }
      // initialize padded elements for vectorization
      for (unsigned pointIndex = numberOfPoints; pointIndex < numPaddedPoints; ++pointIndex) {
        RS_a_array[ltsFace][pointIndex] = 0.0;
        RS_srW_array[ltsFace][pointIndex] = 0.0;
        RS_sl0_array[ltsFace][pointIndex] = 0.0;
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

  seissol::initializers::LTS_RateAndStateFL103TP* ConcreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFL103TP*>(dynRup);
  seissol::dr::friction_law::RateAndStateThermalFL103* SolverFL103 =
      dynamic_cast<seissol::dr::friction_law::RateAndStateThermalFL103*>(FrictionLaw);
  SolverFL103->initializeTP(e_interoperability);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    real(*temperature)[numPaddedPoints] = it->var(ConcreteLts->temperature);
    real(*pressure)[numPaddedPoints] = it->var(ConcreteLts->pressure);

    real(*TP_Theta)[numPaddedPoints][TP_grid_nz] = it->var(ConcreteLts->TP_theta);
    real(*TP_sigma)[numPaddedPoints][TP_grid_nz] = it->var(ConcreteLts->TP_sigma);

    real(*TP_half_width_shear_zone)[numPaddedPoints] =
        it->var(ConcreteLts->TP_half_width_shear_zone);
    real(*alpha_hy)[numPaddedPoints] = it->var(ConcreteLts->alpha_hy);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = m_Params->IniTemp;
        pressure[ltsFace][pointIndex] = m_Params->IniPressure;
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