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

    real(*nucleationStressInFaultCS)[numOfPointsPadded][6] =
        it->var(ConcreteLts->nucleationStressInFaultCS); // get from fortran
    real(*RS_sl0_array)[numOfPointsPadded] =
        it->var(ConcreteLts->RS_sl0_array); // get from faultParameters
    real(*RS_a_array)[numOfPointsPadded] =
        it->var(ConcreteLts->RS_a_array); // get from faultParameters
    real(*RS_srW_array)[numOfPointsPadded] =
        it->var(ConcreteLts->RS_srW_array);                    // get from faultParameters
    bool(*DS)[numOfPointsPadded] = it->var(ConcreteLts->DS);   // par file
    real* averaged_Slip = it->var(ConcreteLts->averaged_Slip); // = 0
    real(*stateVar)[numOfPointsPadded] =
        it->var(ConcreteLts->stateVar); // get from Fortran = EQN%IniStateVar
    real(*dynStress_time)[numOfPointsPadded] = it->var(ConcreteLts->dynStress_time); // = 0

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      e_interoperability.getDynRupStateVar(ltsFace, meshFace, stateVar);
      e_interoperability.getDynRupNucStress(ltsFace, meshFace, nucleationStressInFaultCS);

      for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded;
           ++iBndGP) { // loop includes padded elements
        dynStress_time[ltsFace][iBndGP] = 0.0;
        DS[ltsFace][iBndGP] = m_Params->IsDsOutputOn;
      }
      averaged_Slip[ltsFace] = 0.0;

      for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
        RS_a_array[ltsFace][iBndGP] =
            static_cast<real>(faultParameters["rs_a"][(meshFace)*numberOfPoints + iBndGP]);
        RS_srW_array[ltsFace][iBndGP] =
            static_cast<real>(faultParameters["rs_srW"][(meshFace)*numberOfPoints + iBndGP]);
        RS_sl0_array[ltsFace][iBndGP] =
            static_cast<real>(faultParameters["RS_sl0"][(meshFace)*numberOfPoints + iBndGP]);
      }
      // initialize padded elements for vectorization
      for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
        RS_a_array[ltsFace][iBndGP] = 0.0;
        RS_srW_array[ltsFace][iBndGP] = 0.0;
        RS_sl0_array[ltsFace][iBndGP] = 0.0;
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

    real(*temperature)[numOfPointsPadded] = it->var(ConcreteLts->temperature);
    real(*pressure)[numOfPointsPadded] = it->var(ConcreteLts->pressure);

    real(*TP_Theta)[numOfPointsPadded][TP_grid_nz] = it->var(ConcreteLts->TP_theta);
    real(*TP_sigma)[numOfPointsPadded][TP_grid_nz] = it->var(ConcreteLts->TP_sigma);

    real(*TP_half_width_shear_zone)[numOfPointsPadded] =
        it->var(ConcreteLts->TP_half_width_shear_zone);
    real(*alpha_hy)[numOfPointsPadded] = it->var(ConcreteLts->alpha_hy);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {
        temperature[ltsFace][iBndGP] = m_Params->IniTemp;
        pressure[ltsFace][iBndGP] = m_Params->IniPressure;
        TP_half_width_shear_zone[ltsFace][iBndGP] = static_cast<real>(
            faultParameters["TP_half_width_shear_zone"][(meshFace)*numberOfPoints + iBndGP]);
        alpha_hy[ltsFace][iBndGP] =
            static_cast<real>(faultParameters["alpha_hy"][(meshFace)*numberOfPoints + iBndGP]);
        for (unsigned iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; ++iTP_grid_nz) {
          TP_Theta[ltsFace][iBndGP][iTP_grid_nz] = 0.0;
          TP_sigma[ltsFace][iBndGP][iTP_grid_nz] = 0.0;
        }
      }
    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}
} // namespace seissol::dr::initializers