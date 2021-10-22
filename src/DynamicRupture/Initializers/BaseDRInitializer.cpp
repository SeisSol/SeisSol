#include "BaseDRInitializer.h"

void seissol::dr::initializers::BaseDRInitializer::setInputParam(
    dr::DRParameters* DynRupParameter) {
  m_Params = DynRupParameter;
};

void seissol::dr::initializers::BaseDRInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*iniBulkXX)[numOfPointsPadded] = it->var(dynRup->iniBulkXX);   // get from faultParameters
    real(*iniBulkYY)[numOfPointsPadded] = it->var(dynRup->iniBulkYY);   // get from faultParameters
    real(*iniBulkZZ)[numOfPointsPadded] = it->var(dynRup->iniBulkZZ);   // get from faultParameters
    real(*iniShearXY)[numOfPointsPadded] = it->var(dynRup->iniShearXY); // get from faultParameters
    real(*iniShearXZ)[numOfPointsPadded] = it->var(dynRup->iniShearXZ); // get from faultParameters
    real(*iniShearYZ)[numOfPointsPadded] = it->var(dynRup->iniShearYZ); // get from faultParameters
    real(*initialStressInFaultCS)[numOfPointsPadded][6] =
        it->var(dynRup->initialStressInFaultCS); // get from fortran  EQN%InitialStressInFaultCS
    real(*cohesion)[numOfPointsPadded] = it->var(dynRup->cohesion); // get from faultParameters
    real(*mu)[numOfPointsPadded] = it->var(dynRup->mu);     // get from fortran  EQN%IniMu(:,:)
    real(*slip)[numOfPointsPadded] = it->var(dynRup->slip); // = 0
    real(*slipStrike)[numOfPointsPadded] = it->var(dynRup->slipStrike);               // = 0
    real(*slipDip)[numOfPointsPadded] = it->var(dynRup->slipDip);                     // = 0
    real(*slipRateMagnitude)[numOfPointsPadded] = it->var(dynRup->slipRateMagnitude); // = 0
    real(*slipRateStrike)[numOfPointsPadded] =
        it->var(dynRup->slipRateStrike); // get from fortran  EQN%IniSlipRate1
    real(*slipRateDip)[numOfPointsPadded] =
        it->var(dynRup->slipRateDip); // get from fortran  EQN%IniSlipRate2
    real(*rupture_time)[numOfPointsPadded] = it->var(dynRup->rupture_time); // = 0
    bool(*RF)[numOfPointsPadded] = it->var(dynRup->RF);                     // get from fortran
    real(*peakSR)[numOfPointsPadded] = it->var(dynRup->peakSR);             // = 0
    real(*tractionXY)[numOfPointsPadded] = it->var(dynRup->tractionXY);     // = 0
    real(*tractionXZ)[numOfPointsPadded] = it->var(dynRup->tractionXZ);     // = 0

    dynRup->IsFaultParameterizedByTraction = false;

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      for (unsigned iBndGP = 0; iBndGP < init::QInterpolated::Stop[0];
           ++iBndGP) { // loop includes padded elements
        slip[ltsFace][iBndGP] = 0.0;
        slipStrike[ltsFace][iBndGP] = 0.0;
        slipDip[ltsFace][iBndGP] = 0.0;
        slipRateMagnitude[ltsFace][iBndGP] = 0.0;
        rupture_time[ltsFace][iBndGP] = 0.0;
        peakSR[ltsFace][iBndGP] = 0.0;
        tractionXY[ltsFace][iBndGP] = 0.0;
        tractionXZ[ltsFace][iBndGP] = 0.0;
      }
      // get initial values from fortran
      for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
        if (faultParameters["T_n"] != NULL) {
          iniBulkXX[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["T_n"][(meshFace)*numberOfPoints + iBndGP]);
          dynRup->IsFaultParameterizedByTraction = true;
        } else if (faultParameters["s_xx"] != NULL)
          iniBulkXX[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["s_xx"][(meshFace)*numberOfPoints + iBndGP]);
        else
          iniBulkXX[ltsFace][iBndGP] = 0.0;

        if (faultParameters["s_yy"] != NULL)
          iniBulkYY[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["s_yy"][(meshFace)*numberOfPoints + iBndGP]);
        else
          iniBulkYY[ltsFace][iBndGP] = 0.0;

        if (faultParameters["s_zz"] != NULL)
          iniBulkZZ[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["s_zz"][(meshFace)*numberOfPoints + iBndGP]);
        else
          iniBulkZZ[ltsFace][iBndGP] = 0.0;

        if (faultParameters["T_s"] != NULL)
          iniShearXY[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["T_s"][(meshFace)*numberOfPoints + iBndGP]);
        else if (faultParameters["s_xy"] != NULL)
          iniShearXY[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["s_xy"][(meshFace)*numberOfPoints + iBndGP]);
        else
          iniShearXY[ltsFace][iBndGP] = 0.0;

        if (faultParameters["T_d"] != NULL)
          iniShearXZ[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["T_d"][(meshFace)*numberOfPoints + iBndGP]);
        else if (faultParameters["s_xz"] != NULL)
          iniShearXZ[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["s_xz"][(meshFace)*numberOfPoints + iBndGP]);
        else
          iniShearXZ[ltsFace][iBndGP] = 0.0;

        if (faultParameters["s_yz"] != NULL)
          iniShearYZ[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["s_yz"][(meshFace)*numberOfPoints + iBndGP]);
        else
          iniShearYZ[ltsFace][iBndGP] = 0.0;

        if (faultParameters["cohesion"] != NULL) {
          cohesion[ltsFace][iBndGP] =
              static_cast<real>(faultParameters["cohesion"][(meshFace)*numberOfPoints + iBndGP]);
        } else {
          cohesion[ltsFace][iBndGP] = 0.0;
        }
      }
      // initialize padded elements for vectorization
      for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
        cohesion[ltsFace][iBndGP] = 0.0;
      }
      e_interoperability.getDynRupParameters(
          ltsFace, meshFace, initialStressInFaultCS, mu, slipRateStrike, slipRateDip, RF);

    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}
