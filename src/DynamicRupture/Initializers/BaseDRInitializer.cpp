#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
void BaseDRInitializer::setInputParam(dr::DRParameters* DynRupParameter) {
  m_Params = DynRupParameter;
}

void BaseDRInitializer::initializeFrictionMatrices(
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
    real(*iniBulkXX)[numPaddedPoints] = it->var(dynRup->iniBulkXX);   // get from faultParameters
    real(*iniBulkYY)[numPaddedPoints] = it->var(dynRup->iniBulkYY);   // get from faultParameters
    real(*iniBulkZZ)[numPaddedPoints] = it->var(dynRup->iniBulkZZ);   // get from faultParameters
    real(*iniShearXY)[numPaddedPoints] = it->var(dynRup->iniShearXY); // get from faultParameters
    real(*iniShearXZ)[numPaddedPoints] = it->var(dynRup->iniShearXZ); // get from faultParameters
    real(*iniShearYZ)[numPaddedPoints] = it->var(dynRup->iniShearYZ); // get from faultParameters
    real(*initialStressInFaultCS)[numPaddedPoints][6] =
        it->var(dynRup->initialStressInFaultCS); // get from fortran  EQN%InitialStressInFaultCS
    real(*cohesion)[numPaddedPoints] = it->var(dynRup->cohesion); // get from faultParameters
    real(*mu)[numPaddedPoints] = it->var(dynRup->mu);     // get from fortran  EQN%IniMu(:,:)
    real(*slip)[numPaddedPoints] = it->var(dynRup->slip); // = 0
    real(*slipStrike)[numPaddedPoints] = it->var(dynRup->slipStrike);               // = 0
    real(*slipDip)[numPaddedPoints] = it->var(dynRup->slipDip);                     // = 0
    real(*slipRateMagnitude)[numPaddedPoints] = it->var(dynRup->slipRateMagnitude); // = 0
    real(*slipRateStrike)[numPaddedPoints] =
        it->var(dynRup->slipRateStrike); // get from fortran  EQN%IniSlipRate1
    real(*slipRateDip)[numPaddedPoints] =
        it->var(dynRup->slipRateDip); // get from fortran  EQN%IniSlipRate2
    real(*ruptureTime)[numPaddedPoints] = it->var(dynRup->ruptureTime);   // = 0
    bool(*ruptureFront)[numPaddedPoints] = it->var(dynRup->ruptureFront); // get from fortran
    real(*peakSlipRate)[numPaddedPoints] = it->var(dynRup->peakSlipRate); // = 0
    real(*tractionXY)[numPaddedPoints] = it->var(dynRup->tractionXY);     // = 0
    real(*tractionXZ)[numPaddedPoints] = it->var(dynRup->tractionXZ);     // = 0

    dynRup->IsFaultParameterizedByTraction = false;

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0];
           ++pointIndex) { // loop includes padded elements
        slip[ltsFace][pointIndex] = 0.0;
        slipStrike[ltsFace][pointIndex] = 0.0;
        slipDip[ltsFace][pointIndex] = 0.0;
        slipRateMagnitude[ltsFace][pointIndex] = 0.0;
        ruptureTime[ltsFace][pointIndex] = 0.0;
        peakSlipRate[ltsFace][pointIndex] = 0.0;
        tractionXY[ltsFace][pointIndex] = 0.0;
        tractionXZ[ltsFace][pointIndex] = 0.0;
      }
      // get initial values from fortran
      for (unsigned pointIndex = 0; pointIndex < numberOfPoints; ++pointIndex) {
        if (faultParameters["T_n"] != NULL) {
          iniBulkXX[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["T_n"][(meshFace)*numberOfPoints + pointIndex]);
          dynRup->IsFaultParameterizedByTraction = true;
        } else if (faultParameters["s_xx"] != NULL)
          iniBulkXX[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["s_xx"][(meshFace)*numberOfPoints + pointIndex]);
        else
          iniBulkXX[ltsFace][pointIndex] = 0.0;

        if (faultParameters["s_yy"] != NULL)
          iniBulkYY[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["s_yy"][(meshFace)*numberOfPoints + pointIndex]);
        else
          iniBulkYY[ltsFace][pointIndex] = 0.0;

        if (faultParameters["s_zz"] != NULL)
          iniBulkZZ[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["s_zz"][(meshFace)*numberOfPoints + pointIndex]);
        else
          iniBulkZZ[ltsFace][pointIndex] = 0.0;

        if (faultParameters["T_s"] != NULL)
          iniShearXY[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["T_s"][(meshFace)*numberOfPoints + pointIndex]);
        else if (faultParameters["s_xy"] != NULL)
          iniShearXY[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["s_xy"][(meshFace)*numberOfPoints + pointIndex]);
        else
          iniShearXY[ltsFace][pointIndex] = 0.0;

        if (faultParameters["T_d"] != NULL)
          iniShearXZ[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["T_d"][(meshFace)*numberOfPoints + pointIndex]);
        else if (faultParameters["s_xz"] != NULL)
          iniShearXZ[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["s_xz"][(meshFace)*numberOfPoints + pointIndex]);
        else
          iniShearXZ[ltsFace][pointIndex] = 0.0;

        if (faultParameters["s_yz"] != NULL)
          iniShearYZ[ltsFace][pointIndex] =
              static_cast<real>(faultParameters["s_yz"][(meshFace)*numberOfPoints + pointIndex]);
        else
          iniShearYZ[ltsFace][pointIndex] = 0.0;

        if (faultParameters["cohesion"] != NULL) {
          cohesion[ltsFace][pointIndex] = static_cast<real>(
              faultParameters["cohesion"][(meshFace)*numberOfPoints + pointIndex]);
        } else {
          cohesion[ltsFace][pointIndex] = 0.0;
        }
      }
      // initialize padded elements for vectorization
      for (unsigned pointIndex = numberOfPoints; pointIndex < numPaddedPoints; ++pointIndex) {
        cohesion[ltsFace][pointIndex] = 0.0;
      }
      e_interoperability.getDynRupParameters(
          ltsFace, meshFace, initialStressInFaultCS, mu, slipRateStrike, slipRateDip, ruptureFront);

    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}
} // namespace seissol::dr::initializers