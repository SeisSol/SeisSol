#include "LinearSlipWeakeningInitializer.h"

namespace seissol::dr::initializers {

void seissol::dr::initializers::LinearSlipWeakeningFL2Initializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  BaseDRInitializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  seissol::initializers::LTS_LinearSlipWeakeningFL2* ConcreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningFL2*>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*d_c)[numOfPointsPadded] = it->var(ConcreteLts->d_c);   // from faultParameters
    real(*mu_S)[numOfPointsPadded] = it->var(ConcreteLts->mu_S); // from faultParameters
    real(*mu_D)[numOfPointsPadded] = it->var(ConcreteLts->mu_D); // from faultParameters
    bool(*DS)[numOfPointsPadded] = it->var(ConcreteLts->DS);     // from parameter file
    real* averaged_Slip = it->var(ConcreteLts->averaged_Slip);   // = 0
    real(*dynStress_time)[numOfPointsPadded] = it->var(ConcreteLts->dynStress_time); // = 0

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      // initialize padded elements for vectorization
      for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {
        d_c[ltsFace][iBndGP] = 0.0;
        mu_S[ltsFace][iBndGP] = 0.0;
        mu_D[ltsFace][iBndGP] = 0.0;
      }
      for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
        d_c[ltsFace][iBndGP] =
            static_cast<real>(faultParameters["d_c"][(meshFace)*numberOfPoints + iBndGP]);
        mu_S[ltsFace][iBndGP] =
            static_cast<real>(faultParameters["mu_s"][(meshFace)*numberOfPoints + iBndGP]);
        mu_D[ltsFace][iBndGP] =
            static_cast<real>(faultParameters["mu_d"][(meshFace)*numberOfPoints + iBndGP]);
      }
      averaged_Slip[ltsFace] = 0.0;
      for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded;
           ++iBndGP) { // loop includes padded elements
        dynStress_time[ltsFace][iBndGP] = 0.0;
        DS[ltsFace][iBndGP] = m_Params->IsDsOutputOn;
      }
    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}

void LinearSlipWeakeningFL16Initializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  LinearSlipWeakeningFL2Initializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  seissol::initializers::LTS_LinearSlipWeakeningFL16* ConcreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningFL16*>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    real(*forced_rupture_time)[numOfPointsPadded] =
        it->var(ConcreteLts->forced_rupture_time); // from faultParameters
    real* tn = it->var(ConcreteLts->tn);           // = 0

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      // initialize padded elements for vectorization
      for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {
        forced_rupture_time[ltsFace][iBndGP] = 0.0;
      }
      for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
        if (faultParameters["forced_rupture_time"] != NULL) {
          forced_rupture_time[ltsFace][iBndGP] = static_cast<real>(
              faultParameters["forced_rupture_time"][(meshFace)*numberOfPoints + iBndGP]);
        }
      }
      tn[ltsFace] = 0.0;
    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}

void LinearBimaterialFL6Initializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  LinearSlipWeakeningFL2Initializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  seissol::initializers::LTS_LinearBimaterialFL6* ConcreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearBimaterialFL6*>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*strengthData)[numOfPointsPadded] = it->var(ConcreteLts->strengthData);
    real(*mu)[numOfPointsPadded] = it->var(ConcreteLts->mu);
    real(*initialStressInFaultCS)[numOfPointsPadded][6] = it->var(dynRup->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      // unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {
        strengthData[ltsFace][iBndGP] =
            mu[ltsFace][iBndGP] * initialStressInFaultCS[ltsFace][iBndGP][0];
      }
    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}
} // namespace seissol::dr::initializers
