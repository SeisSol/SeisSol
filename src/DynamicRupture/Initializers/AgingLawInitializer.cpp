#include "AgingLawInitializer.h"

void seissol::dr::initializers::AgingLawInitializer::initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                                                             seissol::initializers::LTSTree* dynRupTree,
                                                                             seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                                                             std::unordered_map<std::string, double*> faultParameters,
                                                                             unsigned* ltsFaceToMeshFace,
                                                                             seissol::Interoperability &e_interoperability) {
  BaseDRInitializer::initializeFrictionMatrices(dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  seissol::initializers::LTS_RateAndStateFL3 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_RateAndStateFL3 *>(dynRup);
  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
    real *RS_f0                                               = it->var(ConcreteLts->RS_f0);
    real *RS_a                                                = it->var(ConcreteLts->RS_a);
    real *RS_b                                                = it->var(ConcreteLts->RS_b);
    real *RS_sl0                                              = it->var(ConcreteLts->RS_sl0);
    real *RS_sr0                                              = it->var(ConcreteLts->RS_sr0);
    real (*stateVar)[numOfPointsPadded]                       = it->var(ConcreteLts->stateVar);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

//get initial values from fortran
//TODO: RS_a, RS_b, rs_f0, RS_sl0, RS_sr0 could be obtained from paramter file
      e_interoperability.getDynRupFL_3(ltsFace, meshFace, RS_f0, RS_a, RS_b, RS_sl0, RS_sr0);
      e_interoperability.getDynRupStateVar(ltsFace, meshFace, stateVar);

    }//lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  }//leaf_iterator loop
}