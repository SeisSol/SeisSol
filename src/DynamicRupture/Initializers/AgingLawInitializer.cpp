#include "AgingLawInitializer.h"

namespace seissol::dr::initializers {
void AgingLawInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  BaseDRInitializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);
  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*RS_a)[numPaddedPoints] = it->var(concreteLts->RS_a);
    real(*RS_sl0)[numPaddedPoints] = it->var(concreteLts->RS_sl0);
    real(*RS_sr0)[numPaddedPoints] = it->var(concreteLts->RS_sr0);
    real(*stateVariable)[numPaddedPoints] = it->var(concreteLts->stateVariable);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      // get initial values from fortran
      // TODO: RS_a, RS_b, rs_f0, RS_sl0, RS_sr0 could be obtained from paramter file
      // e_interoperability.getDynRupFL_3(ltsFace, meshFace, RS_a, RS_sl0, RS_sr0);
      e_interoperability.getDynRupStateVar(ltsFace, meshFace, stateVariable);

    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}
} // namespace seissol::dr::initializers