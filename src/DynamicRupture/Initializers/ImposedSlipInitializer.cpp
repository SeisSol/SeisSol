#include "ImposedSlipInitializer.h"

/*
 * nucleationStressInFaultCS initialized which is used to impose slip rates on the fault surface
 */
void seissol::dr::initializers::ImposedSlipRatesFL33Initializer::initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                                                                       seissol::initializers::LTSTree* dynRupTree,
                                                                                       seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                                                                       std::unordered_map<std::string, double*> faultParameters,
                                                                                       unsigned* ltsFaceToMeshFace,
                                                                                       seissol::Interoperability &e_interoperability) {
  BaseDRInitializer::initializeFrictionMatrices(dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  seissol::initializers::LTS_ImposedSlipRatesFL33 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesFL33 *>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {

    real  (*nucleationStressInFaultCS)[numOfPointsPadded][6]  = it->var(ConcreteLts->nucleationStressInFaultCS); //get from fortran
    real *averaged_Slip                                       = it->var(ConcreteLts->averaged_Slip);      // = 0

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      e_interoperability.getDynRupNucStress(ltsFace, meshFace, nucleationStressInFaultCS);
      averaged_Slip[ltsFace]= 0.0;

    }//lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  }//leaf_iterator loop
}