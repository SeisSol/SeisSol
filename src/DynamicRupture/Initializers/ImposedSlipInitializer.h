#ifndef SEISSOL_IMPOSEDSLIPINITIALIZER_H
#define SEISSOL_IMPOSEDSLIPINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class ImposedSlipRatesFL33Initializer; // imposed slip rates on boundary
}

/*
 * nucleationStressInFaultCS initialized which is used to impose slip rates on the fault surface
 */
class seissol::dr::initializers::ImposedSlipRatesFL33Initializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture* dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability& e_interoperability) override;
};

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
