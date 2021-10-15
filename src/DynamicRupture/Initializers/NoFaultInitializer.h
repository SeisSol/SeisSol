#ifndef SEISSOL_NOFAULTINITIALIZER_H
#define SEISSOL_NOFAULTINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class NoFaultInitializer; // NoFaultFL0
}

class seissol::dr::initializers::NoFaultInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture* dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability& e_interoperability) override;
};

#endif // SEISSOL_NOFAULTINITIALIZER_H
