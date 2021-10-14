#ifndef SEISSOL_AGINGLAWINITIALIZER_H
#define SEISSOL_AGINGLAWINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
  class AgingLawInitializer;         //aging law (need revisit)
}

/*
 * should be revisted !
 * parameter could be obtain with m_Param (but initializer in Fortran does currently not suppport this friction law and FL4,7)
 */
class seissol::dr::initializers::AgingLawInitializer : public seissol::dr::initializers::BaseDRInitializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override;
};

#endif //SEISSOL_AGINGLAWINITIALIZER_H
