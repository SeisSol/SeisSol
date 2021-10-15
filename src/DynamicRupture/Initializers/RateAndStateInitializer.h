#ifndef SEISSOL_RATEANDSTATEINITIALIZER_H
#define SEISSOL_RATEANDSTATEINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class RateAndStateFL103Initializer;   // rate and state with time and space dependent nucleation
                                      // parameter
class RateAndStateFL103TPInitializer; // Fl103 extended with thermal pressurization
} // namespace seissol::dr::initializers

/*
 * time and space dependent nucleation parameter for FL103
 * plus state Variable and dynamic stress required for rate and state friction laws
 */
class seissol::dr::initializers::RateAndStateFL103Initializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture* dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability& e_interoperability) override;
};

/*
 * initialize all thermal pressure parameters
 */
class seissol::dr::initializers::RateAndStateFL103TPInitializer
    : public seissol::dr::initializers::RateAndStateFL103Initializer {
  public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture* dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability& e_interoperability) override;
};

#endif // SEISSOL_RATEANDSTATEINITIALIZER_H
