#ifndef SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
#define SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
  class LinearSlipWeakeningFL2Initializer;  //general initialization for linear slip laws
  class LinearBimaterialFL6Initializer;     //for bimaterial faults
  class LinearSlipWeakeningFL16Initializer; //FL2 extended by forced rupture time
}

class seissol::dr::initializers::LinearSlipWeakeningFL2Initializer : public seissol::dr::initializers::BaseDRInitializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string,
                                              double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override;
};

class seissol::dr::initializers::LinearSlipWeakeningFL16Initializer : public seissol::dr::initializers::LinearSlipWeakeningFL2Initializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override;
};

/*
 * strength data is additionally initialized by computing it from friction and initial normal stress
 */
class seissol::dr::initializers::LinearBimaterialFL6Initializer : public seissol::dr::initializers::LinearSlipWeakeningFL2Initializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string,
                                              double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability);
};

#endif //SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
