#ifndef SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
#define SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class LinearSlipWeakeningInitializer; // general initialization for linear slip weakening laws
class LinearSlipWeakeningBimaterialInitializer;        // bimaterial faults
class LinearSlipWeakeningForcedRuptureTimeInitializer; // forced rupture time
} // namespace seissol::dr::initializers

class seissol::dr::initializers::LinearSlipWeakeningInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;

  protected:
  virtual void
      addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

class seissol::dr::initializers::LinearSlipWeakeningForcedRuptureTimeInitializer
    : public seissol::dr::initializers::LinearSlipWeakeningInitializer {
  public:
  using LinearSlipWeakeningInitializer::LinearSlipWeakeningInitializer;
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;

  protected:
  virtual void
      addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

/*
 * strength data is additionally initialized by computing it from friction and initial normal stress
 */
class seissol::dr::initializers::LinearSlipWeakeningBimaterialInitializer
    : public seissol::dr::initializers::LinearSlipWeakeningInitializer {
  public:
  using LinearSlipWeakeningInitializer::LinearSlipWeakeningInitializer;
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;
};

#endif // SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
