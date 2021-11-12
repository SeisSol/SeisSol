#ifndef SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
#define SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class LinearSlipWeakeningInitializer;
class LinearSlipWeakeningBimaterialInitializer;
class LinearSlipWeakeningForcedRuptureTimeInitializer;
} // namespace seissol::dr::initializers

/**
 * Derived initializer class for the LinearSlipWeakening friction law
 */
class seissol::dr::initializers::LinearSlipWeakeningInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  /**
   * Computes initial friction and slip rates
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;

  protected:
  /**
   * Adds the additional parameters mu_s, mu_d, d_c, cohesion.
   */
  virtual void
      addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

/**
 * Derived initializer class for the LinearSlipWeakening friction law with a forced rupture time
 */
class seissol::dr::initializers::LinearSlipWeakeningForcedRuptureTimeInitializer
    : public seissol::dr::initializers::LinearSlipWeakeningInitializer {
  public:
  using LinearSlipWeakeningInitializer::LinearSlipWeakeningInitializer;
  /**
   * initializes tn to 0
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;

  protected:
  /**
   * Reads the additional parameter forced_ruptre_time
   */
  virtual void
      addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

/**
 * Derived initializer class for the LinearSlipWeakening friction law with bimaterial regularization
 */
class seissol::dr::initializers::LinearSlipWeakeningBimaterialInitializer
    : public seissol::dr::initializers::LinearSlipWeakeningInitializer {
  public:
  using LinearSlipWeakeningInitializer::LinearSlipWeakeningInitializer;
  /**
   * Computes initial value for the regularized strength
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;
};

#endif // SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
