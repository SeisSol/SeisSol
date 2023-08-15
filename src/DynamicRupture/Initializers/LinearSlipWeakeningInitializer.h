#ifndef SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
#define SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {

/**
 * Derived initializer class for the LinearSlipWeakening friction law
 */
template <typename Config>
class LinearSlipWeakeningInitializer : public BaseDRInitializer<Config> {
  public:
  using BaseDRInitializer<Config>::BaseDRInitializer;
  /**
   * Computes initial friction and slip rates
   */
  void initializeFault(seissol::initializers::DynamicRupture<Config> const* const dynRup,
                       seissol::initializers::LTSTree* const dynRupTree) override;

  protected:
  /**
   * Adds the additional parameters mu_s, mu_d, d_c, cohesion and if available forced_rupture_time.
   */
  void addAdditionalParameters(
      std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

/**
 * Derived initializer class for the LinearSlipWeakening friction law with bimaterial regularization
 */
template <typename Config>
class LinearSlipWeakeningBimaterialInitializer : public LinearSlipWeakeningInitializer<Config> {
  public:
  using LinearSlipWeakeningInitializer<Config>::LinearSlipWeakeningInitializer;
  /**
   * Computes initial value for the regularized strength
   */
  void initializeFault(seissol::initializers::DynamicRupture<Config> const* const dynRup,
                       seissol::initializers::LTSTree* const dynRupTree) override;
};

} // namespace seissol::dr::initializers
#endif // SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
