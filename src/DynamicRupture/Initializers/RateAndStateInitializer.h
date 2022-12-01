#ifndef SEISSOL_RATEANDSTATEINITIALIZER_H
#define SEISSOL_RATEANDSTATEINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {

/**
 * Derived initializer class for the common part of RateAndState friction laws
 * For the slip and aging law, this initializer is sufficient
 */
class RateAndStateInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;

  /**
   * Computes initial friction and slip rates
   */
  void initializeFault(seissol::initializers::DynamicRupture const* const dynRup,
                       seissol::initializers::LTSTree* const dynRupTree) override;

  protected: /**
              * Adds the additional parameters sl0, rs_a
              */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               seissol::initializers::DynamicRupture const* const dynRup,
                               seissol::initializers::LTSInternalNode::leaf_iterator& it) override;

  struct StateAndFriction {
    double stateVariable;
    double frictionCoefficient;
  };

  /**
   * Computes the initial stateVariable and frictionCoefficient
   * \f[ \mathbf{\tau} = \sqrt{\tau_{XY}^2 + \tau_{XZ}^2} \f]
   * \f[ \psi = \frac{sl_0}{sr_0} \cdot \exp\left(\frac{a \cdot \log\left(2
   * \sinh\left[\left|\frac{\mathbf{\tau}}{a \cdot p}\right|\right]\right) - f_0 - a \cdot
   * \log\left(\frac{sr_{ini}}{sr_0}\right)}{b}\right) \f] \f[ \mu = a \cdot \sinh^{-1}\left(
   * \frac{rs_{ini}}{2 \cdot sr_0} \cdot \exp\left(\frac{f_0 + b * log\left(\frac{sr_0 \cdot
   * \psi}{sl_0}\right)}{a}\right) \right) \f]
   * @param traction1 \f$ \tau_{XY} \f$
   * @param traction2 \f$ \tau_{XZ} \f$
   * @param pressure \f$ p \f$
   * @param rs_a \f$ a \f$
   * @param rs_b \f$ b \f$
   * @param rs_sl0 \f$ sl_0 \f$
   * @param rs_sr0 \f$ sr_0 \f$
   * @param rs_f0 \f$ f_0 \f$
   * @param initialSlipRate \f$ rs_{ini} \f$
   * @return \f$ \left( \psi, \mu \right) \f$
   */
  virtual StateAndFriction computeInitialStateAndFriction(real traction1,
                                                          real traction2,
                                                          real pressure,
                                                          real rsA,
                                                          real rsB,
                                                          real rsSl0,
                                                          real rsSr0,
                                                          real rsF0,
                                                          real initialSlipRate);
};

/**
 * Derived initializer class for FastVelocityWeakening friction laws
 */
class RateAndStateFastVelocityInitializer : public RateAndStateInitializer {
  public:
  using RateAndStateInitializer::RateAndStateInitializer;

  protected:
  /**
   * Adds the additional parameters rs_srW
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               seissol::initializers::DynamicRupture const* const dynRup,
                               seissol::initializers::LTSInternalNode::leaf_iterator& it) override;

  /**
  \f[ \mathbf{\tau} = \sqrt{\tau_{XY}^2 + \tau_{XZ}^2}; \f]
  \f[ \psi = a \cdot \log\left(\frac{2 \cdot sr_0}{rs_{ini}} \cdot
  \sinh\left(\left|\frac{\mathbf{\tau}}{(a \cdot p)}\right|\right)\right); \f] \f[ \mu = a \cdot
  \sinh^{-1}\left(\frac{sr_{ini}}{2 \cdot sr_0} \cdot \exp\left(\frac{\psi}{a}\right)\right); \f]
   * Computes the initial stateVariable and frictionCoefficient
   * @param traction1 \f$ \tau_{XY} \f$
   * @param traction2 \f$ \tau_{XZ} \f$
   * @param pressure \f$ p \f$
   * @param rs_a \f$ a \f$
   * @param rs_b \f$ b \f$
   * @param rs_sl0 \f$ sl_0 \f$
   * @param rs_sr0 \f$ sr_0 \f$
   * @param rs_f0 \f$ f_0 \f$
   * @param initialSlipRate \f$ rs_{ini} \f$
   * @return \f$ \left( \psi, \mu \right) \f$
   */
  StateAndFriction computeInitialStateAndFriction(real traction1,
                                                  real traction2,
                                                  real pressure,
                                                  real rsA,
                                                  real rsB,
                                                  real rsSl0,
                                                  real rsSr0,
                                                  real rsF0,
                                                  real initialSlipRate) override;
};

/**
 * Derived initializer class for FastVelocityWeakening friction law with additional thermal
 * pressurization
 */
class RateAndStateThermalPressurizationInitializer : public RateAndStateFastVelocityInitializer {
  public:
  using RateAndStateFastVelocityInitializer::RateAndStateFastVelocityInitializer;

  /**
   * Intializes temperature and pressure and sets compute grid to 0
   */
  void initializeFault(seissol::initializers::DynamicRupture const* const dynRup,
                       seissol::initializers::LTSTree* const dynRupTree) override;

  protected:
  /**
   * Adds the additional parameters halfWidthShearZone and hydraulicDiffusivity
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               seissol::initializers::DynamicRupture const* const dynRup,
                               seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

} // namespace seissol::dr::initializers
#endif // SEISSOL_RATEANDSTATEINITIALIZER_H
