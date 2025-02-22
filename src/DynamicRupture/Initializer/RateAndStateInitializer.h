// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_RATEANDSTATEINITIALIZER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_RATEANDSTATEINITIALIZER_H_

#include "BaseDRInitializer.h"

namespace seissol::dr::initializer {

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
  void initializeFault(const seissol::initializer::DynamicRupture* dynRup,
                       seissol::initializer::LTSTree* dynRupTree) override;

  protected: /**
              * Adds the additional parameters sl0, rs_a
              */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;

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
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;

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
class ThermalPressurizationInitializer {
  public:
  explicit ThermalPressurizationInitializer(
      const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters);
  /**
   * Intializes temperature and pressure and sets compute grid to 0
   */
  void initializeFault(const seissol::initializer::DynamicRupture* dynRup,
                       seissol::initializer::LTSTree* dynRupTree);

  protected:
  /**
   * Adds the additional parameters halfWidthShearZone and hydraulicDiffusivity
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer);

  private:
  std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters;
};

/**
 * Derived initializer class for RateAndState (slow, severe) friction law with additional thermal
 * pressurization
 */
class RateAndStateThermalPressurizationInitializer : public RateAndStateInitializer,
                                                     public ThermalPressurizationInitializer {
  public:
  RateAndStateThermalPressurizationInitializer(
      const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
      SeisSol& instance);

  /**
   * Intializes temperature and pressure and sets compute grid to 0
   */
  void initializeFault(const seissol::initializer::DynamicRupture* dynRup,
                       seissol::initializer::LTSTree* dynRupTree) override;

  protected:
  /**
   * Adds the additional parameters halfWidthShearZone and hydraulicDiffusivity
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;
};

/**
 * Derived initializer class for FastVelocityWeakening friction law with additional thermal
 * pressurization
 */
class RateAndStateFastVelocityThermalPressurizationInitializer
    : public RateAndStateFastVelocityInitializer,
      public ThermalPressurizationInitializer {
  public:
  RateAndStateFastVelocityThermalPressurizationInitializer(
      const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
      SeisSol& instance);

  /**
   * Intializes temperature and pressure and sets compute grid to 0
   */
  void initializeFault(const seissol::initializer::DynamicRupture* dynRup,
                       seissol::initializer::LTSTree* dynRupTree) override;

  protected:
  /**
   * Adds the additional parameters halfWidthShearZone and hydraulicDiffusivity
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;
};

} // namespace seissol::dr::initializer

#endif // SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_RATEANDSTATEINITIALIZER_H_
