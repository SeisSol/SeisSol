#ifndef SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H
#define SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {

class FastVelocityWeakeningLaw : public RateAndStateBase<FastVelocityWeakeningLaw> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw>::RateAndStateBase;
  real (*srW)[misc::numPaddedPoints];

  /**
   * Copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  /**
   * Integrates the state variable ODE in time
   * \f[\frac{\partial \Theta}{\partial t} = - \frac{V}{L}\left(\Theta - \Theta_{ss}(V) \right)\f]
   * with steady state variable \f$\Theta_{ss}\f$.
   * Assume \f$V\f$ is constant through the time interval, then the analytic solution is:
   * \f[ \Theta(t) = \Theta_0 \exp\left( -\frac{V}{L} t \right) + \Theta_{ss} \left( 1 - \exp\left(
   * - \frac{V}{L} t\right) \right).\f]
   * @param stateVarReference \f$ \Theta_0 \f$
   * @param timeIncrement \f$ t \f$
   * @param localSlipRate \f$ V \f$
   * @return \f$ \Theta(t) \f$
   */
  real updateStateVariable(unsigned int pointIndex,
                           unsigned int face,
                           real stateVarReference,
                           real timeIncrement,
                           real localSlipRate);

  /**
   * Computes the friction coefficient from the state variable and slip rate
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp
   * \left(\frac{\Theta}{a}\right)\right). \f]
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Theta \f$
   * @return \f$ \mu \f$
   */
  real updateMu(unsigned int ltsFace,
                unsigned int pointIndex,
                real localSlipRateMagnitude,
                real localStateVariable);

  /**
   * Computes the derivative of the friction coefficient with respect to the slip rate.
   * \f[\frac{\partial}{\partial V}\mu = \frac{aC}{\sqrt{ (VC)^2 + 1} \text{ with } C =
   * \frac{1}{2V_0} \cdot \exp \left(\frac{\Theta}{a}\right)\right).\f]
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Theta \f$
   * @return \f$ \mu \f$
   */
  real updateMuDerivative(unsigned int ltsFace,
                          unsigned int pointIndex,
                          real localSlipRateMagnitude,
                          real localStateVariable);

  std::array<real, misc::numPaddedPoints>
      resampleStateVar(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                       unsigned int ltsFace);
};

class RateAndStateThermalPressurizationLaw : public FastVelocityWeakeningLaw {
  public:
  using FastVelocityWeakeningLaw::FastVelocityWeakeningLaw;

  protected:
  real (*temperature)[misc::numPaddedPoints];
  real (*pressure)[misc::numPaddedPoints];
  real (*tpTheta)[misc::numPaddedPoints][numberOfTPGridPoints];
  real (*tpSigma)[misc::numPaddedPoints][numberOfTPGridPoints];
  real (*tpHalfWidthShearZone)[misc::numPaddedPoints];
  real (*alphaHy)[misc::numPaddedPoints];

  real tpGrid[numberOfTPGridPoints];
  real tpDFinv[numberOfTPGridPoints];

  real faultStrength[misc::numPaddedPoints];
  real thetaTmp[numberOfTPGridPoints];
  real sigmaTmp[numberOfTPGridPoints];

  public:
  /*
   * initialize local attributes (used in initializer class respectively)
   */
  void initializeTP(seissol::Interoperability& eInteroperability);

  /*
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  protected:
  /*
   * set initial value of thermal pressure
   */
  void setInitialFluidPressureHook(std::array<real, misc::numPaddedPoints>& fluidPressure,
                                   unsigned int ltsFace);

  /*
   * compute thermal pressure according to Noda and Lapusta 2010
   * bool saveTmpInTP is used to save final thermal pressure values for theta and sigma
   */
  void calcFluidPressureHook(std::array<real, misc::numPaddedPoints>& fluidPressure,
                             FaultStresses& faultStresses,
                             bool saveTmpInTP,
                             unsigned int timeIndex,
                             unsigned int ltsFace);

  /*
   * compute thermal pressure according to Noda and Lapusta 2010
   */
  void updateTemperatureAndPressure(unsigned int pointIndex,
                                    unsigned int timeIndex,
                                    unsigned int ltsFace);

  real heatSource(real tmp, real alpha, unsigned int tpGridPointIndex, unsigned int timeIndex);
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H
