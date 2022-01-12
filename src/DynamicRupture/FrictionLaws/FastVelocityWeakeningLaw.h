#ifndef SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H
#define SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {

class FastVelocityWeakeningLaw : public RateAndStateBase<FastVelocityWeakeningLaw> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw>::RateAndStateBase;
  real (*srW)[numPaddedPoints];

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
  real updateStateVariable(int pointIndex,
                           unsigned int face,
                           real stateVarReference,
                           real timeIncrement,
                           real localSlipRate);

  /**
   * Computes the friction coefficient from the state variable and slip rate
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp
   * \left(\frac{\Theta}{a}\right)\right). \f] \f$V\f$ is taken from the stored values.
   * @param localStateVariable \f$ \Theta \f$
   * @return \f$ \mu \f$
   */
  real updateMu(unsigned int ltsFace, unsigned int pointIndex, real localStateVariable);
};

class RateAndStateThermalPressurizationLaw : public FastVelocityWeakeningLaw {
  public:
  using FastVelocityWeakeningLaw::FastVelocityWeakeningLaw;

  protected:
  real (*temperature)[numPaddedPoints];
  real (*pressure)[numPaddedPoints];
  real (*TP_Theta)[numPaddedPoints][TP_grid_nz];
  real (*TP_sigma)[numPaddedPoints][TP_grid_nz];
  real (*TP_halfWidthShearZone)[numPaddedPoints];
  real (*alphaHy)[numPaddedPoints];

  real TP_grid[TP_grid_nz];
  real TP_DFinv[TP_grid_nz];

  real faultStrength[numPaddedPoints];
  real Theta_tmp[TP_grid_nz];
  real Sigma_tmp[TP_grid_nz];

  public:
  /*
   * initialize local attributes (used in initializer class respectively)
   */
  void initializeTP(seissol::Interoperability& e_interoperability);

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
  void hookSetInitialP_f(std::array<real, numPaddedPoints>& P_f, unsigned int ltsFace);

  /*
   * compute thermal pressure according to Noda and Lapusta 2010
   * bool saveTmpInTP is used to save final thermal pressure values for theta and sigma
   */
  void hookCalcP_f(std::array<real, numPaddedPoints>& P_f,
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

  real heatSource(real tmp, real alpha, unsigned int iTP_grid_nz, unsigned int timeIndex);
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H
