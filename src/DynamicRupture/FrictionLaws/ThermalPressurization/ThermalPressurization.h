#ifndef SEISSOL_THERMALPRESSURIZATION_H
#define SEISSOL_THERMALPRESSURIZATION_H

#include "Kernels/precision.hpp"
#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "DynamicRupture/Parameters.h"

namespace seissol::dr::friction_law {
class ThermalPressurization {
  public:
  ThermalPressurization() : drParameters(nullptr){};
  ThermalPressurization(std::shared_ptr<DRParameters> drParameters)
      : drParameters(std::move(drParameters)){};

  private:
  std::shared_ptr<DRParameters> drParameters;

  protected:
  real (*temperature)[misc::numPaddedPoints];
  real (*pressure)[misc::numPaddedPoints];
  real localPressure[misc::numPaddedPoints];
  real (*tpTheta)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*tpSigma)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*tpHalfWidthShearZone)[misc::numPaddedPoints];
  real (*alphaHy)[misc::numPaddedPoints];

  real tpGrid[misc::numberOfTPGridPoints];
  real tpDFinv[misc::numberOfTPGridPoints];

  real faultStrength[misc::numPaddedPoints];
  real thetaTmp[misc::numberOfTPGridPoints];
  real sigmaTmp[misc::numberOfTPGridPoints];

  public:
  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  /**
   * set initial value of thermal pressure
   */
  void setInitialFluidPressure(unsigned int ltsFace);

  /**
   * compute thermal pressure according to Noda and Lapusta 2010
   * bool saveTmpInTP is used to save final thermal pressure values for theta and sigma
   */
  void calcFluidPressure(FaultStresses const& faultStresses,
                         real (*initialStressInFaultCS)[misc::numPaddedPoints][6],
                         real (*mu)[misc::numPaddedPoints],
                         real (*slipRateMagnitude)[misc::numPaddedPoints],
                         real deltaT,
                         bool saveTmpInTP,
                         unsigned int timeIndex,
                         unsigned int ltsFace);

  /**
   * compute thermal pressure according to Noda and Lapusta 2010
   */
  void updateTemperatureAndPressure(real slipRateMagnitude,
                                    real deltaT,
                                    unsigned int pointIndex,
                                    unsigned int timeIndex,
                                    unsigned int ltsFace);

  real heatSource(
      real tmp, real alpha, real deltaT, unsigned int tpGridPointIndex, unsigned int timeIndex);

  real fluidPressure(unsigned int pointIndex) const { return localPressure[pointIndex]; }
};
}; // namespace seissol::dr::friction_law

#endif // SEISSOL_THERMALPRESSURIZATION_H
