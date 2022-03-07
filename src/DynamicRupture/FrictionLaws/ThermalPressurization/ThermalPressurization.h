#ifndef SEISSOL_THERMALPRESSURIZATION_H
#define SEISSOL_THERMALPRESSURIZATION_H

#include <array>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/DynamicRupture.h"
#include "Kernels/precision.hpp"

namespace seissol::dr::friction_law {

template <size_t N>
struct GridPoints : std::array<real, N> {
  constexpr GridPoints() {
    for (size_t i = 0; i < N; ++i) {
      this->at(i) =
          misc::tpMaxWavenumber * std::exp(-misc::tpLogDz * (misc::numberOfTPGridPoints - i - 1));
    }
  }
};

template <size_t N>
struct InverseFourierCoefficients : std::array<real, N> {
  constexpr InverseFourierCoefficients() {
    GridPoints<N> localGridPoints;
    this->front() = std::sqrt(2 / M_PI) * localGridPoints.front() * (1 + 0.5 * misc::tpLogDz);
    this->back() = std::sqrt(2 / M_PI) * localGridPoints.back() * (1 + 0.5 * misc::tpLogDz);

    for (size_t i = 0; i < N; ++i) {
      this->at(i) = std::sqrt(2 / M_PI) * localGridPoints.at(i) * misc::tpLogDz;
    }
  }
};

class ThermalPressurization {
  public:
  ThermalPressurization(DRParameters& drParameters) : drParameters(drParameters){};

  private:
  DRParameters& drParameters;

  protected:
  real (*temperature)[misc::numPaddedPoints];
  real (*pressure)[misc::numPaddedPoints];
  real localPressure[misc::numPaddedPoints];
  real (*tpTheta)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*tpSigma)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*tpHalfWidthShearZone)[misc::numPaddedPoints];
  real (*alphaHy)[misc::numPaddedPoints];

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
} // namespace seissol::dr::friction_law

#endif // SEISSOL_THERMALPRESSURIZATION_H
