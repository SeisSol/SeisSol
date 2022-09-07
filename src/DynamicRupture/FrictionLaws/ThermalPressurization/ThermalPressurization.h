#ifndef SEISSOL_THERMALPRESSURIZATION_H
#define SEISSOL_THERMALPRESSURIZATION_H

#include <array>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/DynamicRupture.h"
#include "Kernels/precision.hpp"

namespace seissol::dr::friction_law {

/**
 * Logarithmic gridpoints as defined in Noda&Lapusta (14). These are the \f$\hat{l}\f$ for
 * ThermalPressurization.
 */
template <size_t N>
class GridPoints {
  public:
  GridPoints() {
    for (size_t i = 0; i < N; ++i) {
      values[i] =
          misc::tpMaxWaveNumber * std::exp(-misc::tpLogDz * (misc::numberOfTPGridPoints - i - 1));
    }
  }
  real const& operator[](size_t i) const { return values[i]; };

  private:
  std::array<real, N> values;
};

/**
 * Inverse Fourier coefficients on the logarithmic grid.
 */
template <size_t N>
class InverseFourierCoefficients {
  public:
  constexpr InverseFourierCoefficients() {
    GridPoints<N> localGridPoints;

    for (size_t i = 1; i < N - 1; ++i) {
      values[i] = std::sqrt(2 / M_PI) * localGridPoints[i] * misc::tpLogDz;
    }
    values[0] = std::sqrt(2 / M_PI) * localGridPoints[0] * (1 + misc::tpLogDz);
    values[N - 1] = std::sqrt(2 / M_PI) * localGridPoints[N - 1] * 0.5 * misc::tpLogDz;
  }
  real const& operator[](size_t i) const { return values[i]; };

  private:
  std::array<real, N> values;
};

/**
 * Stores the heat generation (without tauV) \f$\exp\left(\hat{l}^2/2\right) / \sqrt{2 \pi}\f$.
 */
template <size_t N>
class GaussianHeatSource {
  public:
  constexpr GaussianHeatSource() {
    GridPoints<N> localGridPoints;
    real factor = 1 / std::sqrt(2.0 * M_PI);

    for (size_t i = 0; i < N; ++i) {
      const real heatGeneration = std::exp(-0.5 * misc::power<2>(localGridPoints[i]));
      values[i] = factor * heatGeneration;
    }
  }
  real const& operator[](size_t i) const { return values[i]; };

  private:
  std::array<real, N> values;
};

/**
 * We follow Noda&Lapusta (2010) doi:10.1029/2010JB007780.
 * Define: \f$p, T\f$ pressure and temperature, \f$\Pi, \Theta\f$ fourier transform of pressure and
 * temperature respectively, \f$\Sigma = \Pi + \Lambda^\prime \Theta\f$. We solve equations (6) and
 * (7) with the method from equation(10).
 * \f[\begin{aligned}\text{Equation 6:} && \frac{\partial \Theta}{\partial t} =& -l^2 \alpha_{th}
 * \Theta + \frac{\Omega}{\rho c}\\ \text{Equation 7:} && \frac{\partial \Sigma}{\partial t} =& -l^2
 * \alpha_{hy} \Theta + (\Lambda + \Lambda^\prime) \frac{\Omega}{\rho c}\\\end{aligned}\f] with \f$
 * \Omega = \tau V \frac{\exp(-l^2 w^2 / 2) }{\sqrt{2\pi}}\f$. We define \f$\hat{l} = lw \in
 * [0,10]\f$ (see comment in [15]). Now, we can apply the solution procedure from equation (10) to
 * get:
 * \f[ \begin{aligned}\Theta(t+\Delta t) &= \frac{\Omega}{\rho c l^2 \alpha_{th}} \left[1 -
 * \exp\left(-l^2 \alpha_{th} \Delta t\right)\right] + \Theta(t)\exp(-l^2\alpha_{th} \Delta t)\\
 * &= \frac{\tau V}{\sqrt{2\pi}\rho c \left(\hat{l}/w\right)^2 \alpha_{th}} \exp(-\hat{l}^2/2)
 * \left[1 - \exp\left(-\left(\hat{l}/w\right)^2 \alpha_{th} \Delta t\right)\right]
 * + \Theta(t)\exp\left(-\left(\hat{l}/w\right)^2\alpha_{hy} \Delta t\right)\\\end{aligned}\f]
 * and
 * \f[ \begin{aligned}\Sigma(t+\Delta t) &= \frac{(\Lambda + \Lambda^\prime)\Omega}{\rho c l^2
 * \alpha_{hy}} \left[1 - \exp\left(-l^2 \alpha_{hy} \Delta t\right)\right] +
 * \Theta(t)\exp(-l^2\alpha_{hy} \Delta t)\\
 * &= \frac{(\Lambda + \Lambda^\prime)\tau V}{\sqrt{2\pi}\rho c \left(\hat{l}/w\right)^2
 * \alpha_{hy}} \exp(-\hat{l}^2/2) \left[1 - \exp\left(-\left(\hat{l}/w\right)^2 \alpha_{hy} \Delta
 * t\right)\right]
 * + \Sigma(t)\exp\left(-\left(\hat{l}/w\right)^2\alpha_{th} \Delta t\right)\\\end{aligned}\f]
 * We then compute the pressure and temperature update with an inverse Fourier transform from
 * \f$\Pi, \Theta\f$.
 */
class ThermalPressurization {
  public:
  explicit ThermalPressurization(DRParameters* drParameters) : drParameters(drParameters){};

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime);

  /**
   * Compute thermal pressure according to Noda&Lapusta (2010) at all Gauss Points within one face
   * bool saveTmpInTP is used to save final values for Theta and Sigma in the LTS tree
   */
  void calcFluidPressure(std::array<real, misc::numPaddedPoints> const& normalStress,
                         real const (*mu)[misc::numPaddedPoints],
                         std::array<real, misc::numPaddedPoints> const& slipRateMagnitude,
                         real deltaT,
                         bool saveTPinLTS,
                         unsigned int timeIndex,
                         unsigned int ltsFace);

  real getFluidPressure(unsigned int ltsFace, unsigned int pointIndex) const {
    return pressure[ltsFace][pointIndex];
  }

  protected:
  real (*temperature)[misc::numPaddedPoints];
  real (*pressure)[misc::numPaddedPoints];
  real (*theta)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*sigma)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*thetaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*sigmaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*halfWidthShearZone)[misc::numPaddedPoints];
  real (*hydraulicDiffusivity)[misc::numPaddedPoints];
  real (*faultStrength)[misc::numPaddedPoints];

  private:
  DRParameters* drParameters;

  /**
   * Compute temperature and pressure update according to Noda&Lapusta (2010) on one Gaus point.
   */
  void updateTemperatureAndPressure(real slipRateMagnitude,
                                    real deltaT,
                                    unsigned int pointIndex,
                                    unsigned int timeIndex,
                                    unsigned int ltsFace);
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_THERMALPRESSURIZATION_H
