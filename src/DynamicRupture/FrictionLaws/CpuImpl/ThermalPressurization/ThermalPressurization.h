// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_

#include <array>

#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::friction_law::cpu {

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
  explicit ThermalPressurization(seissol::initializer::parameters::DRParameters* drParameters)
      : drParameters(drParameters) {};

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* dynRup,
                          real fullUpdateTime);

  /**
   * Compute thermal pressure according to Noda&Lapusta (2010) at all Gauss Points within one face
   * bool saveTmpInTP is used to save final values for Theta and Sigma in the LTS tree
   */
  void calcFluidPressure(const std::array<real, misc::NumPaddedPoints>& normalStress,
                         const real (*mu)[misc::NumPaddedPoints],
                         const std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
                         real deltaT,
                         bool saveTPinLTS,
                         unsigned int timeIndex,
                         unsigned int ltsFace);

  [[nodiscard]] real getFluidPressure(unsigned int ltsFace, unsigned int pointIndex) const {
    return pressure[ltsFace][pointIndex];
  }

  protected:
  real (*temperature)[misc::NumPaddedPoints]{};
  real (*pressure)[misc::NumPaddedPoints]{};
  real (*theta)[misc::NumPaddedPoints][misc::NumTpGridPoints]{};
  real (*sigma)[misc::NumPaddedPoints][misc::NumTpGridPoints]{};
  real (*thetaTmpBuffer)[misc::NumPaddedPoints][misc::NumTpGridPoints]{};
  real (*sigmaTmpBuffer)[misc::NumPaddedPoints][misc::NumTpGridPoints]{};
  real (*halfWidthShearZone)[misc::NumPaddedPoints]{};
  real (*hydraulicDiffusivity)[misc::NumPaddedPoints]{};
  real (*faultStrength)[misc::NumPaddedPoints]{};

  private:
  seissol::initializer::parameters::DRParameters* drParameters;

  /**
   * Compute temperature and pressure update according to Noda&Lapusta (2010) on one Gaus point.
   */
  void updateTemperatureAndPressure(real slipRateMagnitude,
                                    real deltaT,
                                    unsigned int pointIndex,
                                    unsigned int timeIndex,
                                    unsigned int ltsFace);
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_
