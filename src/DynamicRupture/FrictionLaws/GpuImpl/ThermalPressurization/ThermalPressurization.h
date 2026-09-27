// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"

#include <array>
#include <cstddef>

namespace seissol::dr::friction_law::gpu {

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
  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->temperature = layerData.var<LTSThermalPressurization::Temperature>(place);
    data->pressure = layerData.var<LTSThermalPressurization::Pressure>(place);
    data->theta = layerData.var<LTSThermalPressurization::Theta>(place);
    data->sigma = layerData.var<LTSThermalPressurization::Sigma>(place);
    data->halfWidthShearZone = layerData.var<LTSThermalPressurization::HalfWidthShearZone>(place);
    data->hydraulicDiffusivity =
        layerData.var<LTSThermalPressurization::HydraulicDiffusivity>(place);
  }

  SEISSOL_DEVICE static real getFluidPressure(FrictionLawContext& __restrict ctx) {
    return ctx.data->pressure[ctx.ltsFace][ctx.pointIndex];
  }

  /**
   * Compute thermal pressure according to Noda&Lapusta (2010) at all Gauss Points within one face
   * bool saveTmpInTP is used to save final values for Theta and Sigma in the storage.
   * Compute temperature and pressure update according to Noda&Lapusta (2010) on one Gaus point.
   */
  SEISSOL_DEVICE static void
      calcFluidPressure(FrictionLawContext& __restrict ctx, uint32_t timeIndex, bool saveTmpInTP) {
    real temperatureUpdate = 0.0;
    real pressureUpdate = 0.0;

    const real faultStrength =
        -ctx.data->mu[ctx.ltsFace][ctx.pointIndex] * ctx.initialVariables.normalStress;

    const real tauV = faultStrength * ctx.initialVariables.localSlipRate;
    const real lambdaPrime = ctx.data->drParameters.undrainedTPResponse *
                             ctx.data->drParameters.thermalDiffusivity /
                             (ctx.data->hydraulicDiffusivity[ctx.ltsFace][ctx.pointIndex] -
                              ctx.data->drParameters.thermalDiffusivity);

    for (uint32_t tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
         tpGridPointIndex++) {
      // Gaussian shear zone in spectral domain, normalized by w
      // \hat{l} / w
      const real squaredNormalizedTpGrid =
          misc::power<2>(ctx.args->tpGridPoints[tpGridPointIndex] /
                         ctx.data->halfWidthShearZone[ctx.ltsFace][ctx.pointIndex]);

      // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
      const real thetaTpGrid = ctx.data->drParameters.thermalDiffusivity * squaredNormalizedTpGrid;
      const real sigmaTpGrid =
          ctx.data->hydraulicDiffusivity[ctx.ltsFace][ctx.pointIndex] * squaredNormalizedTpGrid;
      const real preExpTheta = -thetaTpGrid * ctx.args->deltaT[timeIndex];
      const real preExpSigma = -sigmaTpGrid * ctx.args->deltaT[timeIndex];
      const real exp1mTheta = -std::expm1(preExpTheta);
      const real exp1mSigma = -std::expm1(preExpSigma);

      // Noda & Lapusta (2010) equation (10), F(t) exp(-A dt) + B/A (1 - exp(-A dt)), written as
      // F(t) + (B/A - F(t)) (1 - exp(-A dt)) with expm1 alone: a mode at its steady state B/A stays
      // there in any precision, rather than drifting by the mismatch of exp and expm1 in every
      // step.
      //
      // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
      const real omega = tauV * ctx.args->heatSource[tpGridPointIndex];
      const real thetaSteady = omega / (ctx.data->drParameters.heatCapacity * thetaTpGrid);
      const real sigmaSteady = omega * (ctx.data->drParameters.undrainedTPResponse + lambdaPrime) /
                               (ctx.data->drParameters.heatCapacity * sigmaTpGrid);
      const real thetaOld = ctx.data->theta[ctx.ltsFace][tpGridPointIndex][ctx.pointIndex];
      const real sigmaOld = ctx.data->sigma[ctx.ltsFace][tpGridPointIndex][ctx.pointIndex];
      const auto thetaNew = thetaOld + (thetaSteady - thetaOld) * exp1mTheta;
      const auto sigmaNew = sigmaOld + (sigmaSteady - sigmaOld) * exp1mSigma;

      // Recover temperature and altered pressure using inverse Fourier transformation from the new
      // contribution
      const real scaledInverseFourierCoefficient =
          ctx.args->tpInverseFourierCoefficients[tpGridPointIndex] /
          ctx.data->halfWidthShearZone[ctx.ltsFace][ctx.pointIndex];
      temperatureUpdate += scaledInverseFourierCoefficient * thetaNew;
      pressureUpdate += scaledInverseFourierCoefficient * sigmaNew;

      if (saveTmpInTP) {
        ctx.data->theta[ctx.ltsFace][tpGridPointIndex][ctx.pointIndex] = thetaNew;
        ctx.data->sigma[ctx.ltsFace][tpGridPointIndex][ctx.pointIndex] = sigmaNew;
      }
    }
    // Update pore pressure change: sigma = pore pressure + lambda' * temperature
    pressureUpdate -= lambdaPrime * temperatureUpdate;

    // Temperature and pore pressure change at single GP on the fault + initial values
    ctx.data->temperature[ctx.ltsFace][ctx.pointIndex] =
        temperatureUpdate + ctx.data->drParameters.initialTemperature;
    ctx.data->pressure[ctx.ltsFace][ctx.pointIndex] =
        -pressureUpdate + ctx.data->drParameters.initialPressure;
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_
