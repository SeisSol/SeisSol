// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_

#include <DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h>
#include <array>
#include <cstddef>

#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"

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
  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::ThermalPressurization*>(dynRup);
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->temperature = layerData.var(concreteLts->temperature, place);
    data->pressure = layerData.var(concreteLts->pressure, place);
    data->theta = layerData.var(concreteLts->theta, place);
    data->sigma = layerData.var(concreteLts->sigma, place);
    data->thetaTmpBuffer = layerData.var(concreteLts->thetaTmpBuffer, place);
    data->sigmaTmpBuffer = layerData.var(concreteLts->sigmaTmpBuffer, place);
    data->faultStrength = layerData.var(concreteLts->faultStrength, place);
    data->halfWidthShearZone = layerData.var(concreteLts->halfWidthShearZone, place);
    data->hydraulicDiffusivity = layerData.var(concreteLts->hydraulicDiffusivity, place);
  }

  /**
   * Compute thermal pressure according to Noda&Lapusta (2010) at all Gauss Points within one face
   * bool saveTmpInTP is used to save final values for Theta and Sigma in the LTS tree
   */
  SEISSOL_DEVICE static void
      calcFluidPressure(FrictionLawContext& ctx, int timeIndex, bool saveTmpInTP) {
    ctx.data->faultStrength[ctx.ltsFace][ctx.pointIndex] =
        -ctx.data->mu[ctx.ltsFace][ctx.pointIndex] * ctx.initialVariables.normalStress;

    for (int i = 0; i < misc::NumTpGridPoints; ++i) {
      ctx.data->thetaTmpBuffer[ctx.ltsFace][ctx.pointIndex][i] =
          ctx.data->theta[ctx.ltsFace][ctx.pointIndex][i];
      ctx.data->sigmaTmpBuffer[ctx.ltsFace][ctx.pointIndex][i] =
          ctx.data->sigma[ctx.ltsFace][ctx.pointIndex][i];
    }

    // use Theta/Sigma from last timestep
    updateTemperatureAndPressure(ctx, timeIndex);

    // copy back to LTS tree, if necessary
    if (saveTmpInTP) {
      for (int i = 0; i < misc::NumTpGridPoints; ++i) {
        ctx.data->theta[ctx.ltsFace][ctx.pointIndex][i] =
            ctx.data->thetaTmpBuffer[ctx.ltsFace][ctx.pointIndex][i];
        ctx.data->sigma[ctx.ltsFace][ctx.pointIndex][i] =
            ctx.data->sigmaTmpBuffer[ctx.ltsFace][ctx.pointIndex][i];
      }
    }
  }

  SEISSOL_DEVICE static real getFluidPressure(FrictionLawContext& ctx) {
    return ctx.data->pressure[ctx.ltsFace][ctx.pointIndex];
  }

  /**
   * Compute temperature and pressure update according to Noda&Lapusta (2010) on one Gaus point.
   */
  SEISSOL_DEVICE static void updateTemperatureAndPressure(FrictionLawContext& ctx,
                                                          unsigned int timeIndex) {
    real temperatureUpdate = 0.0;
    real pressureUpdate = 0.0;

    const real tauV = ctx.data->faultStrength[ctx.ltsFace][ctx.pointIndex] *
                      ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];
    const real lambdaPrime = ctx.data->drParameters.undrainedTPResponse *
                             ctx.data->drParameters.thermalDiffusivity /
                             (ctx.data->hydraulicDiffusivity[ctx.ltsFace][ctx.pointIndex] -
                              ctx.data->drParameters.thermalDiffusivity);

    for (int tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints; tpGridPointIndex++) {
      // Gaussian shear zone in spectral domain, normalized by w
      // \hat{l} / w
      const real squaredNormalizedTpGrid =
          misc::power<2>(ctx.TpGridPoints[tpGridPointIndex] /
                         ctx.data->halfWidthShearZone[ctx.ltsFace][ctx.pointIndex]);

      // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
      const real thetaTpGrid = ctx.data->drParameters.thermalDiffusivity * squaredNormalizedTpGrid;
      const real sigmaTpGrid =
          ctx.data->hydraulicDiffusivity[ctx.ltsFace][ctx.pointIndex] * squaredNormalizedTpGrid;
      const real preExpTheta = -thetaTpGrid * ctx.data->deltaT[timeIndex];
      const real preExpSigma = -sigmaTpGrid * ctx.data->deltaT[timeIndex];
      const real expTheta = std::exp(preExpTheta);
      const real expSigma = std::exp(preExpSigma);
      const real exp1mTheta = -std::expm1(preExpTheta);
      const real exp1mSigma = -std::expm1(preExpSigma);

      // Temperature and pressure diffusion in spectral domain over timestep
      // This is + F(t) exp(-A dt) in equation (10)
      const real thetaDiffusion =
          ctx.data->thetaTmpBuffer[ctx.ltsFace][ctx.pointIndex][tpGridPointIndex] * expTheta;
      const real sigmaDiffusion =
          ctx.data->sigmaTmpBuffer[ctx.ltsFace][ctx.pointIndex][tpGridPointIndex] * expSigma;

      // Heat generation during timestep
      // This is B/A * (1 - exp(-A dt)) in Noda & Lapusta (2010) equation (10)
      // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
      const real omega = tauV * ctx.HeatSource[tpGridPointIndex];
      const real thetaGeneration =
          omega / (ctx.data->drParameters.heatCapacity * thetaTpGrid) * exp1mTheta;
      const real sigmaGeneration = omega *
                                   (ctx.data->drParameters.undrainedTPResponse + lambdaPrime) /
                                   (ctx.data->drParameters.heatCapacity * sigmaTpGrid) * exp1mSigma;

      // Sum both contributions up
      ctx.data->thetaTmpBuffer[ctx.ltsFace][ctx.pointIndex][tpGridPointIndex] =
          thetaDiffusion + thetaGeneration;
      ctx.data->sigmaTmpBuffer[ctx.ltsFace][ctx.pointIndex][tpGridPointIndex] =
          sigmaDiffusion + sigmaGeneration;

      // Recover temperature and altered pressure using inverse Fourier transformation from the new
      // contribution
      const real scaledInverseFourierCoefficient =
          ctx.TpInverseFourierCoefficients[tpGridPointIndex] /
          ctx.data->halfWidthShearZone[ctx.ltsFace][ctx.pointIndex];
      temperatureUpdate += scaledInverseFourierCoefficient *
                           ctx.data->thetaTmpBuffer[ctx.ltsFace][ctx.pointIndex][tpGridPointIndex];
      pressureUpdate += scaledInverseFourierCoefficient *
                        ctx.data->sigmaTmpBuffer[ctx.ltsFace][ctx.pointIndex][tpGridPointIndex];
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
