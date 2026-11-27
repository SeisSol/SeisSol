// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law {

/**
 * The affine law that one step of the thermal pressurization reduces to at the points of a range;
 * see ThermalPressurization::prepareFluidPressure. It belongs to the face that is evaluated, not to
 * the law, which serves all faces of the layer at once.
 */
template <typename Cfg, common::RangeType Type>
struct FluidPressureStep {
  PointValues<Cfg, Type> temperatureOffset{};
  PointValues<Cfg, Type> temperatureSlope{};
  PointValues<Cfg, Type> pressureOffset{};
  PointValues<Cfg, Type> pressureSlope{};
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
 *
 * Once for the host and the device: the wavenumber grid and the coefficients on it come with the
 * arguments of the step (FrictionLawArgs::tpGridPoints, ...).
 */
template <typename Cfg>
class ThermalPressurization {
  public:
  using real = Real<Cfg>;

  template <common::RangeType Type>
  using StepCoefficients = FluidPressureStep<Cfg, Type>;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  static void copyStorageToLocal(FrictionLawData<Cfg>* data,
                                 DynamicRupture::Layer& layerData,
                                 seissol::initializer::AllocationPlace place) {
    data->temperature = layerData.var<LTSThermalPressurization::Temperature>(Cfg(), place);
    data->pressure = layerData.var<LTSThermalPressurization::Pressure>(Cfg(), place);
    data->theta = layerData.var<LTSThermalPressurization::Theta>(Cfg(), place);
    data->sigma = layerData.var<LTSThermalPressurization::Sigma>(Cfg(), place);
    data->halfWidthShearZone =
        layerData.var<LTSThermalPressurization::HalfWidthShearZone>(Cfg(), place);
    data->hydraulicDiffusivity =
        layerData.var<LTSThermalPressurization::HydraulicDiffusivity>(Cfg(), place);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real getFluidPressure(const FrictionLawContext<Cfg, Type>& ctx,
                                                  std::uint32_t pointIndex) {
    return ctx.data->pressure[ctx.ltsFace][pointIndex];
  }

  /**
   * Walk the wavenumber grid once and reduce the step to two coefficients per point.
   *
   * Everything the grid contributes is either carried over from the stored spectral state or
   * generated by this step's shear heating, and the generation is linear in tau * V. The step is
   * therefore an affine function of that one scalar,
   *   p(tau V) = pressureOffset + pressureSlope * tau V,
   * and the same for the temperature. Whatever slip rate the solve settles on, the pressure
   * follows from these coefficients without touching the grid again.
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void prepareFluidPressure(FrictionLawContext<Cfg, Type>& ctx,
                                                      std::uint32_t timeIndex,
                                                      StepCoefficients<Type>& step) {
    const auto& drParameters = ctx.data->drParameters;
    const real deltaT = ctx.args->deltaT[timeIndex];

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      real temperatureDiffusion = 0.0;
      real temperatureSlope = 0.0;
      real pressureDiffusion = 0.0;
      real pressureSlope = 0.0;

      const GridLoop grid(ctx, deltaT, pointIndex);

      const real lambdaPrime = grid.undrainedTPResponse * grid.thermalDiffusivity /
                               (grid.hydraulicDiffusivity - grid.thermalDiffusivity);

      for (std::uint32_t tpGridPointIndex = 0; tpGridPointIndex < grid.count; ++tpGridPointIndex) {
        // Gaussian shear zone in spectral domain, normalized by w
        // \hat{l} / w
        const real squaredNormalizedTpGrid =
            misc::power<2>(grid.gridPoints[tpGridPointIndex] / grid.halfWidthShearZone);

        // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
        const real thetaTpGrid = grid.thermalDiffusivity * squaredNormalizedTpGrid;
        const real sigmaTpGrid = grid.hydraulicDiffusivity * squaredNormalizedTpGrid;
        const real preExpTheta = -thetaTpGrid * grid.deltaT;
        const real preExpSigma = -sigmaTpGrid * grid.deltaT;
        const real exp1mTheta = -std::expm1(preExpTheta);
        const real exp1mSigma = -std::expm1(preExpSigma);

        const real scaledInverseFourierCoefficient =
            grid.inverseFourierCoefficients[tpGridPointIndex] / grid.halfWidthShearZone;

        // The spectral state carries over as + F(t) exp(-A dt) in equation (10), which knows
        // nothing of this step's shear heating.
        temperatureDiffusion += scaledInverseFourierCoefficient *
                                grid.theta[tpGridPointIndex * GridLoop::Stride] *
                                std::exp(preExpTheta);
        pressureDiffusion += scaledInverseFourierCoefficient *
                             grid.sigma[tpGridPointIndex * GridLoop::Stride] *
                             std::exp(preExpSigma);

        // The generation B/A * (1 - exp(-A dt)) is linear in tau * V, so the whole dependence of
        // the step on the slip rate sits in one scalar and the grid does not have to be walked
        // again. heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
        const real omegaSlope = grid.heatSource[tpGridPointIndex];
        temperatureSlope += scaledInverseFourierCoefficient * omegaSlope /
                            (grid.heatCapacity * thetaTpGrid) * exp1mTheta;
        pressureSlope += scaledInverseFourierCoefficient * omegaSlope *
                         (grid.undrainedTPResponse + lambdaPrime) /
                         (grid.heatCapacity * sigmaTpGrid) * exp1mSigma;
      }

      // Update pore pressure change: sigma = pore pressure + lambda' * temperature
      step.temperatureOffset[pointIndex] = temperatureDiffusion + drParameters.initialTemperature;
      step.temperatureSlope[pointIndex] = temperatureSlope;
      step.pressureOffset[pointIndex] =
          -(pressureDiffusion - lambdaPrime * temperatureDiffusion) + drParameters.initialPressure;
      step.pressureSlope[pointIndex] = -(pressureSlope - lambdaPrime * temperatureSlope);
    });
  }

  /// the affine law of this step: p(tau V) = offset + slope * tau V; see prepareFluidPressure
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real fluidPressureOffset(const StepCoefficients<Type>& step,
                                                     std::uint32_t pointIndex) {
    return step.pressureOffset[pointIndex];
  }
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real fluidPressureSlope(const StepCoefficients<Type>& step,
                                                    std::uint32_t pointIndex) {
    return step.pressureSlope[pointIndex];
  }

  /// Evaluate the affine law of the preceding pass at the shear heating of the given state, but at
  /// the points that `frozen` marks.
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      applyShearHeating(FrictionLawContext<Cfg, Type>& ctx,
                        const PointValues<Cfg, Type>& normalStress,
                        const PointValues<Cfg, Type>& slipRateMagnitude,
                        const StepCoefficients<Type>& step,
                        const PointArray<std::int32_t, Cfg, Type>* frozen = nullptr) {
    const auto ltsFace = ctx.ltsFace;
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      if (frozen != nullptr && (*frozen)[pointIndex] != 0) {
        return;
      }
      const real tauV = shearHeating(ctx.data->mu[ltsFace][pointIndex],
                                     normalStress[pointIndex],
                                     slipRateMagnitude[pointIndex]);
      ctx.data->temperature[ltsFace][pointIndex] =
          step.temperatureOffset[pointIndex] + step.temperatureSlope[pointIndex] * tauV;
      ctx.data->pressure[ltsFace][pointIndex] =
          step.pressureOffset[pointIndex] + step.pressureSlope[pointIndex] * tauV;
    });
  }

  /// Advance the spectral state by the shear heating the solve settled on, and store it.
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      finalizeFluidPressure(FrictionLawContext<Cfg, Type>& ctx,
                            const PointValues<Cfg, Type>& normalStress,
                            const PointValues<Cfg, Type>& slipRateMagnitude,
                            std::uint32_t timeIndex,
                            const StepCoefficients<Type>& step) {
    const real deltaT = ctx.args->deltaT[timeIndex];
    const auto ltsFace = ctx.ltsFace;

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      const real tauV = shearHeating(ctx.data->mu[ltsFace][pointIndex],
                                     normalStress[pointIndex],
                                     slipRateMagnitude[pointIndex]);

      const GridLoop grid(ctx, deltaT, pointIndex);

      const real lambdaPrime = grid.undrainedTPResponse * grid.thermalDiffusivity /
                               (grid.hydraulicDiffusivity - grid.thermalDiffusivity);

      for (std::uint32_t tpGridPointIndex = 0; tpGridPointIndex < grid.count; ++tpGridPointIndex) {
        const real squaredNormalizedTpGrid =
            misc::power<2>(grid.gridPoints[tpGridPointIndex] / grid.halfWidthShearZone);
        const real thetaTpGrid = grid.thermalDiffusivity * squaredNormalizedTpGrid;
        const real sigmaTpGrid = grid.hydraulicDiffusivity * squaredNormalizedTpGrid;
        const real preExpTheta = -thetaTpGrid * grid.deltaT;
        const real preExpSigma = -sigmaTpGrid * grid.deltaT;

        // Noda & Lapusta (2010) equation (10), F(t) exp(-A dt) + B/A (1 - exp(-A dt)), written as
        // F(t) + (B/A - F(t)) (1 - exp(-A dt)) with expm1 alone: a mode at its steady state B/A
        // stays there in any precision, rather than drifting by the mismatch of exp and expm1 in
        // every step.
        const real omega = tauV * grid.heatSource[tpGridPointIndex];
        const real thetaSteady = omega / (grid.heatCapacity * thetaTpGrid);
        const real sigmaSteady =
            omega * (grid.undrainedTPResponse + lambdaPrime) / (grid.heatCapacity * sigmaTpGrid);
        const real thetaOld = grid.theta[tpGridPointIndex * GridLoop::Stride];
        const real sigmaOld = grid.sigma[tpGridPointIndex * GridLoop::Stride];
        grid.theta[tpGridPointIndex * GridLoop::Stride] =
            thetaOld + (thetaSteady - thetaOld) * -std::expm1(preExpTheta);
        grid.sigma[tpGridPointIndex * GridLoop::Stride] =
            sigmaOld + (sigmaSteady - sigmaOld) * -std::expm1(preExpSigma);
      }

      // the coefficients of the preceding pass describe this very state, so the output follows
      // from them rather than from a second sum over the grid
      ctx.data->temperature[ltsFace][pointIndex] =
          step.temperatureOffset[pointIndex] + step.temperatureSlope[pointIndex] * tauV;
      ctx.data->pressure[ltsFace][pointIndex] =
          step.pressureOffset[pointIndex] + step.pressureSlope[pointIndex] * tauV;
    });
  }

  private:
  /**
   * What the walk over the wavenumber grid of a point reads, taken from the context once before the
   * walk. The walk calls exp and expm1, which, as far as the compiler knows, may change any memory
   * that a pointer reaches, the context included; so it would load the parameters, the tables and
   * the data of the point anew in every step of the grid, through the context, unless they are
   * held in variables of the point's own.
   */
  struct GridLoop {
    /// the distance between the values of two neighboring wavenumbers of a point in the spectral
    /// state (theta, sigma)
    static constexpr std::size_t Stride = misc::NumPaddedPoints<Cfg>;

    std::uint32_t count;
    real deltaT;
    real thermalDiffusivity;
    real heatCapacity;
    real undrainedTPResponse;
    real hydraulicDiffusivity;
    real halfWidthShearZone;
    const real* __restrict gridPoints;
    const real* __restrict inverseFourierCoefficients;
    const real* __restrict heatSource;
    /// the spectral state of the point at the first wavenumber
    real* __restrict theta;
    real* __restrict sigma;

    template <common::RangeType Type>
    SEISSOL_HOSTDEVICE GridLoop(const FrictionLawContext<Cfg, Type>& ctx,
                                real stepDeltaT,
                                std::uint32_t pointIndex)
        : count(ctx.data->drParameters.tpGridPoints), deltaT(stepDeltaT),
          thermalDiffusivity(ctx.data->drParameters.thermalDiffusivity),
          heatCapacity(ctx.data->drParameters.heatCapacity),
          undrainedTPResponse(ctx.data->drParameters.undrainedTPResponse),
          hydraulicDiffusivity(ctx.data->hydraulicDiffusivity[ctx.ltsFace][pointIndex]),
          halfWidthShearZone(ctx.data->halfWidthShearZone[ctx.ltsFace][pointIndex]),
          gridPoints(ctx.args->tpGridPoints),
          inverseFourierCoefficients(ctx.args->tpInverseFourierCoefficients),
          heatSource(ctx.args->heatSource),
          theta(&ctx.data->theta[static_cast<std::size_t>(ctx.ltsFace) * count][pointIndex]),
          sigma(&ctx.data->sigma[static_cast<std::size_t>(ctx.ltsFace) * count][pointIndex]) {}
  };

  /// tau * V, the shear heating a point produces at the given state
  SEISSOL_HOSTDEVICE static real shearHeating(real mu, real normalStress, real slipRate) {
    return -mu * normalStress * slipRate;
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_THERMALPRESSURIZATION_THERMALPRESSURIZATION_H_
