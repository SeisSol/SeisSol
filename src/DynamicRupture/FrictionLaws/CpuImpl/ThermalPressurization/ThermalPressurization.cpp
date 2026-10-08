// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ThermalPressurization.h"

#include "Config.h"
#include "DynamicRupture/FrictionLaws/TPCommon.h"
#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {

namespace {
template <typename RealT>
const tp::GridPoints<misc::NumTpGridPoints, RealT> TpGridPoints{};
template <typename RealT>
const tp::InverseFourierCoefficients<misc::NumTpGridPoints, RealT> TpInverseFourierCoefficients{};
template <typename RealT>
const tp::GaussianHeatSource<misc::NumTpGridPoints, RealT> HeatSource{};
} // namespace

template <typename Cfg>
void ThermalPressurization<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  temperature_ = layerData.var<LTSThermalPressurization::Temperature>(Cfg());
  pressure_ = layerData.var<LTSThermalPressurization::Pressure>(Cfg());
  theta_ = layerData.var<LTSThermalPressurization::Theta>(Cfg());
  sigma_ = layerData.var<LTSThermalPressurization::Sigma>(Cfg());
  halfWidthShearZone_ = layerData.var<LTSThermalPressurization::HalfWidthShearZone>(Cfg());
  hydraulicDiffusivity_ = layerData.var<LTSThermalPressurization::HydraulicDiffusivity>(Cfg());
}

template <typename Cfg>
void ThermalPressurization<Cfg>::calcFluidPressure(
    const std::array<real, misc::NumPaddedPoints<Cfg>>& normalStress,
    const real (*mu)[misc::NumPaddedPoints<Cfg>],
    const std::array<real, misc::NumPaddedPoints<Cfg>>& slipRateMagnitude,
    real deltaT,
    bool saveTPinLTS,
    std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
    real temperatureUpdate = 0.0;
    real pressureUpdate = 0.0;

    const real faultStrength = -mu[ltsFace][pointIndex] * normalStress[pointIndex];
    const real tauV = faultStrength * slipRateMagnitude[pointIndex];
    const real lambdaPrime =
        drParameters_.undrainedTPResponse * drParameters_.thermalDiffusivity /
        (hydraulicDiffusivity_[ltsFace][pointIndex] - drParameters_.thermalDiffusivity);

    for (uint32_t tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
         ++tpGridPointIndex) {
      // Gaussian shear zone in spectral domain, normalized by w
      // \hat{l} / w
      const real squaredNormalizedTpGrid = misc::power<2>(TpGridPoints<real>[tpGridPointIndex] /
                                                          halfWidthShearZone_[ltsFace][pointIndex]);

      // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
      const real thetaTpGrid = drParameters_.thermalDiffusivity * squaredNormalizedTpGrid;
      const real sigmaTpGrid = hydraulicDiffusivity_[ltsFace][pointIndex] * squaredNormalizedTpGrid;
      const real preExpTheta = -thetaTpGrid * deltaT;
      const real preExpSigma = -sigmaTpGrid * deltaT;
      const real expTheta = std::exp(preExpTheta);
      const real expSigma = std::exp(preExpSigma);
      const real exp1mTheta = -std::expm1(preExpTheta);
      const real exp1mSigma = -std::expm1(preExpSigma);

      // Temperature and pressure diffusion in spectral domain over timestep
      // This is + F(t) exp(-A dt) in equation (10)
      const real thetaDiffusion = theta_[ltsFace][tpGridPointIndex][pointIndex] * expTheta;
      const real sigmaDiffusion = sigma_[ltsFace][tpGridPointIndex][pointIndex] * expSigma;

      // Heat generation during timestep
      // This is B/A * (1 - exp(-A dt)) in Noda & Lapusta (2010) equation (10)
      // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
      const real omega = tauV * HeatSource<real>[tpGridPointIndex];
      const real thetaGeneration = omega / (drParameters_.heatCapacity * thetaTpGrid) * exp1mTheta;
      const real sigmaGeneration = omega * (drParameters_.undrainedTPResponse + lambdaPrime) /
                                   (drParameters_.heatCapacity * sigmaTpGrid) * exp1mSigma;

      // Sum both contributions up
      const auto thetaNew = thetaDiffusion + thetaGeneration;
      const auto sigmaNew = sigmaDiffusion + sigmaGeneration;

      // Recover temperature and altered pressure using inverse Fourier transformation from the new
      // contribution
      const real scaledInverseFourierCoefficient =
          TpInverseFourierCoefficients<real>[tpGridPointIndex] /
          halfWidthShearZone_[ltsFace][pointIndex];
      temperatureUpdate += scaledInverseFourierCoefficient * thetaNew;
      pressureUpdate += scaledInverseFourierCoefficient * sigmaNew;

      if (saveTPinLTS) {
        theta_[ltsFace][tpGridPointIndex][pointIndex] = thetaNew;
        sigma_[ltsFace][tpGridPointIndex][pointIndex] = sigmaNew;
      }
    }
    // Update pore pressure change: sigma = pore pressure + lambda' * temperature
    pressureUpdate -= lambdaPrime * temperatureUpdate;

    // Temperature and pore pressure change at single GP on the fault + initial values
    temperature_[ltsFace][pointIndex] = temperatureUpdate + drParameters_.initialTemperature;
    pressure_[ltsFace][pointIndex] = -pressureUpdate + drParameters_.initialPressure;
  }
}

#define SEISSOL_INSTANTIATE(Cfg) template class ThermalPressurization<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::dr::friction_law::cpu
