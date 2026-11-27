// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Misc.h"

#include "Core/Common/Constants.h"
#include "Features/DynamicRupture/DRParameters.h"

#include <cmath>
#include <string_view>
#include <utils/logger.h>

namespace seissol::dr::misc {
void computeStrikeAndDipVectors(const CoordinateT& normal, CoordinateT& strike, CoordinateT& dip) {
  // compute normalized strike vector; std::hypot, as the squares of tiny components underflow
  const auto horizontalLength = std::hypot(normal[0], normal[1]);
  if (horizontalLength > 0.0) {
    strike[0] = normal[1] / horizontalLength;
    strike[1] = -normal[0] / horizontalLength;
  } else {
    // a horizontal fault leaves the strike open; it is taken as the limit of a fault that strikes
    // along x and dips towards y by a vanishing angle, which turns the strike around with the
    // normal, as everywhere else (see docs/left-lateral-right-lateral-normal-reverse.rst)
    strike[0] = normal[2] > 0.0 ? 1.0 : -1.0;
    strike[1] = 0.0;
  }
  strike[2] = 0.0;

  // compute normalized dip vector
  dip[0] = -1.0 * strike[1] * normal[2];
  dip[1] = strike[0] * normal[2];
  dip[2] = strike[1] * normal[0] - strike[0] * normal[1];
  const auto dipInvLength = 1.0 / std::sqrt(dip[0] * dip[0] + dip[1] * dip[1] + dip[2] * dip[2]);
  dip[0] *= dipInvLength;
  dip[1] *= dipInvLength;
  dip[2] *= dipInvLength;
}

std::string_view frictionLawName(seissol::dr::misc::FrictionLawType type) {
  switch (type) {
  case seissol::dr::misc::FrictionLawType::NoFault:
    return "nofault";
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesYoffe:
    return "imposed-yoffe";
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesGaussian:
    return "imposed-gaussian";
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesDelta:
    return "imposed-delta";
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakening:
    return "lsw-base";
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakeningBimaterial:
    return "lsw-bimaterial";
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakeningTPApprox:
    return "lsw-tpapprox";
  case seissol::dr::misc::FrictionLawType::RateAndStateAgingLaw:
    return "rs-slow-aging";
  case seissol::dr::misc::FrictionLawType::RateAndStateSlipLaw:
    return "rs-slow-slip";
  case seissol::dr::misc::FrictionLawType::RateAndStateSevereVelocityWeakening:
    return "rs-severe";
  case seissol::dr::misc::FrictionLawType::RateAndStateFastVelocityWeakening:
    return "rs-fast";
  case seissol::dr::misc::FrictionLawType::RateAndStateAgingNucleation:
    return "rs-slow-aging-nucleation";
  default:
    logError() << "unknown friction law";
    return "unknown";
  }
}

} // namespace seissol::dr::misc

namespace seissol::dr {
template <typename RealT>
FrictionLawParameters<RealT>::FrictionLawParameters(
    const seissol::initializer::parameters::DRParameters& parameters)
    : healingThreshold(parameters.healingThreshold), tpProxyExponent(parameters.tpProxyExponent),
      rsF0(parameters.rsF0), rsB(parameters.rsB), rsSr0(parameters.rsSr0),
      rsInitialSlipRate1(parameters.rsInitialSlipRate1),
      rsInitialSlipRate2(parameters.rsInitialSlipRate2), muW(parameters.muW),
      thermalDiffusivity(parameters.thermalDiffusivity), heatCapacity(parameters.heatCapacity),
      undrainedTPResponse(parameters.undrainedTPResponse),
      initialTemperature(parameters.initialTemperature),
      initialPressure(parameters.initialPressure), vStar(parameters.vStar),
      prakashLength(parameters.prakashLength),
      terminatorSlipRateThreshold(parameters.terminatorSlipRateThreshold),
      etaDamp(parameters.etaDamp), etaDampEnd(parameters.etaDampEnd),
      forcedRuptureRiseTime(static_cast<RealT>(parameters.t0[0])),
      sourceCount(stressSourceCount(parameters)), tpGridPoints(parameters.tpGridPoints),
      rsMaxNumberSlipRateUpdates(parameters.rsMaxNumberSlipRateUpdates),
      rsNumberStateVariableUpdates(parameters.rsNumberStateVariableUpdates),
      rsSlipRateTolerance(parameters.rsSlipRateTolerance),
      rsStateTolerance(parameters.rsStateTolerance),
      isFrictionEnergyRequired(parameters.isFrictionEnergyRequired),
      isCheckAbortCriteraEnabled(parameters.isCheckAbortCriteraEnabled),
      energiesFromAcrossFaultVelocities(parameters.energiesFromAcrossFaultVelocities) {}

template struct FrictionLawParameters<float>;
template struct FrictionLawParameters<double>;
} // namespace seissol::dr
