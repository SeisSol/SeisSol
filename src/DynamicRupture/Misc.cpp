// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Misc.h"

#include "Geometry/MeshDefinition.h"
#include "Initializer/Parameters/DRParameters.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utils/logger.h>

namespace seissol::dr::misc {
void computeStrikeAndDipVectors(const CoordinateT& normal, CoordinateT& strike, CoordinateT& dip) {
  // compute normalized strike vector
  const auto strikeInvLength = 1.0 / std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
  strike[0] = normal[1] * strikeInvLength;
  strike[1] = -normal[0] * strikeInvLength;
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

std::string frictionLawName(seissol::dr::misc::FrictionLawType type) {
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
  default:
    logError() << "unknown friction law";
    return "unknown";
  }
}

} // namespace seissol::dr::misc

namespace seissol::dr {
FrictionLawParameters::FrictionLawParameters(
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
      nucleationCount(parameters.nucleationCount),
      rsMaxNumberSlipRateUpdates(parameters.rsMaxNumberSlipRateUpdates),
      rsNumberStateVariableUpdates(parameters.rsNumberStateVariableUpdates),
      rsSlipRateTolerance(parameters.rsSlipRateTolerance),
      rsStateTolerance(parameters.rsStateTolerance),
      isFrictionEnergyRequired(parameters.isFrictionEnergyRequired),
      isCheckAbortCriteraEnabled(parameters.isCheckAbortCriteraEnabled),
      energiesFromAcrossFaultVelocities(parameters.energiesFromAcrossFaultVelocities) {

  std::copy(parameters.t0.begin(), parameters.t0.end(), this->t0.begin());
  std::copy(parameters.s0.begin(), parameters.s0.end(), this->s0.begin());
}
} // namespace seissol::dr
