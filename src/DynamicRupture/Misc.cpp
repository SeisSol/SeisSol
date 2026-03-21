// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Misc.h"

#include "Geometry/MeshDefinition.h"
#include "Initializer/Parameters/DRParameters.h"

#include <cmath>
#include <string>
#include <utils/logger.h>

namespace seissol::dr::misc {
template <>
double power<8>(double base) {
  const double baseSquared = base * base;
  const double baseToTheFour = baseSquared * baseSquared;
  const double baseToTheEight = baseToTheFour * baseToTheFour;
  return baseToTheEight;
}

void computeStrikeAndDipVectors(const VrtxCoords normal, VrtxCoords strike, VrtxCoords dip) {
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
    const seissol::initializer::parameters::DRParameters& parameters) {
  this->healingThreshold = parameters.healingThreshold;
  this->tpProxyExponent = parameters.tpProxyExponent;
  this->rsF0 = parameters.rsF0;
  this->rsB = parameters.rsB;
  this->rsSr0 = parameters.rsSr0;
  this->rsInitialSlipRate1 = parameters.rsInitialSlipRate1;
  this->rsInitialSlipRate2 = parameters.rsInitialSlipRate2;
  this->muW = parameters.muW;
  this->thermalDiffusivity = parameters.thermalDiffusivity;
  this->heatCapacity = parameters.heatCapacity;
  this->undrainedTPResponse = parameters.undrainedTPResponse;
  this->initialTemperature = parameters.initialTemperature;
  this->initialPressure = parameters.initialPressure;
  this->vStar = parameters.vStar;
  this->prakashLength = parameters.prakashLength;
  this->terminatorSlipRateThreshold = parameters.terminatorSlipRateThreshold;
  this->etaDamp = parameters.etaDamp;
  this->etaDampEnd = parameters.etaDampEnd;
  std::copy(parameters.t0.begin(), parameters.t0.end(), this->t0.begin());
  std::copy(parameters.s0.begin(), parameters.s0.end(), this->s0.begin());
  this->nucleationCount = parameters.nucleationCount;
  this->rsMaxNumberSlipRateUpdates = parameters.rsMaxNumberSlipRateUpdates;
  this->rsNumberStateVariableUpdates = parameters.rsNumberStateVariableUpdates;
  this->rsNewtonTolerance = parameters.rsNewtonTolerance;
  this->rsStateTolerance = parameters.rsStateTolerance;

  this->isFrictionEnergyRequired = parameters.isFrictionEnergyRequired;
  this->isCheckAbortCriteraEnabled = parameters.isCheckAbortCriteraEnabled;
  this->energiesFromAcrossFaultVelocities = parameters.energiesFromAcrossFaultVelocities;
}
} // namespace seissol::dr
