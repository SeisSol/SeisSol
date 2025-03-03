// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Misc.h"

#include "Geometry/MeshDefinition.h"
#include <Initializer/Parameters/DRParameters.h>
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

std::string frictionLawName(seissol::initializer::parameters::FrictionLawType type) {
  switch (type) {
  case seissol::initializer::parameters::FrictionLawType::NoFault:
    return "nofault";
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesYoffe:
    return "imposed-yoffe";
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesGaussian:
    return "imposed-gaussian";
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesDelta:
    return "imposed-delta";
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakening:
    return "lsw-base";
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningBimaterial:
    return "lsw-bimaterial";
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningTPApprox:
    return "lsw-tpapprox";
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingLaw:
    return "rs-slow-aging";
  case seissol::initializer::parameters::FrictionLawType::RateAndStateSlipLaw:
    return "rs-slow-slip";
  case seissol::initializer::parameters::FrictionLawType::RateAndStateSevereVelocityWeakening:
    return "rs-severe";
  case seissol::initializer::parameters::FrictionLawType::RateAndStateFastVelocityWeakening:
    return "rs-fast";
  default:
    logError() << "unknown friction law";
    return "unknown";
  }
}

} // namespace seissol::dr::misc
