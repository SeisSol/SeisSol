// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_PHYSICS_IMPEDANCES_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_PHYSICS_IMPEDANCES_H_

#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "WavePropagation/Physics/CommonDatastructures.h"

namespace seissol::dr {

/**
 * The impedances of the two sides of a fault face, from the wave speeds of their isotropic
 * materials, and the eta they give the Riemann problem of the face (in the reals of the
 * configuration `Cfg`).
 */
template <typename Cfg>
ImpedancesAndEta<Cfg> impedancesAndEta(const model::IsotropicWaveSpeeds& plus,
                                       const model::IsotropicWaveSpeeds& minus) {
  ImpedancesAndEta<Cfg> result{};
  result.zp = plus.density * plus.pWaveVelocity;
  result.zpNeig = minus.density * minus.pWaveVelocity;
  result.zs = plus.density * plus.sWaveVelocity;
  result.zsNeig = minus.density * minus.sWaveVelocity;

  result.invZp = 1 / result.zp;
  result.invZpNeig = 1 / result.zpNeig;
  result.invZs = 1 / result.zs;
  result.invZsNeig = 1 / result.zsNeig;

  result.etaP = 1.0 / (1.0 / result.zp + 1.0 / result.zpNeig);
  result.invEtaS = 1.0 / result.zs + 1.0 / result.zsNeig;
  result.etaS = 1.0 / (1.0 / result.zs + 1.0 / result.zsNeig);
  return result;
}

} // namespace seissol::dr

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_PHYSICS_IMPEDANCES_H_
