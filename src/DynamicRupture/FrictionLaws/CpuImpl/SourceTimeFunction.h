// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SOURCETIMEFUNCTION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SOURCETIMEFUNCTION_H_

#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"

namespace seissol::dr::friction_law::cpu {
class YoffeSTF {
  private:
  real (*onsetTime)[misc::NumPaddedPoints];
  real (*tauS)[misc::NumPaddedPoints];
  real (*tauR)[misc::NumPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime,
                [[maybe_unused]] real timeIncrement,
                size_t ltsFace,
                size_t pointIndex);
};

class GaussianSTF {
  private:
  real (*onsetTime)[misc::NumPaddedPoints];
  real (*riseTime)[misc::NumPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex);
};

class DeltaSTF {
  private:
  real (*onsetTime)[misc::NumPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex);
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SOURCETIMEFUNCTION_H_
