// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_

#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"

namespace seissol::dr::friction_law::cpu {
template <typename Cfg>
class YoffeSTF {
  private:
  using real = Real<Cfg>;
  real (*__restrict onsetTime)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict tauS)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict tauR)[misc::NumPaddedPoints<Cfg>];

  public:
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  real evaluate(real currentTime,
                [[maybe_unused]] real timeIncrement,
                size_t ltsFace,
                uint32_t pointIndex);
};

template <typename Cfg>
class GaussianSTF {
  private:
  using real = Real<Cfg>;
  real (*__restrict onsetTime)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict riseTime)[misc::NumPaddedPoints<Cfg>];

  public:
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, uint32_t pointIndex);
};

template <typename Cfg>
class DeltaSTF {
  private:
  using real = Real<Cfg>;
  real (*__restrict onsetTime)[misc::NumPaddedPoints<Cfg>];

  public:
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, uint32_t pointIndex);
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_
