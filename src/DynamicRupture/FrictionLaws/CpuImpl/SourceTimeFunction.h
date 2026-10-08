// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_

#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"

namespace seissol::dr::friction_law::cpu {
template <typename Cfg>
class YoffeSTF {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// a time function scaling the slip of the point (cf. ScriptedSTF)
  static constexpr bool Prescribed = false;

  private:
  real (*__restrict onsetTime_)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict tauS_)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict tauR_)[misc::NumPaddedPoints<Cfg>];

  public:
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  real evaluate(real currentTime,
                [[maybe_unused]] real timeIncrement,
                size_t ltsFace,
                uint32_t pointIndex);
};

template <typename Cfg>
class GaussianSTF {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// a time function scaling the slip of the point (cf. ScriptedSTF)
  static constexpr bool Prescribed = false;

  private:
  real (*__restrict onsetTime_)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict riseTime_)[misc::NumPaddedPoints<Cfg>];

  public:
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, uint32_t pointIndex);
};

template <typename Cfg>
class DeltaSTF {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// a time function scaling the slip of the point (cf. ScriptedSTF)
  static constexpr bool Prescribed = false;

  private:
  real (*__restrict onsetTime_)[misc::NumPaddedPoints<Cfg>];

  public:
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, uint32_t pointIndex);
};

/**
 * The slip rates of a script (FL 36), along both directions of the face and per sub-step, as
 * dr::friction_law::SlipRateEvaluator wrote them before the step: see CpuImpl/ScriptedSlipRates.h.
 */
template <typename Cfg>
class ScriptedSTF {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// gives the slip rate along both directions itself
  static constexpr bool Prescribed = true;

  private:
  const real* __restrict slipRates_{nullptr};
  std::size_t points_{0};

  public:
  void copyStorageToLocal(DynamicRupture::Layer& layerData) {
    slipRates_ = reinterpret_cast<const real*>(
        layerData.var<LTSImposedSlipRatesScript::ScriptSlipRates>(Cfg()));
    points_ = layerData.size() * misc::NumPaddedPoints<Cfg>;
  }

  void slipRates(std::size_t ltsFace,
                 uint32_t pointIndex,
                 uint32_t timeIndex,
                 real& rate1,
                 real& rate2) const {
    const std::size_t point = ltsFace * misc::NumPaddedPoints<Cfg> + pointIndex;
    const auto first = static_cast<std::size_t>(2 * timeIndex) * points_ + point;
    rate1 = slipRates_[first];
    rate2 = slipRates_[first + points_];
  }
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_
