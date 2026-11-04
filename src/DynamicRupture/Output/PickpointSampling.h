// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_PICKPOINTSAMPLING_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_PICKPOINTSAMPLING_H_

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace seissol::dr::output {

/**
 * When the on-fault receivers take their samples: at the multiples of the sampling interval, and
 * at the end of the simulation. The times do not depend on the time steps. A time step takes the
 * samples in (its start, its end], each evaluated within the step; so the faces of all clusters
 * sample at the same times. The sample at the start of the simulation is taken before the first
 * step.
 */
class PickpointSampling {
  public:
  /**
   * A sample counts as at or before a time if it lies within this fraction of the sampling
   * interval after it. The end of a step and the start of the next one may differ in their last
   * bits; with the tolerance, both agree about a sample at their common time.
   */
  static constexpr double Tolerance = 1e-6;

  PickpointSampling() = default;

  /**
   * @param interval the time between two samples; zero or infinite for no samples but the one at
   *                 the end
   * @param endTime the end of the simulation
   */
  PickpointSampling(double interval, double endTime) : interval_(interval), endTime_(endTime) {}

  [[nodiscard]] double interval() const { return interval_; }

  /**
   * The sample times of the step of size `stepSize` that starts at `stepStart`, in ascending order.
   */
  [[nodiscard]] std::vector<double> times(double stepStart, double stepSize) const {
    std::vector<double> result;
    const double stepEnd = stepStart + stepSize;
    if (hasGrid()) {
      for (auto index = lastIndex(stepStart) + 1; index <= lastIndex(stepEnd); ++index) {
        result.push_back(static_cast<double>(index) * interval_);
      }
    }
    if (endsWithin(stepStart, stepEnd)) {
      result.push_back(endTime_);
    }
    return result;
  }

  /**
   * Whether the step of size `stepSize` that starts at `stepStart` takes any sample; the same as
   * whether times() is not empty.
   */
  [[nodiscard]] bool any(double stepStart, double stepSize) const {
    const double stepEnd = stepStart + stepSize;
    return (hasGrid() && lastIndex(stepEnd) > lastIndex(stepStart)) ||
           endsWithin(stepStart, stepEnd);
  }

  private:
  [[nodiscard]] bool hasGrid() const { return interval_ > 0 && std::isfinite(interval_); }

  /// The index of the last sample of the grid at or before `time`.
  [[nodiscard]] std::int64_t lastIndex(double time) const {
    return static_cast<std::int64_t>(std::floor(time / interval_ + Tolerance));
  }

  /// Whether the end of the simulation lies in (stepStart, stepEnd], and is no sample of the grid.
  [[nodiscard]] bool endsWithin(double stepStart, double stepEnd) const {
    const double tolerance = Tolerance * (stepEnd - stepStart);
    if (endTime_ <= stepStart + tolerance || endTime_ > stepEnd + tolerance) {
      return false;
    }
    if (!hasGrid()) {
      return true;
    }
    const double position = endTime_ / interval_;
    return std::abs(position - std::round(position)) > Tolerance;
  }

  double interval_{std::numeric_limits<double>::infinity()};
  double endTime_{std::numeric_limits<double>::infinity()};
};

/**
 * A part of the friction state that changes continuously, at `position` within a time step (0 at
 * its start, 1 at its end): linear between its values at the start and the end of the step,
 * exactly these at both ends, and bitwise the same where it does not change within the step (the
 * linear formula alone can round it off by a unit in the last place, once the compiler fuses it
 * into a multiply-add).
 */
template <typename T>
T stateInStep(T start, T end, double position) {
  if (start == end || position >= 1) {
    return end;
  }
  const auto weight = static_cast<T>(position);
  return (1 - weight) * start + weight * end;
}

/**
 * The time of an event of the friction state (e.g. the rupture) as seen at `time`: zero, as before
 * any event, if it happens only later.
 */
template <typename T>
T eventBy(T eventTime, double time) {
  return eventTime > time ? static_cast<T>(0) : eventTime;
}

} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_PICKPOINTSAMPLING_H_
