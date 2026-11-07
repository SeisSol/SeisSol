// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_TIMEPOINTSAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_TIMEPOINTSAUX_H_

#include <cstddef>

// Not in seissol::kernels::device, as the other device helpers are: the kernels of the boundaries
// include this in headers, where that namespace would hide ::device from the kernel classes.
namespace seissol::kernels::timepoints {

/// The most points into a step that timePoints() takes.
inline constexpr std::size_t MaxTimeOffsets = 16;

/// Points into a step, relative to its start: passed by value, so that a replay of the work that
/// launches timePoints() (e.g. from a graph) takes them along.
struct TimeOffsets {
  double values[MaxTimeOffsets]{};
  std::size_t count{0};
};

/**
 * The times of the points `offsets` into the step that starts at the time of the clock on the
 * device (solver::ClusterClock::device()), in the order of the stream: times[i] = *clock +
 * offsets.values[i], in device memory. Work after it on the stream reads the times of the step it
 * runs in, also when it is replayed.
 */
void timePoints(const double* clock, const TimeOffsets& offsets, double* times, void* stream);

} // namespace seissol::kernels::timepoints

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_TIMEPOINTSAUX_H_
