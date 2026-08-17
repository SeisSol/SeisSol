// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Alexander Heinecke
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_MONITORING_STOPWATCH_H_
#define SEISSOL_SRC_MONITORING_STOPWATCH_H_

#include "Parallel/MPI.h"

#include <time.h>
#include <utils/logger.h>

namespace seissol {

/** Returns the time difference in nanoseconds. */
inline long long difftime(const timespec& start, const timespec& end) {
  return 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
}

inline double seconds(long long time) { return 1.0e-9 * time; }

/**
 * Stopwatch
 *
 * Part of SeisSol, so you can easily calculate the needed time of SeisSol computations with a high
 * precision
 */
class Stopwatch {
  private:
  struct timespec startTime_{};

  /** Time already spent */
  long long time_{0};

  public:
  /**
   * Constructor
   *
   * resets the Stopwatch
   */
  Stopwatch();

  /**
   * Destructor
   */
  ~Stopwatch() = default;

  /**
   * Reset the stopwatch to zero
   */
  void reset();

  /**
   * starts the time measuring
   */
  void start();

  /**
   * get time measuring
   *
   * @return measured time (until now) in seconds
   */
  double split();

  /**
   * pauses the measuring
   *
   * @return measured time (until now) in seconds
   */
  double pause();

  /**
   * stops time measuring
   *
   * @return measured time in seconds
   */
  double stop();

  /**
   * Collective operation, printing avg, min and max time
   */
  void printTime(const std::string& text) const;

  static void print(const std::string& text, double time);
};

} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_STOPWATCH_H_
