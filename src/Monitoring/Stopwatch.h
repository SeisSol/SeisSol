// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
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
#include "utils/logger.h"
#include <time.h>

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
  struct timespec startTime{};

  /** Time already spent */
  long long time{0};

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
  void printTime(const char* text, MPI_Comm comm = MPI_COMM_NULL) const;

  static void print(const char* text, double time, MPI_Comm comm = MPI_COMM_NULL);
};

} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_STOPWATCH_H_
