// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Stopwatch.h"

#include "Numerical/Statistics.h"
#include "Unit.h"
#include "utils/logger.h"

#include <string>
#include <time.h>

namespace seissol {

Stopwatch::Stopwatch() = default;

/**
 * Reset the stopwatch to zero
 */
void Stopwatch::reset() { time = 0; }

/**
 * starts the time measuring
 */
void Stopwatch::start() { (void)clock_gettime(CLOCK_MONOTONIC, &startTime); }

/**
 * get time measuring
 *
 * @return measured time (until now) in seconds
 */
double Stopwatch::split() {
  struct timespec end{};
  (void)clock_gettime(CLOCK_MONOTONIC, &end);

  return seconds(difftime(startTime, end));
}

/**
 * pauses the measuring
 *
 * @return measured time (until now) in seconds
 */
double Stopwatch::pause() {
  struct timespec end{};
  (void)clock_gettime(CLOCK_MONOTONIC, &end);

  time += difftime(startTime, end);
  return seconds(time);
}

/**
 * stops time measuring
 *
 * @return measured time in seconds
 */
double Stopwatch::stop() {
  const double time = pause();
  reset();
  return time;
}

/**
 * Collective operation, printing avg, min and max time
 */
void Stopwatch::printTime(const std::string& text) const { print(text, seconds(time)); }

void Stopwatch::print(const std::string& text, double time) {
  const auto summary = statistics::parallelSummary(time);

  logInfo() << text.c_str() << UnitTime.formatTime(summary.mean).c_str()
            << "(per rank:" << UnitTime.formatScientific(summary.mean, summary.std).c_str()
            << "; range: [" << UnitTime.formatScientific(summary.min).c_str() << ","
            << UnitTime.formatScientific(summary.max).c_str() << "])";
}

} // namespace seissol
