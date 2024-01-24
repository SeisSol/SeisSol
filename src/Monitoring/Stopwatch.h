/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Stopwatch originally developed by A. Heinecke
 */

#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <time.h>
#include "Parallel/MPI.h"
#include "utils/logger.h"

namespace seissol {

/** Returns the time difference in nanoseconds. */
inline long long difftime(timespec const& start, timespec const& end) {
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
  struct timespec startTime;

  /** Time already spent */
  long long time;

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
};

} // namespace seissol

#endif // STOPWATCH_H
