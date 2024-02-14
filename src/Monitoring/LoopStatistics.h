/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 */

#ifndef MONITORING_LOOPSTATISTICS_H_
#define MONITORING_LOOPSTATISTICS_H_

#include <cassert>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <vector>
#include "Parallel/MPI.h"

namespace seissol {
class LoopStatistics {
  public:
  void enableSampleOutput(bool enabled);

  void addRegion(std::string const& name, bool includeInSummary = true);

  unsigned getRegion(std::string const& name) const;

  void begin(unsigned region);

  void end(unsigned region, unsigned numIterations, unsigned subRegion);

  void addSample(
      unsigned region, unsigned numIterations, unsigned subRegion, timespec begin, timespec end);

  void reset();

  void printSummary(MPI_Comm comm);

  void writeSamples(const std::string& outputPrefix, bool isLoopStatisticsNetcdfOutputOn);

  private:
  struct Sample {
    timespec begin;
    timespec end;
    unsigned numIters;
    unsigned subRegion;
  };

  struct StatisticVariables {
    double x = 0;
    double x2 = 0;
    double xy = 0;
    double y = 0;
    double y2 = 0;
    unsigned long long n = 0;
  };

  struct Region {
    std::string name;
    std::vector<Sample> times;
    bool includeInSummary;
    timespec begin;
    StatisticVariables variables;

    Region(const std::string& name, bool includeInSummary);
  };

  std::vector<Region> regions;
  bool outputSamples = false;
};
} // namespace seissol

#endif // MONITORING_LOOPSTATISTICS_H_
