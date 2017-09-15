/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#include <unordered_map>
#include <fstream>
#include <iomanip>
#include <utils/env.h>

#include "Stopwatch.h"

namespace seissol {
class LoopStatistics {
public:
  void addRegion(std::string const& name) {
    m_regions.push_back(name);
    m_stopwatch.push_back(Stopwatch());
    m_times.push_back(std::vector<Sample>());
  }
  
  unsigned getRegion(std::string const& name) {
    auto first = m_regions.cbegin();
    auto it = std::find(first, m_regions.cend(), name);
    return std::distance(first, it);
  }
  
  void begin(unsigned region) {
    m_stopwatch[region].start();
  }
  
  void end(unsigned region, unsigned numIterations) {
    Sample sample;
    sample.time = m_stopwatch[region].stop();
    sample.numIters = numIterations;
    m_times[region].push_back(sample);
  }

#ifdef USE_MPI  
  void printSummary(MPI_Comm comm);
#endif

  void writeSamples();
  
private:
  struct Sample {
    double time;
    unsigned numIters;
  };
  
  std::vector<Stopwatch> m_stopwatch;
  std::vector<std::string> m_regions;
  std::vector<std::vector<Sample>> m_times;
};
}

#endif // MONITORING_LOOPSTATISTICS_H_
