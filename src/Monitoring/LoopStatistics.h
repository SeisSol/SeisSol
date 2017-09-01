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
  void printSummary(MPI_Comm comm) {
    unsigned nRegions = m_times.size();
    double* sums = new double[5*nRegions];
    for (unsigned region = 0; region < nRegions; ++region) {
      double x = 0.0, x2 = 0.0, xy = 0.0, y = 0.0;
      unsigned N = 0;
    
      for (auto const& sample : m_times[region]) {
        if (sample.numIters > 0) {
          x  += sample.numIters;
          x2 += static_cast<double>(sample.numIters) * static_cast<double>(sample.numIters);
          xy += static_cast<double>(sample.numIters) * sample.time;
          y  += sample.time;
          ++N;
        }
      }

      sums[5*region + 0] = x;
      sums[5*region + 1] = x2;
      sums[5*region + 2] = xy;
      sums[5*region + 3] = y;
      sums[5*region + 4] = N;
    }

    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, sums, 5*nRegions, MPI_DOUBLE, MPI_SUM, 0, comm);
    } else {
      MPI_Reduce(sums, 0L, 5*nRegions, MPI_DOUBLE, MPI_SUM, 0, comm);
    }

    if (rank == 0) {
      logInfo(rank) << "Regression analysis of compute kernels:";
      for (unsigned region = 0; region < nRegions; ++region) {
        double x  = sums[5*region + 0];
        double x2 = sums[5*region + 1];
        double xy = sums[5*region + 2];
        double y  = sums[5*region + 3];
        double N  = sums[5*region + 4];

        double det = N*x2 - x*x;
        double regressionCoeffs[2];
        regressionCoeffs[0] = (x2*y - x*xy) / det;
        regressionCoeffs[1] = (-x*y + N*xy) / det;
      
        char const* names[] = { "constant", "per element"};
        for (unsigned c = 0; c < 2; ++c) {
          logInfo(rank) << m_regions[region]
                        << "(" << names[c] << "):"
                        << regressionCoeffs[c];
        }
      }
    }

    delete[] sums;
  }
#endif
  
  void writeSamplesToCSV() {
    std::string loopStatFile = utils::Env::get("SEISSOL_LOOP_STAT_PREFIX", "");
    if (!loopStatFile.empty()) {
      unsigned nRegions = m_times.size();
      for (unsigned region = 0; region < nRegions; ++region) {
        std::ofstream file;
        std::stringstream ss;
        ss << loopStatFile << "_" << seissol::MPI::mpi.rank() << "_" << m_regions[region] << ".csv";
        file.open(ss.str());
        file << "loopLength,time\n";
        for (auto const& sample : m_times[region]) {
          file << std::setprecision(20) << sample.numIters << "," << sample.time << "\n";
        }
        file.close();
      }
    }
  }
  
private:
  struct Sample {
    double time;
    unsigned numIters;
  };
  
  std::vector<Stopwatch> m_stopwatch;
  std::vector<std::string> m_regions;
  std::vector<std::vector<Sample>> m_times;
};

#endif // MONITORING_LOOPSTATISTICS_H_
