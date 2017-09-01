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
    double* regressionCoeffs = new double[2*nRegions];
    for (unsigned region = 0; region < nRegions; ++region) {
      double x = 0.0, x2 = 0.0, xy = 0.0, y = 0.0;
      unsigned N = m_times[region].size();
    
      for (auto const& sample : m_times[region]) {
        x  += sample.numIters;
        x2 += sample.numIters * sample.numIters;
        xy += sample.numIters * sample.time;
        y  += sample.time;
      }
      double det = N*x2 - x*x;
      if (det != 0.0) {
        regressionCoeffs[2*region + 0] = (x2*y - x*xy) / det;
        regressionCoeffs[2*region + 1] = (-x*y + N*xy) / det;
      } else {
        regressionCoeffs[2*region + 0] = NAN;
        regressionCoeffs[2*region + 1] = NAN;
      }
    }
    
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double* regressionCoeffsWorld = nullptr;
    if (rank == 0) {
      regressionCoeffsWorld = new double[2*nRegions*size];
    }
    
    MPI_Gather(regressionCoeffs, 2*nRegions, MPI_DOUBLE, regressionCoeffsWorld, 2*nRegions, MPI_DOUBLE, 0, comm);    
    delete[] regressionCoeffs;
    
    if (rank == 0) {
      char const* names[] = { "constant", "per element"};
      logInfo(rank) << "Regression analysis of compute kernels:";
      for (unsigned region = 0; region < nRegions; ++region) {
        for (unsigned c = 0; c < 2; ++c) {
          double avg = 0.0, min = std::numeric_limits<double>::infinity(), max = -std::numeric_limits<double>::infinity();
          unsigned nValid = 0;
          for (int rk = 0; rk < size; ++rk) {
            double val = regressionCoeffsWorld[2*nRegions*rk + 2*nRegions*region + c];
            if (!std::isnan(val)) {
              avg += val;
              min = std::min(min, val);
              max = std::max(max, val);
              ++nValid;
            }
          }
          avg /= nValid;
          
          logInfo(rank) << utils::nospace
                        << m_regions[region]
                        << " (" << names[c] << "): "
                        << avg
                        << " (min:"
                        << min
                        << ", max: "
                        << max
                        << ')';
        }
      }
      delete[] regressionCoeffsWorld;
    }
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
