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
 
#include "LoopStatistics.h"

#include <cmath>
#ifdef USE_NETCDF
#include <netcdf.h>
#include <netcdf_par.h>
#endif

#include "Numerical_aux/Statistics.h"

#ifdef USE_MPI  
void seissol::LoopStatistics::printSummary(MPI_Comm comm) {
  unsigned const nRegions = m_times.size();
  auto sums = std::vector<double>(5*nRegions);
  double totalTimePerRank = 0.0;
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
    totalTimePerRank += y;
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  const auto summary = seissol::statistics::parallelSummary(totalTimePerRank);
  logInfo(rank) << "Time spent in compute kernels: mean =" << summary.mean
    << " std =" << summary.std
    << " min =" << summary.min
    << " median =" << summary.median
    << " max =" << summary.max;

  const auto loadImbalance = 1.0 - summary.mean / summary.max;
  logInfo(rank) << "Load imbalance:" << 100.0 * loadImbalance << "%";

  MPI_Allreduce(MPI_IN_PLACE, sums.data(), sums.size(), MPI_DOUBLE, MPI_SUM, comm);

  auto regressionCoeffs = std::vector<double>(2*nRegions);
  auto stderror = std::vector<double>(nRegions, 0.0);
  for (unsigned region = 0; region < nRegions; ++region) {
    double const x = sums[5*region + 0];
    double const x2 = sums[5*region + 1];
    double const xy = sums[5*region + 2];
    double const y = sums[5*region + 3];
    double const N = sums[5*region + 4];

    double const det = N*x2 - x*x;
    double const constant = (x2*y - x*xy) / det;
    double const slope = (-x*y + N*xy) / det;
    regressionCoeffs[2 * region + 0] = constant;
    regressionCoeffs[2 * region + 1] = slope;

    for (auto const& sample : m_times[region]) {
      if (sample.numIters > 0) {
        double const error = sample.time - (constant + slope * sample.numIters);
        stderror[region] += error * error;
      }
    }
  }

  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, stderror.data(), stderror.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  } else {
    MPI_Reduce(stderror.data(), 0L, stderror.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  }

  if (rank == 0) {
    double totalTime = 0.0;
    logInfo(rank) << "Regression analysis of compute kernels:";
    for (unsigned region = 0; region < nRegions; ++region) {
      double const x = sums[5*region + 0];
      double const x2 = sums[5*region + 1];
      double const y = sums[5*region + 3];
      double const N = sums[5*region + 4];

      double const xm = x / N;
      double const xv = x2 - 2*x*xm + xm*xm;

      // https://en.wikipedia.org/wiki/Simple_linear_regression#Normality_assumption
      double const se = std::sqrt((stderror[region] / (N-2)) / xv);

      char const* names[] = { "constant", "per element"};
      for (unsigned c = 0; c < 2; ++c) {
        logInfo(rank) << m_regions[region]
                      << "(" << names[c] << "):"
                      << regressionCoeffs[2 * region + c]
                      << "(sample size:" << N << ", standard error:" << se << ")";
      }
      totalTime += y;
    }

    logInfo(rank) << "Total time spent in compute kernels:" << totalTime;
  }
}
#endif

#ifdef USE_NETCDF
static void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    logError() << "line" << line << "of" << file << ":" << nc_strerror(stat) << std::endl;
  }
}
#endif
  
void seissol::LoopStatistics::writeSamples() {
  std::string loopStatFile = utils::Env::get<std::string>("SEISSOL_LOOP_STAT_PREFIX", "");
  if (!loopStatFile.empty()) {
#if defined(USE_NETCDF) && defined(USE_MPI)
    unsigned nRegions = m_times.size();
    for (unsigned region = 0; region < nRegions; ++region) {
      std::ofstream file;
      std::stringstream ss;
      ss << loopStatFile << m_regions[region] << ".nc";
      std::string fileName = ss.str();
      
      int nSamples = m_times[region].size();
      int sampleOffset;
      MPI_Scan(&nSamples, &sampleOffset, 1, MPI_INT, MPI_SUM, seissol::MPI::mpi.comm());
      
      int ncid, stat;
      stat = nc_create_par(fileName.c_str(), NC_MPIIO | NC_CLOBBER | NC_NETCDF4, seissol::MPI::mpi.comm(), MPI_INFO_NULL, &ncid); check_err(stat,__LINE__,__FILE__);
      
      int sampledim, rankdim, sampletyp, offsetid, sampleid;
      
      stat = nc_def_dim(ncid, "rank", 1+seissol::MPI::mpi.size(), &rankdim);             check_err(stat,__LINE__,__FILE__);
      stat = nc_def_dim(ncid, "sample", NC_UNLIMITED, &sampledim); check_err(stat,__LINE__,__FILE__);
      
      stat = nc_def_compound(ncid, sizeof(Sample), "Sample", &sampletyp); check_err(stat,__LINE__,__FILE__);
      {
        stat = nc_insert_compound(ncid, sampletyp, "time", NC_COMPOUND_OFFSET(Sample,time), NC_DOUBLE);   check_err(stat,__LINE__,__FILE__);
        stat = nc_insert_compound(ncid, sampletyp, "loopLength", NC_COMPOUND_OFFSET(Sample,numIters), NC_UINT); check_err(stat,__LINE__,__FILE__);
      }
      
      stat = nc_def_var(ncid, "offset", NC_INT,   1, &rankdim,   &offsetid); check_err(stat,__LINE__,__FILE__);
      stat = nc_def_var(ncid, "sample", sampletyp, 1, &sampledim, &sampleid); check_err(stat,__LINE__,__FILE__);
      
      stat = nc_enddef(ncid); check_err(stat,__LINE__,__FILE__);
      
      stat = nc_var_par_access(ncid, offsetid, NC_COLLECTIVE); check_err(stat,__LINE__,__FILE__);
      stat = nc_var_par_access(ncid, sampleid, NC_COLLECTIVE); check_err(stat,__LINE__,__FILE__);
  
      size_t start, count;
      int offsetData[2];
      if (seissol::MPI::mpi.rank() == 0) {
        start = 0;
        count = 2;        
        offsetData[0] = 0;
      } else {
        start = 1+seissol::MPI::mpi.rank();
        count = 1;
      }
      offsetData[count-1] = sampleOffset;
      stat = nc_put_vara_int(ncid, offsetid, &start, &count, offsetData);  check_err(stat,__LINE__,__FILE__);
      
      start = sampleOffset-nSamples;
      count = nSamples;
      stat = nc_put_vara(ncid, sampleid, &start, &count, m_times[region].data());  check_err(stat,__LINE__,__FILE__);      
      
      stat = nc_close(ncid); check_err(stat,__LINE__,__FILE__);
    }
#else
    logWarning(seissol::MPI::mpi.rank()) << "Writing loop statistics requires NetCDF and MPI.";
#endif
  }
}
