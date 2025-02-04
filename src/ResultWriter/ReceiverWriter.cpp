/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2019, SeisSol Group
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
 **/

#include "ReceiverWriter.h"
#include "ParallelHdf5ReceiverWriter.h"

#include <Equations/Datastructures.h>
#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include <Initializer/Parameters/OutputParameters.h>
#include <Initializer/PointMapper.h>
#include <Initializer/Tree/Layer.h>
#include <Initializer/Tree/Lut.h>
#include <Kernels/Receiver.h>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <utils/logger.h>
#include <vector>

#include "Modules/Modules.h"
#include "Parallel/MPI.h"

namespace seissol::writer {

Eigen::Vector3d parseReceiverLine(const std::string& line) {
  const std::regex rgx("\\s+");
  std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
  const std::sregex_token_iterator end;
  Eigen::Vector3d coordinates{};
  unsigned numberOfCoordinates = 0;
  for (; iter != end; ++iter, ++numberOfCoordinates) {
    if (numberOfCoordinates >= coordinates.size()) {
      throw std::runtime_error("Too many coordinates in line " + line + ".");
    }
    coordinates[numberOfCoordinates] = std::stod(*iter);
  }
  if (numberOfCoordinates != coordinates.size()) {
    throw std::runtime_error("To few coordinates in line " + line + ".");
  }
  return coordinates;
}

std::vector<Eigen::Vector3d> parseReceiverFile(const std::string& receiverFileName) {
  std::vector<Eigen::Vector3d> points{};

  std::ifstream file{receiverFileName};
  std::string line{};
  while (std::getline(file, line)) {
    const bool onlyWhiteSpace = std::all_of(line.begin(), line.end(), [](auto& c) {
      return std::isspace(c); // This lambda is needed (opt. argument)
    });
    if (!onlyWhiteSpace) {
      points.emplace_back(parseReceiverLine(line));
    }
  }
  return points;
}

// --------------------------------------------------------------------------
// Instead of multiple ASCII files, we keep one HDF5 writer
static std::unique_ptr<ParallelHdf5ReceiverWriter> g_hdf5Writer = nullptr;

// We'll also keep track of how many time samples we've written so far
static hsize_t g_nextTimeOffset = 0;

// Keep track of how many receivers exist globally
static hsize_t g_totalReceivers = 0;

// We'll store each rank's offset in the "receiver" dimension
static hsize_t g_localReceiverOffset = 0;

// This is the new single-file name
static std::string hdf5FileName(const std::string& prefix) { return prefix + "_receivers.h5"; }

// --------------------------------------------------------------------------
void ReceiverWriter::init(
    const std::string& fileNamePrefix,
    double endTime,
    const seissol::initializer::parameters::ReceiverOutputParameters& parameters) {
  m_fileNamePrefix = fileNamePrefix;
  m_receiverFileName = parameters.fileName;
  m_samplingInterval = parameters.samplingInterval;

  if (parameters.computeRotation) {
    derivedQuantities.push_back(std::make_shared<kernels::ReceiverRotation>());
  }
  if (parameters.computeStrain) {
    derivedQuantities.push_back(std::make_shared<kernels::ReceiverStrain>());
  }

  setSyncInterval(std::min(endTime, parameters.interval));
  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  Modules::registerHook(*this, ModuleHook::Shutdown);

  // We do NOT open the HDF5 file yet, because we still need to know
  // total # of receivers, total # of columns, etc.
  // We'll do that in addPoints() after we parse them.
}

// --------------------------------------------------------------------------
void ReceiverWriter::addPoints(const seissol::geometry::MeshReader& mesh,
                               const seissol::initializer::Lut& ltsLut,
                               const seissol::initializer::LTS& lts,
                               const GlobalData* global) {
  // Same as before...
  std::vector<Eigen::Vector3d> points;
  const auto rank = seissol::MPI::mpi.rank();

  if (!m_receiverFileName.empty()) {
    points = parseReceiverFile(m_receiverFileName);
    logInfo(rank) << "Record points read from" << m_receiverFileName;
    logInfo(rank) << "Number of record points =" << points.size();
  } else {
    logInfo(rank) << "No record points read.";
  }

  const unsigned numberOfPoints = points.size();
  std::vector<short> contained(numberOfPoints);
  std::vector<unsigned> meshIds(numberOfPoints);

  // We want to plot all quantities except for the memory variables
  std::vector<unsigned> quantities(seissol::model::MaterialT::Quantities.size());
  std::iota(quantities.begin(), quantities.end(), 0);

  logInfo(rank) << "Finding meshIds for receivers...";
  initializer::findMeshIds(
      points.data(), mesh, numberOfPoints, contained.data(), meshIds.data(), 1e-3);
#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring receivers for MPI...";
  initializer::cleanDoubles(contained.data(), numberOfPoints);
#endif

  // Then reduce to see which points exist globally
  std::vector<short> globalContained(contained);
#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE,
                globalContained.data(),
                numberOfPoints,
                MPI_SHORT,
                MPI_MAX,
                seissol::MPI::mpi.comm());
#endif

  bool receiversMissing = false;
  for (int i = 0; i < numberOfPoints; ++i) {
    if (globalContained[i] == 0) {
      logWarning(rank) << "Receiver point" << i
                       << "could not be found. Coordinates:" << points[i](0) << points[i](1)
                       << points[i](2);
      receiversMissing = true;
    }
  }
  if (receiversMissing) {
    logError() << "Some receivers could not be found. Aborting simulation.";
  }

  logInfo(rank) << "Mapping receivers to LTS cells...";
  m_receiverClusters[Interior].clear();
  m_receiverClusters[Copy].clear();

  // Count how many we have locally
  size_t localReceiverCount = 0;

  for (unsigned point = 0; point < numberOfPoints; ++point) {
    if (contained[point] == 1) {
      const unsigned meshId = meshIds[point];
      const unsigned cluster = ltsLut.cluster(meshId);
      const LayerType layer = ltsLut.layer(meshId);

      auto& clusters = m_receiverClusters[layer];
      // Make sure that needed empty clusters are initialized.
      for (unsigned c = clusters.size(); c <= cluster; ++c) {
        clusters.emplace_back(global,
                              quantities,
                              m_samplingInterval,
                              syncInterval(),
                              derivedQuantities,
                              seissolInstance);
      }

      // For ASCII, we used to call writeHeader(point, points[point]) here,
      // but not needed for HDF5 (unless you want coordinate attributes).
      localReceiverCount++;

        m_receiverClusters[layer][cluster].addReceiver(
            meshId, point, points[point], mesh, ltsLut, lts);
    }
  }

  // -------------------------------------------------------
  // Now, sum up total # of receivers across ranks
  hsize_t localCountH = static_cast<hsize_t>(localReceiverCount);
  MPI_Allreduce(&localCountH,
                &g_totalReceivers,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                seissol::MPI::mpi.comm());

  logInfo(rank) << "Total number of receivers: " << g_totalReceivers;

  // We also need the offset in the "receivers" dimension for this rank
  MPI_Scan(&localCountH,
           &g_localReceiverOffset,
           1,
           MPI_UNSIGNED_LONG_LONG,
           MPI_SUM,
           seissol::MPI::mpi.comm());
  g_localReceiverOffset -= localCountH; // so the rank's chunk starts at "offset"

  // We can now open the single HDF5 file if not done yet:
  if (g_hdf5Writer == nullptr) {

    // Find the local maximum ncols
    unsigned localNcols = 0;
    for (auto& [layer, cvec] : m_receiverClusters) {
      for (auto& cluster : cvec) {
        // Instead of checking cluster.empty(), check if ncols() is nonzero
        if (cluster.ncols() > 0) {
          unsigned c = cluster.ncols();
          if (c > localNcols) {
            localNcols = c;
          }
        }
      }
    }

    // Gather the global maximum ncols
    unsigned globalNcols = 0;
    MPI_Allreduce(&localNcols, &globalNcols, 1, MPI_UNSIGNED, MPI_MAX, seissol::MPI::mpi.comm());

    // Now use globalNcols
    logInfo(rank) << "Global number of columns: " << globalNcols;

    // Create the HDF5 writer
    g_hdf5Writer = std::make_unique<ParallelHdf5ReceiverWriter>(seissol::MPI::mpi.comm(),
                                                                hdf5FileName(m_fileNamePrefix),
                                                                g_totalReceivers,
                                                                globalNcols);
  }

  g_hdf5Writer->writeCoordinates(points);
}

// --------------------------------------------------------------------------

void ReceiverWriter::syncPoint(double /*currentTime*/) {

  const auto rank = seissol::MPI::mpi.rank();

  size_t totalNewSamples = 0;
  size_t localReceiverCount = 0;

  for (auto& [layer, clusters] : m_receiverClusters) {
    for (auto& cluster : clusters) {
      for (auto& receiver : cluster) {
        size_t thisReceiverSamples = receiver.output.size() / cluster.ncols();
        if (thisReceiverSamples > totalNewSamples) {
          totalNewSamples = thisReceiverSamples;
        }
        localReceiverCount++;
      }
    }
  }

  bool noData = (localReceiverCount == 0 || totalNewSamples == 0);
  auto actualOffset = noData ? 0 : g_localReceiverOffset;
  auto actualTimeCount = noData ? 0 : totalNewSamples;
  auto actualRecCount = noData ? 0 : localReceiverCount;

  struct LocalReceiverData {
    kernels::Receiver* rcv;
    kernels::ReceiverCluster* clus;
  };
  std::vector<LocalReceiverData> localReceivers;
  localReceivers.reserve(localReceiverCount);

  std::vector<hsize_t> pointIds;
  pointIds.reserve(localReceiverCount);

  for (auto& [layer, clusters] : m_receiverClusters) {
    for (auto& cluster : clusters) {
      for (auto& receiver : cluster) {
        localReceivers.push_back({&receiver, &cluster});
        pointIds.push_back(receiver.pointId);
      }
    }
  }

  unsigned ncols = localReceivers.empty() ? 0 : localReceivers[0].clus->ncols();
  std::vector<double> hdf5Data(totalNewSamples * localReceiverCount * ncols);

  for (size_t lr = 0; lr < localReceivers.size(); ++lr) {
    auto& rec = *localReceivers[lr].rcv;
    auto& cluster = *localReceivers[lr].clus;
    const size_t nSamples = rec.output.size() / ncols;

    for (size_t t = 0; t < nSamples; ++t) {
      for (size_t v = 0; v < ncols; ++v) {
        double value = rec.output[t * ncols + v];
        size_t idx = t * (localReceiverCount * ncols) + lr * ncols + v;
        hdf5Data[idx] = value;
      }
    }
    rec.output.clear();
  }

  std::vector<double> emptyBuffer;
  g_hdf5Writer->writeChunk(g_nextTimeOffset,
                           actualTimeCount,
                           actualOffset,
                           actualRecCount,
                           noData ? emptyBuffer : hdf5Data);

  std::vector<unsigned long long> emptyPointIds;
  g_hdf5Writer->writePointIds(actualOffset,
                           actualRecCount,
                           noData ? emptyPointIds : pointIds);

  g_nextTimeOffset += totalNewSamples;
}

    // --------------------------------------------------------------------------
  void ReceiverWriter::simulationStart() {
    for (auto& [layer, clusters] : m_receiverClusters) {
      for (auto& cluster : clusters) {
        cluster.allocateData();
      }
    }
  }

  // --------------------------------------------------------------------------
  void ReceiverWriter::shutdown() {
    for (auto& [layer, clusters] : m_receiverClusters) {
      for (auto& cluster : clusters) {
        cluster.freeData();
      }
    }
    g_hdf5Writer.reset();
  }

  // --------------------------------------------------------------------------
  kernels::ReceiverCluster* ReceiverWriter::receiverCluster(unsigned clusterId, LayerType layer) {
    assert(layer != Ghost);
    assert(m_receiverClusters.find(layer) != m_receiverClusters.end());
    auto& clusters = m_receiverClusters[layer];
    if (clusterId < clusters.size()) {
      return &clusters[clusterId];
    }
    return nullptr;
  }

} // namespace seissol::writer
