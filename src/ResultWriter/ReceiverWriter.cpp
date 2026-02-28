// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "ReceiverWriter.h"

#include "Equations/Datastructures.h"
#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/PointMapper.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Receiver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Modules/Modules.h"
#include "Parallel/MPI.h"
#include "ParallelHdf5ReceiverWriter.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <optional>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <utils/logger.h>
#include <vector>

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
// Helper function for HDF5 output file name
std::string ReceiverWriter::hdf5FileName(const std::string& prefix) {
  return prefix + "_receivers.h5";
}

void ReceiverWriter::init(
    const std::string& fileNamePrefix,
    double endTime,
    const seissol::initializer::parameters::ReceiverOutputParameters& parameters) {
  m_fileNamePrefix = fileNamePrefix;
  m_receiverFileName = parameters.fileName;
  m_samplingInterval = parameters.samplingInterval;
  m_endTime = endTime;

  if (parameters.computeRotation) {
    derivedQuantities.push_back(std::make_shared<kernels::ReceiverRotation>());
  }
  if (parameters.computeStrain) {
    derivedQuantities.push_back(std::make_shared<kernels::ReceiverStrain>());
  }

  setSyncInterval(std::min(endTime, parameters.writeInterval));
  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  Modules::registerHook(*this, ModuleHook::Shutdown);
}

void ReceiverWriter::addPoints(const seissol::geometry::MeshReader& mesh,
                               const LTS::Backmap& backmap,
                               const CompoundGlobalData& global) {
  std::vector<Eigen::Vector3d> points;
  const auto rank = seissol::Mpi::mpi.rank();
  // Only parse if we have a receiver file
  if (!m_receiverFileName.empty()) {
    points = parseReceiverFile(m_receiverFileName);
    logInfo() << "Record points read from" << m_receiverFileName;
    logInfo() << "Number of record points =" << points.size();
  } else {
    logInfo() << "No record points read.";
  }

  const unsigned numberOfPoints = points.size();
  std::vector<short> contained(numberOfPoints);
  std::vector<std::size_t> meshIds(numberOfPoints);

  // We want to plot all quantities except for the memory variables
  std::vector<unsigned> quantities(seissol::model::MaterialT::Quantities.size());
  std::iota(quantities.begin(), quantities.end(), 0);

  logInfo() << "Finding meshIds for receivers...";
  initializer::findMeshIds(
      points.data(), mesh, numberOfPoints, contained.data(), meshIds.data(), 1e-3);

  logInfo() << "Cleaning possible double occurring receivers for multi-rank setups...";
  initializer::cleanDoubles(contained.data(), numberOfPoints);

  // Then reduce to see which points exist globally
  std::vector<short> globalContained(contained.begin(), contained.end());

  MPI_Allreduce(MPI_IN_PLACE,
                globalContained.data(),
                globalContained.size(),
                MPI_SHORT,
                MPI_MAX,
                seissol::Mpi::mpi.comm());

  bool receiversMissing = false;
  for (std::size_t i = 0; i < numberOfPoints; ++i) {
    if (globalContained[i] == 0) {
      logWarning() << "Receiver point" << i << "could not be found. Coordinates:" << points[i](0)
                   << points[i](1) << points[i](2);
      receiversMissing = true;
    }
  }
  if (receiversMissing) {
    logError() << "Some receivers could not be found. Aborting simulation.";
  }

  logInfo() << "Mapping receivers to LTS cells...";
  m_receiverClusters.clear();

  size_t localReceiverCount = 0;

  for (std::size_t point = 0; point < numberOfPoints; ++point) {
    if (contained[point] == 1) {
      const std::size_t meshId = meshIds[point];
      const auto id = backmap.get(meshId).color;

      // Make sure that needed empty clusters are initialized.
      for (std::size_t c = m_receiverClusters.size(); c <= id; ++c) {
        m_receiverClusters.emplace_back(
            std::make_shared<kernels::ReceiverCluster>(global,
                                                       quantities,
                                                       m_samplingInterval,
                                                       syncInterval(),
                                                       derivedQuantities,
                                                       seissolInstance));
      }

      // For ASCII, we used to call writeHeader(point, points[point]) here,
      // but not needed for HDF5 (unless you want coordinate attributes).
      localReceiverCount++;

      m_receiverClusters[id]->addReceiver(meshId, point, points[point], mesh, backmap);
    }
  }

  // -------------------------------------------------------
  // Now, sum up total # of receivers across ranks
  hsize_t localCountH = static_cast<hsize_t>(localReceiverCount);
  MPI_Allreduce(&localCountH,
                &m_totalReceivers,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                seissol::Mpi::mpi.comm());

  // We also need the offset in the "receivers" dimension for this rank
  MPI_Scan(&localCountH,
           &m_localReceiverOffset,
           1,
           MPI_UNSIGNED_LONG_LONG,
           MPI_SUM,
           seissol::Mpi::mpi.comm());
  m_localReceiverOffset -= localCountH; // so the rank's chunk starts at "offset"

  // We can now open the single HDF5 file if not done yet:
  if (m_hdf5Writer == nullptr) {

    // Find the local maximum ncols
    unsigned localNcols = 0;
    for (auto& cluster : m_receiverClusters) {
      if (cluster->ncols() > 0) {
        unsigned c = cluster->ncols();
        if (c > localNcols) {
          localNcols = c;
        }
      }
    }

    // Gather the global maximum ncols
    unsigned globalNcols = 0;
    MPI_Allreduce(&localNcols, &globalNcols, 1, MPI_UNSIGNED, MPI_MAX, seissol::Mpi::mpi.comm());

    // Create the HDF5 writer
    // Compute total number of time steps from simulation end time and sampling interval
    hsize_t totalTimeSteps = static_cast<hsize_t>(std::ceil(m_endTime / m_samplingInterval));

    m_hdf5Writer = std::make_unique<ParallelHdf5ReceiverWriter>(seissol::Mpi::mpi.comm(),
                                                                hdf5FileName(m_fileNamePrefix),
                                                                m_totalReceivers,
                                                                globalNcols,
                                                                totalTimeSteps);
  }

  m_hdf5Writer->writeCoordinates(points);
}

// --------------------------------------------------------------------------

void ReceiverWriter::syncPoint(double /*currentTime*/) {

  size_t totalNewSamples = 0;
  size_t localReceiverCount = 0;

  for (auto& cluster : m_receiverClusters) {
    for (auto& receiver : *cluster) {
      size_t thisReceiverSamples = receiver.output.size() / cluster->ncols();
      if (thisReceiverSamples > totalNewSamples) {
        totalNewSamples = thisReceiverSamples;
      }
      localReceiverCount++;
    }
  }

  bool noData = (localReceiverCount == 0 || totalNewSamples == 0);
  auto actualOffset = noData ? 0 : m_localReceiverOffset;
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

  for (auto& cluster : m_receiverClusters) {
    for (auto& receiver : *cluster) {
      localReceivers.push_back({&receiver, cluster.get()});
      pointIds.push_back(receiver.pointId);
    }
  }

  unsigned ncols = localReceivers.empty() ? 0 : localReceivers[0].clus->ncols();
  std::vector<double> hdf5Data(totalNewSamples * localReceiverCount * ncols);

  for (size_t lr = 0; lr < localReceivers.size(); ++lr) {
    auto& rec = *localReceivers[lr].rcv;
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
  m_hdf5Writer->writeChunk(m_nextTimeOffset,
                           actualTimeCount,
                           actualOffset,
                           actualRecCount,
                           noData ? emptyBuffer : hdf5Data);

  std::vector<hsize_t> emptyPointIds;
  m_hdf5Writer->writePointIds(actualOffset, actualRecCount, noData ? emptyPointIds : pointIds);

  m_nextTimeOffset += totalNewSamples;
}

// --------------------------------------------------------------------------
void ReceiverWriter::simulationStart(std::optional<double> checkpointTime) {
  for (auto& cluster : m_receiverClusters) {
    cluster->allocateData();
  }
}

// --------------------------------------------------------------------------
void ReceiverWriter::shutdown() {
  for (auto& cluster : m_receiverClusters) {
    cluster->freeData();
  }
  m_hdf5Writer.reset();
}

// --------------------------------------------------------------------------
kernels::ReceiverCluster* ReceiverWriter::receiverCluster(std::size_t id) {
  if (id < m_receiverClusters.size()) {
    return m_receiverClusters[id].get();
  }
  return nullptr;
}

} // namespace seissol::writer
