// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "ReceiverWriter.h"

#include <Equations/Datastructures.h>
#include <Geometry/MeshReader.h>
#include <Initializer/Parameters/OutputParameters.h>
#include <Initializer/PointMapper.h>
#include <Kernels/Receiver.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Memory/Tree/Lut.h>
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

std::string ReceiverWriter::fileName(unsigned pointId) const {
  std::stringstream fns;
  fns << std::setfill('0') << m_fileNamePrefix << "-receiver-" << std::setw(5) << (pointId + 1);
#ifdef PARALLEL
  fns << "-" << std::setw(5) << seissol::MPI::mpi.rank();
#endif
  fns << ".dat";
  return fns.str();
}

void ReceiverWriter::writeHeader(unsigned pointId, const Eigen::Vector3d& point) {
  auto name = fileName(pointId);

  std::vector<std::string> names(seissol::model::MaterialT::Quantities.begin(),
                                 seissol::model::MaterialT::Quantities.end());
  for (const auto& derived : derivedQuantities) {
    auto derivedNames = derived->quantities();
    names.insert(names.end(), derivedNames.begin(), derivedNames.end());
  }

  /// \todo Find a nicer solution that is not so hard-coded.
  struct stat fileStat;
  // Write header if file does not exist
  if (stat(name.c_str(), &fileStat) != 0) {
    std::ofstream file;
    file.open(name);
    file << "TITLE = \"Temporal Signal for receiver number " << std::setfill('0') << std::setw(5)
         << (pointId + 1) << "\"" << std::endl;
    file << "VARIABLES = \"Time\"";
#ifdef MULTIPLE_SIMULATIONS
    for (unsigned sim = init::QAtPoint::Start[0]; sim < init::QAtPoint::Stop[0]; ++sim) {
      for (const auto& name : names) {
        file << ",\"" << name << sim << "\"";
      }
    }
#else
    for (auto const& name : names) {
      file << ",\"" << name << "\"";
    }
#endif
    file << std::endl;
    for (int d = 0; d < 3; ++d) {
      file << "# x" << (d + 1) << "       " << std::scientific << std::setprecision(12) << point[d]
           << std::endl;
    }
    file.close();
  }
}

void ReceiverWriter::syncPoint(double /*currentTime*/) {
  if (m_receiverClusters.empty()) {
    return;
  }

  m_stopwatch.start();

  for (auto& [layer, clusters] : m_receiverClusters) {
    for (auto& cluster : clusters) {
      auto ncols = cluster.ncols();
      for (auto& receiver : cluster) {
        assert(receiver.output.size() % ncols == 0);
        const size_t nSamples = receiver.output.size() / ncols;

        std::ofstream file;
        file.open(fileName(receiver.pointId), std::ios::app);
        file << std::scientific << std::setprecision(15);
        for (size_t i = 0; i < nSamples; ++i) {
          for (size_t q = 0; q < ncols; ++q) {
            file << "  " << receiver.output[q + i * ncols];
          }
          file << std::endl;
        }
        file.close();
        receiver.output.clear();
      }
    }
  }

  auto time = m_stopwatch.stop();
  logInfo() << "Wrote receivers in" << time << "seconds.";
}
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
}

void ReceiverWriter::addPoints(const seissol::geometry::MeshReader& mesh,
                               const seissol::initializer::Lut& ltsLut,
                               const seissol::initializer::LTS& lts,
                               const GlobalData* global) {
  std::vector<Eigen::Vector3d> points;
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
  std::vector<unsigned> meshIds(numberOfPoints);

  // We want to plot all quantities except for the memory variables
  std::vector<unsigned> quantities(seissol::model::MaterialT::Quantities.size());
  std::iota(quantities.begin(), quantities.end(), 0);

  logInfo() << "Finding meshIds for receivers...";
  initializer::findMeshIds(
      points.data(), mesh, numberOfPoints, contained.data(), meshIds.data(), 1e-3);
  std::vector<short> globalContained(contained.begin(), contained.end());
#ifdef USE_MPI
  logInfo() << "Cleaning possible double occurring receivers for MPI...";
  initializer::cleanDoubles(contained.data(), numberOfPoints);
  MPI_Allreduce(MPI_IN_PLACE,
                globalContained.data(),
                globalContained.size(),
                MPI_SHORT,
                MPI_MAX,
                seissol::MPI::mpi.comm());
#endif

  bool receiversMissing = false;
  for (int i = 0; i < numberOfPoints; ++i) {
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
  m_receiverClusters[Interior].clear();
  m_receiverClusters[Copy].clear();
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

      writeHeader(point, points[point]);
      m_receiverClusters[layer][cluster].addReceiver(
          meshId, point, points[point], mesh, ltsLut, lts);
    }
  }
}

void ReceiverWriter::simulationStart() {
  for (auto& [layer, clusters] : m_receiverClusters) {
    for (auto& cluster : clusters) {
      cluster.allocateData();
    }
  }
}

void ReceiverWriter::shutdown() {
  for (auto& [layer, clusters] : m_receiverClusters) {
    for (auto& cluster : clusters) {
      cluster.freeData();
    }
  }
}

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
