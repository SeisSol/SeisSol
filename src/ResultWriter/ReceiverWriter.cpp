/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#include <cctype>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <Parallel/MPI.h>
#include <Modules/Modules.h>

#include <sstream>
#include <string>
#include <fstream>
#include <regex>

Eigen::Vector3d seissol::writer::parseReceiverLine(const std::string& line) {
  std::regex rgx("\\s+");
  std::sregex_token_iterator iter(line.begin(),
                                  line.end(),
                                  rgx,
                                  -1);
  std::sregex_token_iterator end;
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

std::vector<Eigen::Vector3d> seissol::writer::parseReceiverFile(const std::string& receiverFileName) {
  std::vector<Eigen::Vector3d> points{};

  std::ifstream file{receiverFileName};
  std::string line{};
  while (std::getline(file, line)) {
    bool onlyWhiteSpace = std::all_of(line.begin(), line.end(), [](auto& c) {
      return std::isspace(c); // This lambda is needed (opt. argument)
    });
    if (!onlyWhiteSpace) points.emplace_back(parseReceiverLine(line));
  }
  return points;
}

std::string seissol::writer::ReceiverWriter::fileName(unsigned pointId) const {
  std::stringstream fns;
  fns << std::setfill('0') << m_fileNamePrefix << "-receiver-" << std::setw(5) << (pointId+1);
#ifdef PARALLEL
  fns << "-" << std::setw(5) << seissol::MPI::mpi.rank();
#endif
  fns << ".dat";
  return fns.str();
}

void seissol::writer::ReceiverWriter::writeHeader( unsigned               pointId,
                                                   Eigen::Vector3d const& point   ) {
  auto name = fileName(pointId);

  std::vector<std::string> names({"xx", "yy", "zz", "xy", "yz", "xz", "u", "v", "w"});
#ifdef USE_POROELASTIC
  std::array<std::string, 4> additionalNames({"p", "u_f", "v_f", "w_f"});
  names.insert(names.end() ,additionalNames.begin(), additionalNames.end());
#endif

  /// \todo Find a nicer solution that is not so hard-coded.
  struct stat fileStat;
  // Write header if file does not exist
  if (stat(name.c_str(), &fileStat) != 0) {
    std::ofstream file;
    file.open(name);
    file << "TITLE = \"Temporal Signal for receiver number " << std::setfill('0') << std::setw(5) << (pointId+1) << "\"" << std::endl;
    file << "VARIABLES = \"Time\"";
#ifdef MULTIPLE_SIMULATIONS
    for (unsigned sim = init::QAtPoint::Start[0]; sim < init::QAtPoint::Stop[0]; ++sim) {
      for (auto const& name : names) {
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
      file << "# x" << (d+1) << "       " << std::scientific << std::setprecision(12) << point[d] << std::endl;
    }
    file.close();
  }
}

void seissol::writer::ReceiverWriter::syncPoint(double)
{
  if (m_receiverClusters.empty()) {
    return;
  }

  m_stopwatch.start();

  for (auto& [layer, clusters] : m_receiverClusters) {
    for (auto& cluster : clusters) {
      auto ncols = cluster.ncols();
      for (auto &receiver : cluster) {
        assert(receiver.output.size() % ncols == 0);
        size_t nSamples = receiver.output.size() / ncols;

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
  int const rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Wrote receivers in" << time << "seconds.";
}
void seissol::writer::ReceiverWriter::init(std::string receiverFileName, std::string fileNamePrefix,
                                           double syncPointInterval, double samplingInterval)
{
  m_receiverFileName = std::move(receiverFileName);
  m_fileNamePrefix = std::move(fileNamePrefix);
  m_samplingInterval = samplingInterval;
  setSyncInterval(syncPointInterval);
  Modules::registerHook(*this, SYNCHRONIZATION_POINT);
}

void seissol::writer::ReceiverWriter::addPoints(MeshReader const& mesh,
                                                const seissol::initializers::Lut& ltsLut,
                                                const seissol::initializers::LTS& lts,
                                                const GlobalData* global ) {
  std::vector<Eigen::Vector3d> points;
  const auto rank = seissol::MPI::mpi.rank();
  // Only parse if we have a receiver file
  if (!m_receiverFileName.empty()) {
    points = parseReceiverFile(m_receiverFileName);
    logInfo(rank) << "Record points read from" << m_receiverFileName;
    logInfo(rank) << "Number of record points =" << points.size();
  } else {
    logInfo(rank) << "No record points read.";
  }

  unsigned numberOfPoints = points.size();
  std::vector<short> contained(numberOfPoints);
  std::vector<unsigned> meshIds(numberOfPoints);
  
  // We want to plot all quantities except for the memory variables
  const int n = NUMBER_OF_QUANTITIES - 6*NUMBER_OF_RELAXATION_MECHANISMS;
  std::vector<unsigned> quantities(n);
  std::iota(quantities.begin(), quantities.end(), 0);

  logInfo(rank) << "Finding meshIds for receivers...";
  initializers::findMeshIds(points.data(), mesh, numberOfPoints, contained.data(), meshIds.data());
#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring receivers for MPI...";
  initializers::cleanDoubles(contained.data(), numberOfPoints);
#endif

  logInfo(rank) << "Mapping receivers to LTS cells...";
  m_receiverClusters[Interior].clear();
  m_receiverClusters[Copy].clear();
  for (unsigned point = 0; point < numberOfPoints; ++point) {
    if (contained[point] == 1) {
      unsigned meshId = meshIds[point];
      unsigned cluster = ltsLut.cluster(meshId);
      LayerType layer = ltsLut.layer(meshId);

      auto& clusters = m_receiverClusters[layer];
      // Make sure that needed empty clusters are initialized.
      for (unsigned c = clusters.size(); c <= cluster; ++c) {
        clusters.emplace_back(global, quantities, m_samplingInterval, syncInterval());
      }

      writeHeader(point, points[point]);
      m_receiverClusters[layer][cluster].addReceiver(meshId, point, points[point], mesh, ltsLut, lts);
    }
  }
}
