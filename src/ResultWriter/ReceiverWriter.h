/**
 * @file
 * This file is part of SeisSol.
 *
 * @author
 *    Carsten Uphoff (c.uphoff AT tum.de,
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
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
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

#ifndef RESULTWRITER_RECEIVERWRITER_H_
#define RESULTWRITER_RECEIVERWRITER_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "Geometry/MeshReader.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/Lut.h"
#include "Kernels/Receiver.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"

namespace seissol {
struct LocalIntegrationData;
struct GlobalData;
class SeisSol;
namespace initializer::parameters {
struct ReceiverOutputParameters; // forward declaration
} // namespace initializer::parameters
} // namespace seissol

namespace seissol::writer {

// Forward declarations for helper functions
Eigen::Vector3d parseReceiverLine(const std::string& line);
std::vector<Eigen::Vector3d> parseReceiverFile(const std::string& receiverFileName);

/**
 * \brief Writes out receiver data in a single parallel HDF5 file
 * instead of multiple ASCII files.
 */
class ReceiverWriter : public seissol::Module {
  private:
  seissol::SeisSol& seissolInstance;

  public:
  explicit ReceiverWriter(seissol::SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

  /**
   * \brief Initializes internal data.
   *
   * @param fileNamePrefix Prefix for the output file name.
   * @param endTime        Final simulation time (used to estimate max # of steps).
   * @param parameters     Receiver output parameters (e.g., sampling interval).
   */
  void init(const std::string& fileNamePrefix,
            double endTime,
            const seissol::initializer::parameters::ReceiverOutputParameters& parameters);

  /**
   * \brief Registers receivers by reading their positions and mapping them to
   *        mesh cells, then (optionally) creates the single parallel HDF5 file.
   */
  void addPoints(const seissol::geometry::MeshReader& mesh,
                 const seissol::initializer::Lut& ltsLut,
                 const seissol::initializer::LTS& lts,
                 const seissol::GlobalData* global);

  /**
   * \brief Returns the ReceiverCluster for a given cluster ID and layer type.
   */
  kernels::ReceiverCluster* receiverCluster(unsigned clusterId, LayerType layer);

  //
  // Hooks
  //
  /// Called at each synchronization point (i.e., sampling) to flush data.
  void syncPoint(double currentTime) override;
  /// Called at simulation start.
  void simulationStart() override;
  /// Called at shutdown.
  void shutdown() override;

  private:
  /** Not used for HDF5, but remains for legacy or potential metadata. */
  [[nodiscard]] std::string fileName(unsigned pointId) const;

  /** Not used for HDF5, but remains for potential coordinate attributes, etc. */
  void writeHeader(unsigned pointId, const Eigen::Vector3d& point);

  // -- Members --
  std::string m_receiverFileName;
  std::string m_fileNamePrefix;
  double m_samplingInterval{0.0};

  /// Additional derived quantities (e.g., rotation, strain)
  std::vector<std::shared_ptr<kernels::DerivedReceiverQuantity>> derivedQuantities;

  /// Map layer type -> array of ReceiverClusters
  std::unordered_map<LayerType, std::vector<kernels::ReceiverCluster>> m_receiverClusters;

  /// Stopwatch for timing I/O
  Stopwatch m_stopwatch;
};

} // namespace seissol::writer

#endif // RESULTWRITER_RECEIVERWRITER_H_
