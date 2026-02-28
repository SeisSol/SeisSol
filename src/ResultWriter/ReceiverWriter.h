// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_

#include "Geometry/MeshReader.h"
#include "Kernels/Receiver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"

#include <Eigen/Dense>
#include <hdf5.h>
#include <memory>
#include <string_view>
#include <vector>

namespace seissol {
struct LocalIntegrationData;
struct GlobalData;
class SeisSol;
namespace initializer::parameters {
struct ReceiverOutputParameters; // forward declaration
} // namespace initializer::parameters
} // namespace seissol

class ParallelHdf5ReceiverWriter;

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
                 const LTS::Backmap& backmap,
                 const CompoundGlobalData& global);

  /**
   * \brief Returns the ReceiverCluster for a given cluster ID and layer type.
   */
  kernels::ReceiverCluster* receiverCluster(std::size_t id);

  //
  // Hooks
  //
  /// Called at each synchronization point (i.e., sampling) to flush data.
  void syncPoint(double currentTime) override;
  /// Called at simulation start.
  void simulationStart(std::optional<double> checkpointTime) override;

  /// Called at shutdown.
  void shutdown() override;

  private:
  static std::string hdf5FileName(const std::string& prefix);

  // -- Members --
  std::string m_receiverFileName;
  std::string m_fileNamePrefix;
  double m_samplingInterval{0.0};
  double m_endTime{0.0};

  /// Additional derived quantities (e.g., rotation, strain)
  std::vector<std::shared_ptr<kernels::DerivedReceiverQuantity>> derivedQuantities;

  std::vector<std::shared_ptr<kernels::ReceiverCluster>> m_receiverClusters;

  /// Parallel HDF5 writer for receiver data
  std::unique_ptr<ParallelHdf5ReceiverWriter> m_hdf5Writer;

  /// Current time offset for HDF5 writes
  hsize_t m_nextTimeOffset{0};

  /// Total number of receivers across all ranks
  hsize_t m_totalReceivers{0};

  /// This rank's offset in the receiver dimension
  hsize_t m_localReceiverOffset{0};

  /// Stopwatch for timing I/O
  Stopwatch m_stopwatch;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
