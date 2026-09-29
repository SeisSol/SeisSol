// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: John Rekoske

#ifndef SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_

#include "Geometry/MeshReader.h"
#include "IO/Instance/Point/Hdf5Table.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Kernels/Receiver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"

#include <Eigen/Dense>
#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace seissol {
struct LocalIntegrationData;
class SeisSol;
} // namespace seissol

namespace seissol::writer {

Eigen::Vector3d parseReceiverLine(const std::string& line);
std::vector<Eigen::Vector3d> parseReceiverFile(const std::string& receiverFileName);

/**
 * \brief Writes out receiver data, either as one file per receiver or as a single HDF5 table.
 */
class ReceiverWriter : public seissol::Module {
  private:
  seissol::SeisSol& seissolInstance_;

  public:
  explicit ReceiverWriter(seissol::SeisSol& seissolInstance);
  ~ReceiverWriter() override;

  /**
   * \brief Initializes receiver output state and registers hooks required for recording.
   *
   * @param fileNamePrefix Prefix for the output file name including the folder it's stored in.
   * @param endTime        Final simulation time (used to estimate max # of steps).
   * @param parameters     Receiver output parameters (e.g., sampling interval).
   */
  void init(const std::string& fileNamePrefix,
            double endTime,
            const seissol::initializer::parameters::ReceiverOutputParameters& parameters);

  /**
   * \brief Registers receivers by reading their positions and mapping them to
   *        mesh cells, then creates/initializes the parallel HDF5 output.
   *
   *        This setup is required for receiver recording; without it, no
   *        receiver traces are produced.
   */
  void addPoints(const seissol::geometry::MeshReader& mesh,
                 const LTS::Backmap& backmap,
                 const CompoundGlobalData& global);

  /**
   * \brief Returns the ReceiverCluster for a given cluster ID.
   */
  kernels::ReceiverCluster* receiverCluster(std::size_t id);

  //
  // Hooks
  //
  /// Called at each synchronization point to flush data.
  void syncPoint(double currentTime) override;
  /// Called at simulation start.
  void simulationStart(std::optional<double> checkpointTime) override;

  /// Called at shutdown.
  void shutdown() override;

  private:
  [[nodiscard]] std::string fileName(std::size_t pointId) const;
  [[nodiscard]] std::vector<std::string> variableNames() const;
  void writeHeader(std::size_t pointId, const Eigen::Vector3d& point, std::size_t globalId);

  //! @brief A receiver together with the number of columns one of its samples takes.
  struct OrderedReceiver {
    kernels::Receiver* receiver{nullptr};
    std::size_t columns{0};
  };

  //! @brief The receivers this rank holds, ordered the way their rows are written.
  [[nodiscard]] std::vector<OrderedReceiver> orderedReceivers();

  //! @brief Moves the samples collected since the last write into the table.
  void collectSamples();

  // -- Members --
  seissol::initializer::parameters::ReceiverOutputFormat format_{
      seissol::initializer::parameters::ReceiverOutputFormat::Hdf5};
  std::string receiverFileName_;
  std::string fileNamePrefix_;
  double samplingInterval_{0.0};
  double endTime_{0.0};

  /// Additional derived quantities (e.g., rotation, strain)
  std::vector<std::shared_ptr<kernels::DerivedReceiverQuantity>> derivedQuantities_;

  std::vector<std::shared_ptr<kernels::ReceiverCluster>> receiverClusters_;

  /// One HDF5 table per quantity set, for the receivers this rank holds
  std::unique_ptr<io::instance::point::Hdf5Table> table_;

  /// How far a storage chunk of the table reaches along the sample axis
  std::size_t sampleChunk_{0};

  /// Stopwatch for timing receiver I/O only
  Stopwatch stopwatch_;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
