// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_

#include "DynamicRupture/Output/Builders/ElementWiseBuilder.h"
#include "DynamicRupture/Output/Builders/PickPointBuilder.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "IO/Instance/Point/Hdf5Table.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Tree/Backmap.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol {
class SeisSol;

namespace dr::output {

class OutputManager {
  public:
  ~OutputManager();
  auto operator=(const OutputManager&) = delete;
  auto operator=(OutputManager&&) = delete;
  OutputManager(const OutputManager&) = delete;
  OutputManager(OutputManager&&) = delete;

  OutputManager() = delete;
  OutputManager(std::unique_ptr<ReceiverOutput> concreteImpl, seissol::SeisSol& seissolInstance);
  void setInputParam(seissol::geometry::MeshReader& userMesher);
  void setLtsData(LTS::Storage& userWpStorage,
                  LTS::Backmap& userWpBackmap,
                  DynamicRupture::Storage& userDrStorage);
  void setBackupTimeStamp(const std::string& stamp) { this->backupTimeStamp_ = stamp; }

  void init();
  void initFaceToLtsMap();

  /**
   * Records the fault receivers of all layers at the start of the simulation, `time`, before the
   * first time step; the time steps take the samples after it (see PickpointSampling).
   */
  void recordInitialPickpointOutput(double time);

  /**
   * Whether the layer has fault receivers.
   */
  [[nodiscard]] bool hasPickpoints(std::size_t layerId) const;

  /**
   * The stream the fault receivers of the layer take their samples on, if it is one of their own;
   * null otherwise.
   */
  [[nodiscard]] void* pickpointStream(std::size_t layerId);

  /**
   * For fault receivers whose samples get decided when the work of a step runs: waits until
   * the fault receivers of the last step have been recorded, takes the start of the step (from the
   * clock on the device, into pinned memory at `stepStartOnDevice` there; as `hostTime` into
   * `stepStart` on the host), and copies the data of the output points to the host. Returns the
   * runtime to enqueue the decision on.
   */
  parallel::runtime::StreamRuntime& gatherPickpointData(std::size_t layerId,
                                                        const double* deviceClock,
                                                        double* stepStart,
                                                        double* stepStartOnDevice,
                                                        double hostTime,
                                                        parallel::runtime::StreamRuntime& runtime);

  /**
   * Records the fault receivers of the layer at `time` right away, from data gathered before; the
   * times are those of recordPickpointOutput().
   */
  void evaluatePickpointOutput(std::size_t layerId, double time, double stepSize, double inStep);

  /**
   * Records the fault receivers of the layer at `time`, in the time step of size `stepSize` that
   * contains it, `inStep` after the start of that step.
   */
  void recordPickpointOutput(std::size_t layerId,
                             double time,
                             double stepSize,
                             double inStep,
                             parallel::runtime::StreamRuntime& runtime);

  /**
   * Keeps the friction state at the end of a time step of the layer, after its samples, as the
   * start of the next one; `sampled` tells whether the step recorded any sample, which copied the
   * state from the device already.
   */
  void finishPickpointStep(std::size_t layerId,
                           bool sampled,
                           parallel::runtime::StreamRuntime& runtime);

  /**
   * As finishPickpointStep(), right away, from data gathered before; for work that already runs in
   * the order of the stream.
   */
  void finishPickpointStepNow(std::size_t layerId);
  /**
   * Formats the samples the fault receivers cached since the last write (text output), which frees
   * their cache for the next ones.
   */
  void takePickpointText();
  /** Appends what takePickpointText() formatted to the files. */
  void writePickpointText();
  /** Copies the data of the elementwise output from the device to the host, and waits for it. */
  void gatherElementwiseOutput();
  /** Evaluates the elementwise output from the data gathered, on the host. */
  void evaluateElementwiseOutput(double time);

  private:
  seissol::SeisSol& seissolInstance_;

  protected:
  /// makes room in the cache for one more sample
  static void ensureCacheLevel(ReceiverOutputData& outputData);
  void initElementwiseOutput();
  void initPickpointOutput();
  void initPickpointTable();

  //! @brief Moves the samples cached since the last write into the tables.
  void collectPickpointSamples();

  //! @brief The stress the fault starts out under at an on-fault receiver, aligned with strike
  //! and dip, as the six components of a symmetric tensor.
  std::array<double, 6> initialStress(const ReceiverOutputData& outputData, std::size_t index);

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder_{nullptr};
  std::unique_ptr<PickPointBuilder> ppOutputBuilder_{nullptr};

  std::shared_ptr<ReceiverOutputData> ewOutputData_{nullptr};
  std::unordered_map<std::size_t, std::shared_ptr<ReceiverOutputData>> ppOutputData_;

  struct PickpointFile {
    std::string fileName;

    // the receiver the file holds: its receivers, one per simulation, in the order of their index
    std::vector<std::size_t> indices;
  };

  std::unordered_map<std::size_t, std::vector<PickpointFile>> ppFiles_;

  //! One table per quantity set, holding the on-fault receivers this rank samples.
  std::unique_ptr<io::instance::point::Hdf5Table> ppTable_;
  //! what takePickpointText() formatted, by file
  std::vector<std::pair<std::string, std::string>> takenPickpointText_;
  //! The layer and the receiver a row of the tables belongs to, in the order of the rows.
  std::vector<std::pair<std::size_t, std::size_t>> ppTableRows_;

  LTS::Storage* wpStorage_{nullptr};
  LTS::Backmap* wpBackmap_{nullptr};
  DynamicRupture::Storage* drStorage_{nullptr};

  ::seissol::initializer::StorageBackmap<1> faceToLtsMap_;
  seissol::geometry::MeshReader* meshReader_{nullptr};

  std::string backupTimeStamp_;

  std::unique_ptr<ReceiverOutput> impl_{nullptr};

  parallel::runtime::StreamRuntime runtime_;
};
} // namespace dr::output
} // namespace seissol

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
