// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_PICKPOINTOUTPUT_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_PICKPOINTOUTPUT_H_

#include "Core/Geometry/MeshReader.h"
#include "Core/IO/Instance/Point/Hdf5Table.h"
#include "Core/IO/PickpointWriter.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Run/Context.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Output/Builders/PickPointBuilder.h"
#include "Features/DynamicRupture/Output/DataTypes.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Features/DynamicRupture/Output/ReceiverBasedOutput.h"

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::dr::output {

/// What the outputs of the fault take from the run once its layout is known.
struct FaultOutputSetup {
  //! the simulations of the run, and those fused into a cell
  std::size_t simulationCount{1};
  std::size_t fusedCount{1};
  //! the shortest time step of the run, which sizes the cache of the samples of on-fault receivers
  //! without a sampling interval
  double minimumTimestep{0};
  double endTime{0};
};

/**
 * One output of receivers on the fault (pickpoints): the samples of a set of points on the fault
 * at the multiples of an interval, into a text file per receiver or into an HDF5 table. The
 * clusters of the fault faces take the samples of their layers (see DynamicRuptureCluster), the
 * friction law evaluates them (ReceiverOutput).
 */
class PickpointOutput {
  public:
  /// The output `name`, with the settings `parameters`, evaluated by `impl`. The names of its files
  /// and of its table follow `name`; those of the output named "faultreceivers" keep the names
  /// SeisSol has always given them.
  PickpointOutput(std::string name,
                  initializer::parameters::PickpointParameters parameters,
                  ReceiverOutput& impl,
                  RunContext& run,
                  const initializer::parameters::DRParameters& drParameters);
  ~PickpointOutput();

  PickpointOutput(const PickpointOutput&) = delete;
  PickpointOutput(PickpointOutput&&) = delete;
  auto operator=(const PickpointOutput&) = delete;
  auto operator=(PickpointOutput&&) = delete;

  [[nodiscard]] const std::string& name() const { return name_; }

  /// The time between two samples: pickdt, or else printtimeinterval time steps of size
  /// `stepSize`, that of the finest time cluster with fault faces (for a run without them, as
  /// their defaults say).
  [[nodiscard]] double interval(double stepSize) const;

  /// Reads the receivers in `mesh`, for the run `setup`.
  void setInputParam(geometry::MeshReader& mesh, const FaultOutputSetup& setup);

  /// The storages the receivers lie in, and the map from fault faces to their faces in them.
  void setLtsData(LTS::Storage& wpStorage,
                  LTS::Backmap& wpBackmap,
                  DynamicRupture::Storage& drStorage,
                  initializer::StorageBackmap<1>* faceToLtsMap);

  /// Builds the receivers, and prepares their files (with backups stamped `backupTimeStamp`) or
  /// their table.
  void init(const std::string& backupTimeStamp);

  /// Records the receivers of all layers at the start of the simulation, `time`, before the first
  /// time step, on `runtime`; the time steps take the samples after it (see PickpointSampling).
  void recordInitial(double time, parallel::runtime::StreamRuntime& runtime);

  /// Whether the layer has receivers of this output.
  [[nodiscard]] bool hasReceivers(std::size_t layerId) const;

  /// The stream runtime the receivers of the layer take their samples on, if it is one of their
  /// own; null otherwise.
  [[nodiscard]] parallel::runtime::StreamRuntime* runtime(std::size_t layerId);

  /**
   * For samples that get decided when the work of a step runs: waits until the samples of the
   * last step have been recorded, takes the start of the step (from the clock on the device, into
   * pinned memory at `stepStartOnDevice` there; as `hostTime` into `stepStart` on the host), and
   * copies the data of the output points to the host. Returns the runtime to enqueue the decision
   * on.
   */
  parallel::runtime::StreamRuntime& gatherData(std::size_t layerId,
                                               const double* deviceClock,
                                               double* stepStart,
                                               double* stepStartOnDevice,
                                               double hostTime,
                                               parallel::runtime::StreamRuntime& runtime);

  /**
   * For samples that get decided when the work of a step runs, with the host actor (see
   * parallel::runtime::HostActor): gathers the data of the output points and the start of the
   * step (from the clock on the device) into a buffer of their own, and has the actor call
   * `record` with that start once they are there, also after each replay of a recording;
   * meanwhile, the data is at hand for evaluate() and finishStepNow(). Returns false, and does
   * nothing, without the host actor.
   */
  bool postAtRunTime(std::size_t layerId,
                     const double* deviceClock,
                     parallel::runtime::StreamRuntime& runtime,
                     std::function<void(double)> record);

  /// Records the receivers of the layer at `time` right away, from data gathered before; the
  /// times are those of record().
  void evaluate(std::size_t layerId, double time, double stepSize, double inStep);

  /// Records the receivers of the layer at `time`, in the time step of size `stepSize` that
  /// contains it, `inStep` after the start of that step.
  void record(std::size_t layerId,
              double time,
              double stepSize,
              double inStep,
              parallel::runtime::StreamRuntime& runtime);

  /// Keeps the friction state at the end of a time step of the layer, after its samples, as the
  /// start of the next one; `sampled` tells whether the step recorded any sample, which copied the
  /// state from the device already.
  void finishStep(std::size_t layerId, bool sampled, parallel::runtime::StreamRuntime& runtime);

  /// As finishStep(), right away, from data gathered before; for work that already runs in the
  /// order of the stream.
  void finishStepNow(std::size_t layerId);

  private:
  /// makes room in the cache for one more sample
  static void ensureCacheLevel(ReceiverOutputData& outputData);

  void initFiles(const std::string& backupTimeStamp);
  void initTable();

  /// Formats the samples cached since the last write (text output), which frees their cache for
  /// the next ones.
  void takeText();
  /// Appends what takeText() formatted to the files.
  void writeText();
  /// Moves the samples cached since the last write into the tables.
  void collectSamples();

  /// The stress the fault starts out under at the receiver `index`, aligned with strike and dip,
  /// as the six components of a symmetric tensor.
  std::array<double, 6> initialStress(const ReceiverOutputData& outputData, std::size_t index);

  std::string name_;
  initializer::parameters::PickpointParameters parameters_;
  ReceiverOutput& impl_;
  RunContext& run_;
  const initializer::parameters::DRParameters& drParameters_;
  FaultOutputSetup setup_;

  std::unique_ptr<PickPointBuilder> builder_;
  std::unordered_map<std::size_t, std::shared_ptr<ReceiverOutputData>> outputData_;

  struct File {
    std::string fileName;

    // the receiver the file holds: its receivers, one per simulation, in the order of their index
    std::vector<std::size_t> indices;
  };
  std::unordered_map<std::size_t, std::vector<File>> files_;

  //! one table per quantity set, holding the receivers this rank samples
  std::unique_ptr<io::instance::point::Hdf5Table> table_;
  //! what takeText() formatted, by file
  std::vector<std::pair<std::string, std::string>> takenText_;
  //! the layer and the receiver a row of the tables belongs to, in the order of the rows
  std::vector<std::pair<std::size_t, std::size_t>> tableRows_;

  DynamicRupture::Storage* drStorage_{nullptr};
  initializer::StorageBackmap<1>* faceToLtsMap_{nullptr};

  //! writes the text files at its synchronization points
  writer::PickpointWriter writer_;
};

} // namespace seissol::dr::output

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_PICKPOINTOUTPUT_H_
