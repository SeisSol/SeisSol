// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_

#include "Core/Memory/Tree/Backmap.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Features/DynamicRupture/Output/Builders/ElementWiseBuilder.h"
#include "Features/DynamicRupture/Output/DataTypes.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Features/DynamicRupture/Output/PickpointOutput.h"
#include "Features/DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Simulation/Parameters/SeisSolParameters.h"

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
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
  /// The settings of the outputs of the fault: of its whole surface (elementwise), none if the
  /// run does not write it, and of each output of receivers on it (pickpoints), by its name.
  void
      setOutputs(const seissol::initializer::parameters::ElementwiseFaultParameters* elementwise,
                 const std::vector<
                     std::pair<std::string, seissol::initializer::parameters::PickpointParameters>>&
                     pickpoints);
  void setInputParam(seissol::geometry::MeshReader& userMesher);
  void setLtsData(LTS::Storage& userWpStorage,
                  LTS::Backmap& userWpBackmap,
                  DynamicRupture::Storage& userDrStorage);
  void setBackupTimeStamp(const std::string& stamp) { this->backupTimeStamp_ = stamp; }

  void init();
  void initFaceToLtsMap();

  /// The outputs of receivers on the fault, in the order of the parameters.
  [[nodiscard]] std::vector<PickpointOutput*> pickpointOutputs();

  /**
   * Records the receivers of all outputs on the fault at the start of the simulation, `time`,
   * before the first time step; the time steps take the samples after it (see PickpointSampling).
   */
  void recordInitialPickpointOutput(double time);

  /** Copies the data of the elementwise output from the device to the host, and waits for it. */
  void gatherElementwiseOutput();
  /** Evaluates the elementwise output from the data gathered, on the host. */
  void evaluateElementwiseOutput(double time);

  private:
  seissol::SeisSol& seissolInstance_;

  protected:
  void initElementwiseOutput();

  //! the settings of the elementwise output, none if the run does not write it
  std::optional<seissol::initializer::parameters::ElementwiseFaultParameters> elementwise_;

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder_{nullptr};
  std::shared_ptr<ReceiverOutputData> ewOutputData_{nullptr};

  //! the outputs of receivers on the fault
  std::vector<std::unique_ptr<PickpointOutput>> pickpointOutputs_;

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

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
