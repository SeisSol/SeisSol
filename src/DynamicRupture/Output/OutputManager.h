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
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Tree/Backmap.h"
#include "Parallel/Runtime/Stream.h"

#include <memory>
#include <unordered_map>

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
  void writePickpointOutput(double time, double dt);

  /**
   * Counts an output step of the layer; returns whether it records its fault receivers in it. Only
   * on the host, without touching any data.
   */
  bool beginPickpointStep(std::size_t layerId, double time, double dt);

  /**
   * Records the fault receivers of the layer, for an output step that
   * `beginPickpointStep()` has found due.
   */
  void recordPickpointOutput(std::size_t layerId,
                             double stateTime,
                             double time,
                             double meshDt,
                             double meshInDt,
                             parallel::runtime::StreamRuntime& runtime);
  void flushPickpointDataToFile();
  void updateElementwiseOutput(double time);

  private:
  seissol::SeisSol& seissolInstance_;

  protected:
  /**
   * Whether the layer records its fault receivers in its current step: in its first step, in every
   * `printTimeInterval`-th step after it, and close to the end of the simulation.
   */
  bool isAtPickpoint(std::size_t layerId, double time, double dt);
  void initElementwiseOutput();
  void initPickpointOutput();

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder_{nullptr};
  std::unique_ptr<PickPointBuilder> ppOutputBuilder_{nullptr};

  std::shared_ptr<ReceiverOutputData> ewOutputData_{nullptr};
  std::unordered_map<std::size_t, std::shared_ptr<ReceiverOutputData>> ppOutputData_;

  struct PickpointFile {
    std::string fileName;

    // all receivers to be printed into this file
    std::vector<std::size_t> indices;
  };

  std::unordered_map<std::size_t, std::vector<PickpointFile>> ppFiles_;

  LTS::Storage* wpStorage_{nullptr};
  LTS::Backmap* wpBackmap_{nullptr};
  DynamicRupture::Storage* drStorage_{nullptr};

  ::seissol::initializer::StorageBackmap<1> faceToLtsMap_;
  seissol::geometry::MeshReader* meshReader_{nullptr};

  //! the number of steps each layer has taken so far
  std::unordered_map<std::size_t, std::size_t> iterationSteps_;
  static constexpr double TimeMargin{1.005};
  std::string backupTimeStamp_;

  std::unique_ptr<ReceiverOutput> impl_{nullptr};

  parallel::runtime::StreamRuntime runtime_;
};
} // namespace dr::output
} // namespace seissol

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
