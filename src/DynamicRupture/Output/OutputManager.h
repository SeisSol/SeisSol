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
  void setBackupTimeStamp(const std::string& stamp) { this->backupTimeStamp = stamp; }

  void init();
  void initFaceToLtsMap();
  void writePickpointOutput(double time, double dt, parallel::runtime::StreamRuntime& runtime);
  void flushPickpointDataToFile();
  void updateElementwiseOutput();

  private:
  seissol::SeisSol& seissolInstance;

  protected:
  bool isAtPickpoint(double time, double dt);
  void initElementwiseOutput();
  void initPickpointOutput();

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder{nullptr};
  std::unique_ptr<PickPointBuilder> ppOutputBuilder{nullptr};

  std::shared_ptr<ReceiverOutputData> ewOutputData{nullptr};
  std::shared_ptr<ReceiverOutputData> ppOutputData{nullptr};

  struct PickpointFile {
    std::string fileName;

    // all receivers to be printed into this file
    std::vector<std::size_t> indices;
  };

  std::vector<PickpointFile> ppFiles;

  LTS::Storage* wpStorage{nullptr};
  LTS::Backmap* wpBackmap{nullptr};
  DynamicRupture::Storage* drStorage{nullptr};

  FaceToLtsMapType faceToLtsMap;
  std::vector<std::size_t> globalFaceToLtsMap;
  seissol::geometry::MeshReader* meshReader{nullptr};

  size_t iterationStep{0};
  static constexpr double TimeMargin{1.005};
  std::string backupTimeStamp;

  std::unique_ptr<ReceiverOutput> impl{nullptr};

  parallel::runtime::StreamRuntime runtime;
};
} // namespace dr::output
} // namespace seissol

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
