// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_

#include "DynamicRupture/Output/Builders/ElementWiseBuilder.h"
#include "DynamicRupture/Output/Builders/PickPointBuilder.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include <memory>

namespace seissol {
class SeisSol;

namespace dr::output {

class OutputManager {
  public:
  ~OutputManager();
  OutputManager() = delete;
  OutputManager(std::unique_ptr<ReceiverOutput> concreteImpl, seissol::SeisSol& seissolInstance);
  void setInputParam(seissol::geometry::MeshReader& userMesher);
  void setLtsData(seissol::initializer::LTSTree* userWpTree,
                  seissol::initializer::LTS* userWpDescr,
                  seissol::initializer::Lut* userWpLut,
                  seissol::initializer::LTSTree* userDrTree,
                  seissol::initializer::DynamicRupture* userDrDescr);
  void setBackupTimeStamp(const std::string& stamp) { this->backupTimeStamp = stamp; }

  void init();
  void initFaceToLtsMap();
  void writePickpointOutput(double time, double dt);
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

  seissol::initializer::LTS* wpDescr{nullptr};
  seissol::initializer::LTSTree* wpTree{nullptr};
  seissol::initializer::Lut* wpLut{nullptr};
  seissol::initializer::LTSTree* drTree{nullptr};
  seissol::initializer::DynamicRupture* drDescr{nullptr};

  FaceToLtsMapType faceToLtsMap{};
  std::vector<std::size_t> globalFaceToLtsMap;
  seissol::geometry::MeshReader* meshReader{nullptr};

  size_t iterationStep{0};
  static constexpr double timeMargin{1.005};
  std::string backupTimeStamp{};

  std::unique_ptr<ReceiverOutput> impl{nullptr};
};
} // namespace dr::output
} // namespace seissol

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTMANAGER_H_
