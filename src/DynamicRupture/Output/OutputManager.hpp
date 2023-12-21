#ifndef SEISSOL_DR_OUTPUT_MANAGER_HPP
#define SEISSOL_DR_OUTPUT_MANAGER_HPP

#include "DynamicRupture/Output/Builders/ElementWiseBuilder.hpp"
#include "DynamicRupture/Output/Builders/PickPointBuilder.hpp"
#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"
#include <memory>

namespace seissol {
class SeisSol;

namespace dr::output {

class OutputManager {
  public:
  ~OutputManager();
  OutputManager() = delete;
  OutputManager(std::unique_ptr<ReceiverOutput> concreteImpl, seissol::SeisSol& seissolInstance);
  void setInputParam(const YAML::Node& inputData, seissol::geometry::MeshReader& userMesher);
  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut,
                  seissol::initializers::LTSTree* userDrTree,
                  seissol::initializers::DynamicRupture* userDrDescr);
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

  GeneralParams generalParams;
  ElementwiseFaultParams elementwiseParams{};
  PickpointParams pickpointParams{};

  seissol::initializers::LTS* wpDescr{nullptr};
  seissol::initializers::LTSTree* wpTree{nullptr};
  seissol::initializers::Lut* wpLut{nullptr};
  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* drDescr{nullptr};

  FaceToLtsMapType faceToLtsMap{};
  seissol::geometry::MeshReader* meshReader{nullptr};

  size_t iterationStep{0};
  static constexpr double timeMargin{1.005};
  std::string backupTimeStamp{};

  std::unique_ptr<ReceiverOutput> impl{nullptr};
};
} // namespace dr::output
} // namespace seissol

#endif // SEISSOL_DR_OUTPUT_MANAGER_HPP
