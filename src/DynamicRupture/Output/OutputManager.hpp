#ifndef SEISSOL_DR_OUTPUT_MANAGER_HPP
#define SEISSOL_DR_OUTPUT_MANAGER_HPP

#include "DynamicRupture/Output/Builders/ElementWiseBuilder.hpp"
#include "DynamicRupture/Output/Builders/PickPointBuilder.hpp"
#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"
#include <memory>

namespace seissol::dr::output {
class OutputManager {
  public:
  ~OutputManager();
  OutputManager() = delete;
  OutputManager(ReceiverOutput* concreteImpl) : impl(concreteImpl){};
  void setInputParam(const YAML::Node& inputData, MeshReader& userMesher);
  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut,
                  seissol::initializers::LTSTree* userDrTree,
                  seissol::initializers::DynamicRupture* userDrDescr);

  void init();
  void initFaceToLtsMap();
  void writePickpointOutput(double time, double dt);
  void flushPickpointDataToFile();
  void updateElementwiseOutput();

  protected:
  bool isAtPickpoint(double time, double dt);
  void initElementwiseOutput();
  void initPickpointOutput();

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder{nullptr};
  std::unique_ptr<PickPointBuilder> ppOutputBuilder{nullptr};

  ReceiverOutputData ewOutputData{};
  ReceiverOutputData ppOutputData{};

  GeneralParams generalParams;
  ElementwiseFaultParams elementwiseParams{};
  PickpointParams pickpointParams{};

  seissol::initializers::LTS* wpDescr{nullptr};
  seissol::initializers::LTSTree* wpTree{nullptr};
  seissol::initializers::Lut* wpLut{nullptr};
  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* drDescr{nullptr};

  FaceToLtsMapType faceToLtsMap{};
  MeshReader* meshReader{nullptr};

  size_t iterationStep{0};
  static constexpr double timeMargin{1.005};

  std::unique_ptr<ReceiverOutput> impl{nullptr};
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_MANAGER_HPP
