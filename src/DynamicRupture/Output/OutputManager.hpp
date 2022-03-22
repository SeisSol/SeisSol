#ifndef SEISSOL_DR_OUTPUT_MANAGER_HPP
#define SEISSOL_DR_OUTPUT_MANAGER_HPP

#include "DynamicRupture/Output/Builders/ElementWiseBuilder.hpp"
#include "DynamicRupture/Output/Builders/PickPointBuilder.hpp"
#include "DynamicRupture/Output/Base.hpp"
#include <iostream>
#include <memory>

namespace seissol::dr::output {
class Base;

class OutputManager {
  public:
  ~OutputManager();
  OutputManager() = delete;
  OutputManager(Base* concreteImpl) : impl(concreteImpl){};

  void setInputParam(const YAML::Node& inputData, MeshReader& userMesher);

  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut,
                  seissol::initializers::LTSTree* userDrTree,
                  seissol::initializers::DynamicRupture* userDrDescr);

  void init();
  void initFaceToLtsMap();

  void writePickpointOutput(double time, double dt);
  bool isAtPickpoint(double time, double dt);
  void updateElementwiseOutput();

  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* description,
                   seissol::Interoperability& eInteroperability);

  protected:
  void initElementwiseOutput();
  void initPickpointOutput();

  [[nodiscard]] std::string constructPickpointReceiverFileName(int receiverGlobalIndex) const;

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder{nullptr};
  std::unique_ptr<PickPointBuilder> ppOutputBuilder{nullptr};

  OutputData ewOutputData{};
  OutputData ppOutputData{};

  GeneralParamsT generalParams;
  ElementwiseFaultParamsT elementwiseParams{};
  PickpointParamsT pickpointParams{};

  size_t iterationStep{0};
  static constexpr double timeMargin{1.005};

  std::unique_ptr<Base> impl{nullptr};
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_MANAGER_HPP
