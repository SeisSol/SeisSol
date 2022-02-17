#ifndef SEISSOL_DR_OUTPUT_BASE_HPP
#define SEISSOL_DR_OUTPUT_BASE_HPP

#include "DynamicRupture/Output/Builders/ElementWiseBuilder.hpp"
#include "DynamicRupture/Output/Builders/PickPointBuilder.hpp"
#include "DynamicRupture/Output/ParametersInitializer.hpp"
#include "Initializer/DynamicRupture.h"
#include "Initializer/InputAux.hpp"
#include <iostream>
#include <memory>

namespace seissol::dr::output {
class Base {
  public:
  virtual ~Base() = default;

  void setInputParam(const YAML::Node& inputData, const MeshReader& mesher) {
    using namespace initializers;

    ParametersInitializer reader(inputData);
    generalParams = reader.getDrGeneralParams();

    // adjust general output parameters
    generalParams.isRfTimeOn = generalParams.isRfOutputOn;
    if (generalParams.isDsOutputOn && !generalParams.isRfOutputOn) {
      generalParams.isRfOutputOn = true;
      generalParams.isRfTimeOn = true;
    }

    bool bothEnabled = generalParams.outputPointType == OutputType::AtPickpointAndElementwise;
    bool pointEnabled = generalParams.outputPointType == OutputType::AtPickpoint || bothEnabled;
    bool elementwiseEnabled =
        generalParams.outputPointType == OutputType::Elementwise || bothEnabled;
    if (pointEnabled) {
      logInfo() << "Enabling on-fault receiver output";
      ppOutputBuilder = std::make_unique<PickPointBuilder>();
      ppOutputBuilder->setMeshReader(&mesher);
      ppOutputBuilder->setParams(reader.getPickPointParams());
    }
    if (elementwiseEnabled) {
      logInfo() << "Enabling 2D fault output";
      ewOutputBuilder = std::make_unique<ElementWiseBuilder>();
      ewOutputBuilder->setMeshReader(&mesher);
      ewOutputBuilder->setParams(reader.getElementwiseFaultParams());
    }
    if (!elementwiseEnabled && !pointEnabled) {
      logInfo() << "No dynamic rupture output enabled";
    }
  }
  void setDrData(seissol::initializers::LTSTree* userDrTree,
                 seissol::initializers::DynamicRupture* drDescription) {
    drTree = userDrTree;
    dynRup = drDescription;
  }

  void init();
  void initFaceToLtsMap();

  void writePickpointOutput(double time, double dt);
  bool isAtPickpoint(double time, double dt);
  void updateElementwiseOutput();

  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* description,
                           seissol::Interoperability& eInteroperability);

  virtual void postCompute(seissol::initializers::DynamicRupture& dynRup) = 0;

  protected:
  void initElementwiseOutput();
  void initPickpointOutput();

  [[nodiscard]] std::string constructPickpointReceiverFileName(int receiverGlobalIndex) const;
  void calcFaultOutput(OutputType type, OutputData& state, double time = 0.0);

  GeneralParamsT generalParams;

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder{nullptr};
  std::unique_ptr<PickPointBuilder> ppOutputBuilder{nullptr};

  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* dynRup{nullptr};

  std::vector<std::pair<seissol::initializers::Layer*, size_t>> faceToLtsMap{};
  size_t iterationStep{0};

  static constexpr double timeMargin{1.005};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_BASE_HPP
