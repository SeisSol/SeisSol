#ifndef SEISSOL_DROUTOUT_DR_BASE_HPP
#define SEISSOL_DROUTOUT_DR_BASE_HPP

#include "Initializer/InputAux.hpp"
#include "Initializer/DynamicRupture.h"
#include "DynamicRupture/Output/ParametersInitializer.hpp"
#include "DynamicRupture/Output/Builders/ElementWiseOutput.hpp"
#include "DynamicRupture/Output/Builders/PickpointOutput.hpp"
#include <iostream>
#include <memory>

namespace seissol::dr::output {
class Base {
  public:
  virtual ~Base() = default;

  void setInputParam(const YAML::Node& InputData, const MeshReader& Mesher) {
    using namespace initializers;

    ParametersInitializer Reader(InputData);
    generalParams = Reader.getDrGeneralParams();

    // adjust general output parameters
    generalParams.isRfTimeOn = generalParams.isRfOutputOn;
    if (generalParams.isDsOutputOn && !generalParams.isRfOutputOn) {
      generalParams.isRfOutputOn = true;
      generalParams.isRfTimeOn = true;
    }

    PickpointParamsT PpParams;
    ElementwiseFaultParamsT EwParams;
    switch (generalParams.outputPointType) {
    case OutputType::None:
      break;

    case OutputType::AtPickpoint:
      ppOutputBuilder = std::make_unique<PickpointOutput>();
      ppOutputBuilder->setParams(Reader.getPickPointParams(), &Mesher);
      break;

    case OutputType::Elementwise:
      ewOutputBuilder = std::make_unique<ElementWiseOutput>();
      ewOutputBuilder->setParams(Reader.getElementwiseFaultParams(), &Mesher);
      break;

    case OutputType::AtPickpointAndElementwise:
      ppOutputBuilder = std::make_unique<PickpointOutput>();
      ppOutputBuilder->setParams(Reader.getPickPointParams(), &Mesher);

      ewOutputBuilder = std::make_unique<ElementWiseOutput>();
      ewOutputBuilder->setParams(Reader.getElementwiseFaultParams(), &Mesher);
      break;

    default:
      throw std::runtime_error("Unknown fault output type (not 3,4,5)");
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
                           seissol::Interoperability& e_interoperability);

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) = 0;

  protected:
  void initEwOutput();
  void initPickpointOutput();

  [[nodiscard]] std::string constructPpReceiverFileName(int receiverGlobalIndex) const;
  void calcFaultOutput(OutputType type, OutputData& state, double time = 0.0);

  GeneralParamsT generalParams;

  std::unique_ptr<ElementWiseOutput> ewOutputBuilder{nullptr};
  OutputData ewOutputData{};
  // std::vector<std::pair<initializers::Layer*, size_t>> faceToLtsMap{};

  std::unique_ptr<PickpointOutput> ppOutputBuilder{nullptr};
  OutputData ppOutputState{};

  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* dynRup{nullptr};

  std::vector<std::pair<seissol::initializers::Layer*, size_t>> faceToLtsMap{};
  size_t iterationStep{0};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DROUTOUT_DR_BASE_HPP
