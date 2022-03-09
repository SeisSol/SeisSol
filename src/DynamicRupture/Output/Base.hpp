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

  void setInputParam(const YAML::Node& inputData, MeshReader& userMesher) {
    using namespace initializers;

    ParametersInitializer reader(inputData);
    mesher = &userMesher;
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
      ppOutputBuilder->setMeshReader(&userMesher);
      ppOutputBuilder->setParams(reader.getPickPointParams());
    }
    if (elementwiseEnabled) {
      logInfo() << "Enabling 2D fault output";
      ewOutputBuilder = std::make_unique<ElementWiseBuilder>();
      ewOutputBuilder->setMeshReader(&userMesher);
      ewOutputBuilder->setParams(reader.getElementwiseFaultParams());
    }
    if (!elementwiseEnabled && !pointEnabled) {
      logInfo() << "No dynamic rupture output enabled";
    }
  }

  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut,
                  seissol::initializers::LTSTree* userDrTree,
                  seissol::initializers::DynamicRupture* userDrDescr) {
    wpTree = userWpTree;
    wpDescr = userWpDescr;
    wpLut = userWpLut;
    drTree = userDrTree;
    drDescr = userDrDescr;
  }

  void init();
  void initFaceToLtsMap();

  void writePickpointOutput(double time, double dt);
  bool isAtPickpoint(double time, double dt);
  void updateElementwiseOutput();

  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* description,
                           seissol::Interoperability& eInteroperability);

  virtual void postCompute(seissol::initializers::DynamicRupture& drDescr) = 0;

  protected:
  void getDofs(real dofsPlus[tensor::Q::size()], int meshId, int side);
  void computeLocalStresses();
  virtual real computeLocalStrength() = 0;
  virtual real computePf() { return 0.0; }
  void computeLocalTraction(real strength);
  virtual void computeSlipAndRate(std::array<real, 6>&, std::array<real, 6>&);
  virtual void outputSpecifics(OutputData& data, size_t level, size_t receiverIdx) {}
  real computeRuptureVelocity();

  void initElementwiseOutput();
  void initPickpointOutput();

  [[nodiscard]] std::string constructPickpointReceiverFileName(int receiverGlobalIndex) const;
  void calcFaultOutput(OutputType type, OutputData& state, double time = 0.0);

  GeneralParamsT generalParams;

  std::unique_ptr<ElementWiseBuilder> ewOutputBuilder{nullptr};
  std::unique_ptr<PickPointBuilder> ppOutputBuilder{nullptr};

  seissol::initializers::LTS* wpDescr{nullptr};
  seissol::initializers::LTSTree* wpTree{nullptr};
  seissol::initializers::Lut* wpLut{nullptr};

  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* drDescr{nullptr};

  MeshReader* mesher{nullptr};

  std::vector<std::pair<seissol::initializers::Layer*, size_t>> faceToLtsMap{};
  size_t iterationStep{0};
  static constexpr double timeMargin{1.005};

  struct LocalInfo {
    seissol::initializers::Layer* layer{};
    size_t ltsId{};
    int nearestGpIndex{};

    real pf{};
    real mu{};
    real sXY{};
    real sXZ{};
    real p0{};

    real p{};
    real u{};
    real yyStress{};
    real zzStress{};
    real xyStress{};
    real xzStress{};
    real yzStress{};
    real tracEla{};

    real xyTraction{};
    real xzTraction{};

    real srS{};
    real srD{};

    real faceAlignedValuesPlus[tensor::QAtPoint::size()]{};
    real faceAlignedValuesMinus[tensor::QAtPoint::size()]{};

    model::IsotropicWaveSpeeds* waveSpeedsPlus{};
    model::IsotropicWaveSpeeds* waveSpeedsMinus{};
  } local{};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_BASE_HPP
