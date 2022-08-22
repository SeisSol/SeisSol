#ifndef SEISSOL_DR_RECEIVER_BASED_OUTPUT_HPP
#define SEISSOL_DR_RECEIVER_BASED_OUTPUT_HPP

#include "DynamicRupture/Output/ParametersInitializer.hpp"
#include "Geometry/MeshReader.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/LTS.h"
#include "Initializer/tree/Lut.hpp"

namespace seissol::dr::output {
class ReceiverBasedOutput {
  public:
  virtual ~ReceiverBasedOutput() = default;

  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut,
                  seissol::initializers::LTSTree* userDrTree,
                  seissol::initializers::DynamicRupture* userDrDescr);

  void setMeshReader(MeshReader* userMeshReader) { meshReader = userMeshReader; }
  void setFaceToLtsMap(FaceToLtsMapType* map) { faceToLtsMap = map; }
  void calcFaultOutput(OutputType type,
                       ReceiverOutputData& state,
                       const GeneralParams& generalParams,
                       double time = 0.0);

  protected:
  void getDofs(real dofs[tensor::Q::size()], int meshId);
  void getNeighbourDofs(real dofs[tensor::Q::size()], int meshId, int side);
  void computeLocalStresses();
  virtual real computeLocalStrength() = 0;
  virtual real computeFluidPressure() { return 0.0; }
  virtual real computeStateVariable() { return 0.0; }
  void updateLocalTractions(real strength);
  virtual void computeSlipRate(std::array<real, 6>&, std::array<real, 6>&);
  void computeSlipRate(const double* tangent1,
                       const double* tangent2,
                       const double* strike,
                       const double* dip);

  virtual void adjustRotatedUpdatedStress(std::array<real, 6>& rotatedUpdatedStress,
                                          const std::array<real, 6>& rotatedStress){};

  virtual void
      outputSpecifics(ReceiverOutputData& data, size_t outputSpecifics, size_t receiverIdx) {}
  real computeRuptureVelocity(Eigen::Matrix<real, 2, 2>& jacobiT2d);

  seissol::initializers::LTS* wpDescr{nullptr};
  seissol::initializers::LTSTree* wpTree{nullptr};
  seissol::initializers::Lut* wpLut{nullptr};
  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* drDescr{nullptr};
  MeshReader* meshReader{nullptr};
  FaceToLtsMapType* faceToLtsMap{nullptr};

  struct LocalInfo {
    seissol::initializers::Layer* layer{};
    size_t ltsId{};
    int nearestGpIndex{};
    int nearestInternalGpIndex{};

    real iniTraction1{};
    real iniTraction2{};

    real transientNormalTraction{};
    real iniNormalTraction{};
    real fluidPressure{};

    real frictionCoefficient{};
    real stateVariable{};

    real faultNormalVelocity{};

    real faceAlignedStress22{};
    real faceAlignedStress33{};
    real faceAlignedStress12{};
    real faceAlignedStress13{};
    real faceAlignedStress23{};

    real updatedTraction1{};
    real updatedTraction2{};

    real slipRateStrike{};
    real slipRateDip{};

    real faceAlignedValuesPlus[tensor::QAtPoint::size()]{};
    real faceAlignedValuesMinus[tensor::QAtPoint::size()]{};

    model::IsotropicWaveSpeeds* waveSpeedsPlus{};
    model::IsotropicWaveSpeeds* waveSpeedsMinus{};
  } local{};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_RECEIVER_BASED_OUTPUT_HPP
