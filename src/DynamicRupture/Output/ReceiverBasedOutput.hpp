#ifndef SEISSOL_DR_RECEIVER_BASED_OUTPUT_HPP
#define SEISSOL_DR_RECEIVER_BASED_OUTPUT_HPP

#include "DynamicRupture/Output/ParametersInitializer.hpp"
#include "Geometry/MeshReader.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/LTS.h"
#include "Initializer/tree/Lut.hpp"

#include <memory>
#include <vector>

namespace seissol::dr::output {
class ReceiverOutput {
  public:
  virtual ~ReceiverOutput();

  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut,
                  seissol::initializers::LTSTree* userDrTree,
                  seissol::initializers::DynamicRupture* userDrDescr);

  void allocateMemory(const std::vector<std::shared_ptr<ReceiverOutputData>>& states);
  void setMeshReader(seissol::geometry::MeshReader* userMeshReader) { meshReader = userMeshReader; }
  void setFaceToLtsMap(FaceToLtsMapType* map) { faceToLtsMap = map; }
  void calcFaultOutput(OutputType type,
                       std::shared_ptr<ReceiverOutputData> state,
                       const GeneralParams& generalParams,
                       double time = 0.0);

  protected:
  seissol::initializers::LTS* wpDescr{nullptr};
  seissol::initializers::LTSTree* wpTree{nullptr};
  seissol::initializers::Lut* wpLut{nullptr};
  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* drDescr{nullptr};
  seissol::geometry::MeshReader* meshReader{nullptr};
  FaceToLtsMapType* faceToLtsMap{nullptr};
  real* deviceCopyMemory{nullptr};

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
  };

  void getDofs(real dofs[tensor::Q::size()], int meshId);
  void getNeighbourDofs(real dofs[tensor::Q::size()], int meshId, int side);
  void computeLocalStresses(LocalInfo& local);
  virtual real computeLocalStrength(LocalInfo& local) = 0;
  virtual real computeFluidPressure(LocalInfo& local) { return 0.0; }
  virtual real computeStateVariable(LocalInfo& local) { return 0.0; }
  void updateLocalTractions(LocalInfo& local, real strength);
  real computeRuptureVelocity(Eigen::Matrix<real, 2, 2>& jacobiT2d, const LocalInfo& local);
  virtual void
      computeSlipRate(LocalInfo& local, const std::array<real, 6>&, const std::array<real, 6>&);
  void computeSlipRate(LocalInfo& local,
                       const std::array<double, 3>& tangent1,
                       const std::array<double, 3>& tangent2,
                       const std::array<double, 3>& strike,
                       const std::array<double, 3>& dip);
  virtual void outputSpecifics(std::shared_ptr<ReceiverOutputData>& data,
                               const LocalInfo& local,
                               size_t outputSpecifics,
                               size_t receiverIdx) {}
  virtual void adjustRotatedUpdatedStress(std::array<real, 6>& rotatedUpdatedStress,
                                          const std::array<real, 6>& rotatedStress){};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_RECEIVER_BASED_OUTPUT_HPP
