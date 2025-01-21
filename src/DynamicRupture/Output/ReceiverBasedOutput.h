// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_

#include "DynamicRupture/Output/ParametersInitializer.h"
#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Lut.h"

#include <DynamicRupture/Misc.h>
#include <memory>
#include <vector>

namespace seissol::dr::output {
class ReceiverOutput {
  public:
  virtual ~ReceiverOutput() = default;

  void setLtsData(seissol::initializer::LTSTree* userWpTree,
                  seissol::initializer::LTS* userWpDescr,
                  seissol::initializer::Lut* userWpLut,
                  seissol::initializer::LTSTree* userDrTree,
                  seissol::initializer::DynamicRupture* userDrDescr);

  void setMeshReader(seissol::geometry::MeshReader* userMeshReader) { meshReader = userMeshReader; }
  void setFaceToLtsMap(FaceToLtsMapType* map) { faceToLtsMap = map; }
  void calcFaultOutput(seissol::initializer::parameters::OutputType outputType,
                       seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
                       std::shared_ptr<ReceiverOutputData> outputData,
                       double time = 0.0);

  [[nodiscard]] virtual std::vector<std::size_t> getOutputVariables() const;

  protected:
  seissol::initializer::LTS* wpDescr{nullptr};
  seissol::initializer::LTSTree* wpTree{nullptr};
  seissol::initializer::Lut* wpLut{nullptr};
  seissol::initializer::LTSTree* drTree{nullptr};
  seissol::initializer::DynamicRupture* drDescr{nullptr};
  seissol::geometry::MeshReader* meshReader{nullptr};
  FaceToLtsMapType* faceToLtsMap{nullptr};
  real* deviceCopyMemory{nullptr};

  struct LocalInfo {
    seissol::initializer::Layer* layer{};
    size_t ltsId{};
    int nearestGpIndex{};
    int nearestInternalGpIndex{};

    std::size_t index{};

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

    ReceiverOutputData* state{};
  };

  template <typename T>
  std::remove_extent_t<T>* getCellData(const LocalInfo& local,
                                       const seissol::initializer::Variable<T>& variable) {
    auto devVar = local.state->deviceVariables.find(variable.index);
    if (devVar != local.state->deviceVariables.end()) {
      return reinterpret_cast<std::remove_extent_t<T>*>(
          devVar->second->get(local.state->deviceIndices[local.index]));
    } else {
      return local.layer->var(variable)[local.ltsId];
    }
  }

  void getDofs(real dofs[tensor::Q::size()], int meshId);
  void getNeighbourDofs(real dofs[tensor::Q::size()], int meshId, int side);
  void computeLocalStresses(LocalInfo& local);
  virtual real computeLocalStrength(LocalInfo& local) = 0;
  virtual real computeFluidPressure(LocalInfo& local) { return 0.0; }
  virtual real computeStateVariable(LocalInfo& local) { return 0.0; }
  static void updateLocalTractions(LocalInfo& local, real strength);
  real computeRuptureVelocity(Eigen::Matrix<real, 2, 2>& jacobiT2d, const LocalInfo& local);
  virtual void computeSlipRate(LocalInfo& local,
                               const std::array<real, 6>& /*rotatedUpdatedStress*/,
                               const std::array<real, 6>& /*rotatedStress*/);
  static void computeSlipRate(LocalInfo& local,
                              const std::array<double, 3>& tangent1,
                              const std::array<double, 3>& tangent2,
                              const std::array<double, 3>& strike,
                              const std::array<double, 3>& dip);
  virtual void outputSpecifics(std::shared_ptr<ReceiverOutputData>& data,
                               const LocalInfo& local,
                               size_t outputSpecifics,
                               size_t receiverIdx) {}
  virtual void adjustRotatedUpdatedStress(std::array<real, 6>& rotatedUpdatedStress,
                                          const std::array<real, 6>& rotatedStress) {};
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
