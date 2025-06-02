// SPDX-FileCopyrightText: 2022 SeisSol Group
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

    real iniTraction1[seissol::multisim::NumSimulations]{};
    real iniTraction2[seissol::multisim::NumSimulations]{};

    real transientNormalTraction[seissol::multisim::NumSimulations]{};
    real iniNormalTraction[seissol::multisim::NumSimulations]{};
    real fluidPressure[seissol::multisim::NumSimulations]{};

    real frictionCoefficient[seissol::multisim::NumSimulations]{};
    real stateVariable[seissol::multisim::NumSimulations]{};

    real faultNormalVelocity[seissol::multisim::NumSimulations]{};

    real faceAlignedStress22[seissol::multisim::NumSimulations]{};
    real faceAlignedStress33[seissol::multisim::NumSimulations]{};
    real faceAlignedStress12[seissol::multisim::NumSimulations]{};
    real faceAlignedStress13[seissol::multisim::NumSimulations]{};
    real faceAlignedStress23[seissol::multisim::NumSimulations]{};

    real updatedTraction1[seissol::multisim::NumSimulations]{};
    real updatedTraction2[seissol::multisim::NumSimulations]{};

    real slipRateStrike[seissol::multisim::NumSimulations]{};
    real slipRateDip[seissol::multisim::NumSimulations]{};

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
  void getNeighborDofs(real dofs[tensor::Q::size()], int meshId, int side);
  void computeLocalStresses(LocalInfo& local);
  virtual std::array<real, seissol::multisim::NumSimulations> computeLocalStrength(LocalInfo& local) = 0;
  virtual real computeFluidPressure(LocalInfo& local, unsigned int index) { return 0.0; }
  virtual real computeStateVariable(LocalInfo& local, unsigned int index) { return 0.0; }
  static void updateLocalTractions(LocalInfo& local, const std::array<real, seissol::multisim::NumSimulations>& strength);
  real computeRuptureVelocity(Eigen::Matrix<real, 2, 2>& jacobiT2d, const LocalInfo& local);
  virtual void computeSlipRate(LocalInfo& local,
                               const std::array<std::array<real, 6>, seissol::multisim::NumSimulations>& /*rotatedUpdatedStress*/,
                               const std::array<std::array<real, 6>, seissol::multisim::NumSimulations>& /*rotatedStress*/);
  static void computeSlipRate(LocalInfo& local,
                              const std::array<double, 3>& tangent1,
                              const std::array<double, 3>& tangent2,
                              const std::array<double, 3>& strike,
                              const std::array<double, 3>& dip);
  virtual void outputSpecifics(std::shared_ptr<ReceiverOutputData>& data,
                               const LocalInfo& local,
                               size_t outputSpecifics,
                               size_t receiverIdx) {}
  virtual void adjustRotatedUpdatedStress(std::array<std::array<real, 6>, seissol::multisim::NumSimulations>& rotatedUpdatedStress,
                                          const std::array<std::array<real, 6>, seissol::multisim::NumSimulations>& rotatedStress) {};
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
