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

#include <DynamicRupture/Misc.h>
#include <Memory/Tree/Backmap.h>
#include <Parallel/Runtime/Stream.h>
#include <memory>
#include <vector>

namespace seissol::dr::output {
class ReceiverOutput {
  public:
  virtual ~ReceiverOutput() = default;

  void setLtsData(LTS::Storage& userWpStorage,
                  LTS::Backmap& userWpBackmap,
                  DynamicRupture::Storage& userDrStorage);

  void setMeshReader(seissol::geometry::MeshReader* userMeshReader) { meshReader = userMeshReader; }
  void setFaceToLtsMap(FaceToLtsMapType* map) { faceToLtsMap = map; }
  virtual void
      calcFaultOutput(seissol::initializer::parameters::OutputType outputType,
                      seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
                      const std::shared_ptr<ReceiverOutputData>& outputData,
                      parallel::runtime::StreamRuntime& runtime,
                      double time = 0.0) = 0;

  [[nodiscard]] virtual std::vector<std::size_t> getOutputVariables() const = 0;

  protected:
  LTS::Storage* wpStorage{nullptr};
  LTS::Backmap* wpBackmap{nullptr};
  DynamicRupture::Storage* drStorage{nullptr};
  seissol::geometry::MeshReader* meshReader{nullptr};
  FaceToLtsMapType* faceToLtsMap{nullptr};
};

template <typename Derived>
class ReceiverOutputImpl : public ReceiverOutput {
  public:
  void calcFaultOutput(seissol::initializer::parameters::OutputType outputType,
                       seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
                       const std::shared_ptr<ReceiverOutputData>& outputData,
                       parallel::runtime::StreamRuntime& runtime,
                       double time = 0.0) override;

  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const override;

  protected:
  template <typename Cfg>
  struct LocalInfo {
    using real = Real<Cfg>;

    DynamicRupture::Layer* layer{};
    size_t ltsId{};
    int nearestGpIndex{};
    int nearestInternalGpIndex{};
    int gpIndex{};
    int internalGpIndexFused{};

    std::size_t index{};
    std::size_t fusedIndex{};

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

    real faceAlignedValuesPlus
        [tensor::QAtPoint<Cfg>::Shape[seissol::multisim::BasisFunctionDimension]]{};
    real faceAlignedValuesMinus
        [tensor::QAtPoint<Cfg>::Shape[seissol::multisim::BasisFunctionDimension]]{};

    model::IsotropicWaveSpeeds* waveSpeedsPlus{};
    model::IsotropicWaveSpeeds* waveSpeedsMinus{};

    ReceiverOutputData* state{};
  };

  template <typename StorageT, typename Cfg>
  std::remove_extent_t<typename StorageT::template VariantType<Cfg>>*
      getCellData(Cfg cfg, const LocalInfo<Cfg>& local) {
    auto devVar = local.state->deviceVariables.find(drStorage->info<StorageT>().index);
    if (devVar != local.state->deviceVariables.end()) {
      return reinterpret_cast<std::remove_extent_t<typename StorageT::template VariantType<Cfg>>*>(
          devVar->second->get(local.state->deviceIndices[local.index]));
    } else {
      return local.layer->template var<StorageT>(cfg)[local.ltsId];
    }
  }

  template <typename Cfg>
  void getDofs(Real<Cfg> dofs[tensor::Q<Cfg>::size()], int meshId);

  template <typename Cfg>
  void getNeighborDofs(Real<Cfg> dofs[tensor::Q<Cfg>::size()], int meshId, int side);

  template <typename Cfg>
  void computeLocalStresses(LocalInfo<Cfg>& local);

  template <typename Cfg>
  static void updateLocalTractions(LocalInfo<Cfg>& local, Real<Cfg> strength);

  template <typename Cfg>
  Real<Cfg> computeRuptureVelocity(const Eigen::Matrix<Real<Cfg>, 2, 2>& jacobiT2d,
                                   const LocalInfo<Cfg>& local);

  template <typename Cfg>
  void computeSlipRate(LocalInfo<Cfg>& local,
                       const std::array<Real<Cfg>, 6>& /*rotatedUpdatedStress*/,
                       const std::array<Real<Cfg>, 6>& /*rotatedStress*/);

  template <typename Cfg>
  static void computeSlipRate(LocalInfo<Cfg>& local,
                              const std::array<double, 3>& tangent1,
                              const std::array<double, 3>& tangent2,
                              const std::array<double, 3>& strike,
                              const std::array<double, 3>& dip);

  template <typename Cfg>
  Real<Cfg> computeLocalStrength(LocalInfo<Cfg>& local) {
    return 0;
  }

  template <typename Cfg>
  Real<Cfg> computeFluidPressure(LocalInfo<Cfg>& local) {
    return 0;
  }

  template <typename Cfg>
  Real<Cfg> computeStateVariable(LocalInfo<Cfg>& local) {
    return 0;
  }

  template <typename Cfg>
  void outputSpecifics(const std::shared_ptr<ReceiverOutputData>& data,
                       const LocalInfo<Cfg>& local,
                       size_t outputSpecifics,
                       size_t receiverIdx) {}

  template <typename Cfg>
  void adjustRotatedUpdatedStress(std::array<Real<Cfg>, 6>& rotatedUpdatedStress,
                                  const std::array<Real<Cfg>, 6>& rotatedStress) {}
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
