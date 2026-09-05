// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/ParametersInitializer.h"
#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Parallel/Runtime/Stream.h"

#include <vector>

namespace seissol::dr::output {
/**
  Type-erased interface of the on-fault receiver output, plus everything which is shared by all
  friction laws.

  Only the entry points which the OutputManager calls through a `unique_ptr<ReceiverOutput>` are
  virtual. The per-point customisation points live in ReceiverOutputImpl below and are resolved
  statically, so that they can be turned into templates over the cell configuration later on.
 */
class ReceiverOutput {
  public:
  virtual ~ReceiverOutput() = default;

  void setLtsData(LTS::Storage& userWpStorage,
                  LTS::Backmap& userWpBackmap,
                  DynamicRupture::Storage& userDrStorage);

  void setMeshReader(seissol::geometry::MeshReader* userMeshReader) {
    meshReader_ = userMeshReader;
  }
  void setFaceToLtsMap(::seissol::initializer::StorageBackmap<1>* map) { faceToLtsMap_ = map; }

  virtual void
      calcFaultOutput(seissol::initializer::parameters::OutputType outputType,
                      seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
                      const std::shared_ptr<ReceiverOutputData>& outputData,
                      parallel::runtime::StreamRuntime& runtime,
                      double time = 0.0,
                      double dt = 1.0,
                      double indt = 0.0) = 0;

  [[nodiscard]] virtual std::vector<std::size_t> getOutputVariables() const;

  protected:
  LTS::Storage* wpStorage_{nullptr};
  LTS::Backmap* wpBackmap_{nullptr};
  DynamicRupture::Storage* drStorage_{nullptr};
  seissol::geometry::MeshReader* meshReader_{nullptr};
  ::seissol::initializer::StorageBackmap<1>* faceToLtsMap_{nullptr};
  real* deviceCopyMemory_{nullptr};

  kernels::Time timeKernel_;

  bool printRSFWarning_{false};

  struct LocalInfo {
    DynamicRupture::Layer* layer{};
    size_t ltsId{};
    int nearestGpIndex{};
    int nearestInternalGpIndex{};
    int gpIndex{};
    int internalGpIndexFused{};

    double time{};
    bool* printWarning{nullptr};

    std::size_t index{};
    std::size_t faceId{};
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

    real
        faceAlignedValuesPlus[tensor::QAtPoint::Shape[seissol::multisim::BasisFunctionDimension]]{};
    real faceAlignedValuesMinus
        [tensor::QAtPoint::Shape[seissol::multisim::BasisFunctionDimension]]{};

    model::IsotropicWaveSpeeds* waveSpeedsPlus{};
    model::IsotropicWaveSpeeds* waveSpeedsMinus{};

    ReceiverOutputData* state{};
  };

  /**
    Gets the cell data defined by the type StorageT.
    (we cannot just access the storage data structure in case we need to sparsely copy data for the
    onfault receiver output on GPUs)
   */
  template <typename StorageT>
  [[nodiscard]] const std::remove_extent_t<typename StorageT::Type>*
      getCellData(const LocalInfo& local) const {
    const auto devVar = local.state->deviceVariables.find(drStorage_->info<StorageT>().index);
    if (devVar != local.state->deviceVariables.end()) {
      return reinterpret_cast<const std::remove_extent_t<typename StorageT::Type>*>(
          devVar->second->get(local.faceId));
    } else {
      return local.layer->var<StorageT>()[local.ltsId];
    }
  }

  void getDofs(const real*(&derivatives), std::size_t meshId);
  void getNeighborDofs(const real*(&derivatives), std::size_t meshId, std::size_t side);
  void computeLocalStresses(LocalInfo& local);
  static void updateLocalTractions(LocalInfo& local, real strength);
  real computeRuptureVelocity(const Eigen::Matrix<real, 2, 2>& jacobiT2d, const LocalInfo& local);
  void computeSlipRate(LocalInfo& local,
                       const std::array<real, 6>& rotatedUpdatedStress,
                       const std::array<real, 6>& rotatedStress);
  static void computeSlipRate(LocalInfo& local,
                              const std::array<double, 3>& tangent1,
                              const std::array<double, 3>& tangent2,
                              const std::array<double, 3>& strike,
                              const std::array<double, 3>& dip);
};

/**
  Implements the output loop for one friction law, calling back into `Derived` for the parts which
  differ between them.

  The customisation points are ordinary member functions, not virtual ones: a derived class simply
  declares the ones it needs and thereby shadows the default below. Two consequences: the hooks
  have to be public in the derived class, since a base cannot reach a protected member of its own
  derived class; and there is no `override` to catch a misspelled hook, so the name has to match
  exactly.

  computeLocalStrength has no default on purpose -- a friction law which does not provide it fails
  to compile, the same way the pure virtual did.
 */
template <typename Derived>
class ReceiverOutputImpl : public ReceiverOutput {
  public:
  void calcFaultOutput(seissol::initializer::parameters::OutputType outputType,
                       seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
                       const std::shared_ptr<ReceiverOutputData>& outputData,
                       parallel::runtime::StreamRuntime& runtime,
                       double time = 0.0,
                       double dt = 1.0,
                       double indt = 0.0) override;

  protected:
  real computeFluidPressure(LocalInfo& /*local*/) { return 0.0; }
  real computeStateVariable(LocalInfo& /*local*/) { return 0.0; }
  void outputSpecifics(const std::shared_ptr<ReceiverOutputData>& /*data*/,
                       const LocalInfo& /*local*/,
                       size_t /*cacheLevel*/,
                       size_t /*receiverIdx*/) {}
  void adjustRotatedUpdatedStress(std::array<real, 6>& /*rotatedUpdatedStress*/,
                                  const std::array<real, 6>& /*rotatedStress*/) {}
  void handleNonConvergence(LocalInfo& /*local*/) {}

  private:
  Derived& derived() { return static_cast<Derived&>(*this); }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
