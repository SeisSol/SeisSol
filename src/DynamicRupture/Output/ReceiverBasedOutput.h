// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_

#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/ParametersInitializer.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Parallel/Runtime/Stream.h"

#include <vector>

namespace seissol::dr::output {
/**
  The output of the on-fault receivers, whatever the configurations of the fault faces.
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
  void setDrParameters(const seissol::initializer::parameters::DRParameters* userDrParameters) {
    drParameters_ = userDrParameters;
  }
  /**
   * @param stateTime the time the stored friction state belongs to, which is the end of the dynamic
   *                  rupture time step that computed it; the stress sources are evaluated there,
   *                  where the friction law evaluated them last, so that the tractions rebuilt
   *                  here are consistent with the state they are rebuilt from
   * @param time the time the output is recorded under
   */
  virtual void
      calcFaultOutput(seissol::initializer::parameters::OutputType outputType,
                      seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
                      const std::shared_ptr<ReceiverOutputData>& outputData,
                      parallel::runtime::StreamRuntime& runtime,
                      double stateTime,
                      double time = 0.0,
                      double dt = 1.0,
                      double indt = 0.0) = 0;

  [[nodiscard]] virtual std::vector<std::size_t> getOutputVariables() const;

  protected:
  LTS::Storage* wpStorage_{nullptr};
  LTS::Backmap* wpBackmap_{nullptr};
  DynamicRupture::Storage* drStorage_{nullptr};
  seissol::geometry::MeshReader* meshReader_{nullptr};
  const seissol::initializer::parameters::DRParameters* drParameters_{nullptr};
  ::seissol::initializer::StorageBackmap<1>* faceToLtsMap_{nullptr};

  bool printRSFWarning_{false};
};

/**
  The output of the on-fault receivers of a friction law `Derived` (CRTP). `Derived` provides
  `computeLocalStrength` and may provide its own versions of the other hooks below; each hook is a
  template of the configuration of the fault face it is evaluated for.
 */
template <typename Derived>
class ReceiverOutputImpl : public ReceiverOutput {
  public:
  void calcFaultOutput(seissol::initializer::parameters::OutputType outputType,
                       seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
                       const std::shared_ptr<ReceiverOutputData>& outputData,
                       parallel::runtime::StreamRuntime& runtime,
                       double stateTime,
                       double time = 0.0,
                       double dt = 1.0,
                       double indt = 0.0) override;

  template <typename Cfg>
  struct LocalInfo {
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

    DynamicRupture::Layer* layer{};
    size_t ltsId{};
    int nearestGpIndex{};
    int nearestInternalGpIndex{};
    int gpIndex{};
    int internalGpIndexFused{};

    double time{};
    /// width of the last sub time step of the friction solve, which is the one the stored friction
    /// state belongs to
    double deltaT{};
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
    /// slip rate in the fault-local frame, filled where the friction reconstruction resolves the
    /// slip direction itself instead of inheriting it from the trial traction
    real slipRateTangent1{};
    real slipRateTangent2{};

    real faceAlignedValuesPlus[tensor::QAtPoint<Cfg>::Shape[seissol::multisim::BasisDim<Cfg>]]{};
    real faceAlignedValuesMinus[tensor::QAtPoint<Cfg>::Shape[seissol::multisim::BasisDim<Cfg>]]{};

    model::IsotropicWaveSpeeds* waveSpeedsPlus{};
    model::IsotropicWaveSpeeds* waveSpeedsMinus{};

    ReceiverOutputData* state{};
  };

  /**
    Gets the cell data defined by the type StorageT.
    (we cannot just access the storage data structure in case we need to sparsely copy data for the
    onfault receiver output on GPUs)
   */
  template <typename StorageT, typename Cfg>
  [[nodiscard]] const std::remove_extent_t<seissol::initializer::StorageType<StorageT, Cfg>>*
      getCellData(const LocalInfo<Cfg>& local) const {
    using ValueT = std::remove_extent_t<seissol::initializer::StorageType<StorageT, Cfg>>;
    const auto devVar = local.state->deviceVariables.find(drStorage_->info<StorageT>().index);
    if (devVar != local.state->deviceVariables.end()) {
      return reinterpret_cast<const ValueT*>(devVar->second->get(local.faceId));
    } else {
      return local.layer->template var<StorageT>(Cfg())[local.ltsId];
    }
  }

  /**
    d(strength) / d(-sigma_eff), the counterpart of the friction laws' strengthSlope. Only read for
    materials whose impedance couples shear slip to the fault-normal traction; zero means that the
    strength does not follow the normal stress.
   */
  template <typename Cfg>
  Real<Cfg> computeLocalStrengthSlope(LocalInfo<Cfg>& /*local*/) {
    return 0.0;
  }
  template <typename Cfg>
  Real<Cfg> computeFluidPressure(LocalInfo<Cfg>& /*local*/) {
    return 0.0;
  }
  template <typename Cfg>
  Real<Cfg> computeStateVariable(LocalInfo<Cfg>& /*local*/) {
    return 0.0;
  }
  template <typename Cfg>
  void computeSlipRate(LocalInfo<Cfg>& local,
                       const std::array<Real<Cfg>, 6>& /*rotatedUpdatedStress*/,
                       const std::array<Real<Cfg>, 6>& /*rotatedStress*/,
                       const std::array<double, 3>& /*tangent1*/,
                       const std::array<double, 3>& /*tangent2*/,
                       const std::array<double, 3>& /*strike*/,
                       const std::array<double, 3>& /*dip*/);
  template <typename Cfg>
  void outputSpecifics(const std::shared_ptr<ReceiverOutputData>& data,
                       const LocalInfo<Cfg>& local,
                       size_t outputSpecifics,
                       size_t receiverIdx) {}
  template <typename Cfg>
  void adjustRotatedUpdatedStress(std::array<Real<Cfg>, 6>& rotatedUpdatedStress,
                                  const std::array<Real<Cfg>, 6>& rotatedStress) {}
  template <typename Cfg>
  void handleNonConvergence(LocalInfo<Cfg>& local) {}

  protected:
  template <typename Cfg>
  void getDofs(const Real<Cfg>*(&derivatives), std::size_t meshId);
  template <typename Cfg>
  void getNeighborDofs(const Real<Cfg>*(&derivatives), std::size_t meshId, std::size_t side);
  template <typename Cfg>
  void computeLocalStresses(LocalInfo<Cfg>& local);
  template <typename Cfg>
  static void
      updateLocalTractions(LocalInfo<Cfg>& local, Real<Cfg> strength, Real<Cfg> strengthSlope);
  template <typename Cfg>
  Real<Cfg> computeRuptureVelocity(const Eigen::Matrix<Real<Cfg>, 2, 2>& jacobiT2d,
                                   const LocalInfo<Cfg>& local);
  template <typename Cfg>
  static void computeSlipRate(LocalInfo<Cfg>& local,
                              const std::array<double, 3>& tangent1,
                              const std::array<double, 3>& tangent2,
                              const std::array<double, 3>& strike,
                              const std::array<double, 3>& dip);
  /// Writes a fault plane vector given in the (tangent1, tangent2) frame to the slip rate along
  /// strike and dip.
  template <typename Cfg>
  static void projectOntoStrikeAndDip(LocalInfo<Cfg>& local,
                                      Real<Cfg> alongTangent1,
                                      Real<Cfg> alongTangent2,
                                      const std::array<double, 3>& tangent1,
                                      const std::array<double, 3>& tangent2,
                                      const std::array<double, 3>& strike,
                                      const std::array<double, 3>& dip);

  private:
  template <typename Cfg>
  void calcFaultOutputOfConfig(
      seissol::initializer::parameters::OutputType outputType,
      seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
      const std::shared_ptr<ReceiverOutputData>& outputData,
      parallel::runtime::StreamRuntime& runtime,
      double stateTime,
      double dt,
      double indt);

  Derived& derived() { return static_cast<Derived&>(*this); }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RECEIVERBASEDOUTPUT_H_
