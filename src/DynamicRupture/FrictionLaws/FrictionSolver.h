// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_

#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Parallel/Runtime/Stream.h"
#include <vector>

namespace seissol::dr::friction_law {
/**
 * Abstract Base for friction solver class with the public interface
 * Only needed to be able to store a shared_ptr<FrictionSolver> in MemoryManager and TimeCluster.
 * BaseFrictionLaw has a template argument for CRTP, hence, we can't store a pointer to any
 * BaseFrictionLaw.
 */
class FrictionSolver {
  public:
  // Note: FrictionSolver must be trivially copyable. It is important for GPU offloading
  FrictionSolver() = default;
  virtual ~FrictionSolver() = default;

  struct FrictionTime {
    double sumDt;
    std::vector<double> deltaT;
  };

  virtual void setupLayer(DynamicRupture::Layer& layerData,
                          seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void evaluate(double fullUpdateTime,
                        const FrictionTime& frictionTime,
                        const double* timeWeights,
                        seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  /**
   * compute the DeltaT from the current timePoints call this function before evaluate
   * to set the correct DeltaT
   */
  static FrictionTime computeDeltaT(const std::vector<double>& timePoints);

  /**
   * copies all common parameters from the DynamicRupture LTS to the local attributes
   */
  virtual void copyStorageToLocal(DynamicRupture::Layer& layerData) = 0;

  virtual void allocateAuxiliaryMemory(GlobalData* globalData) {}

  virtual seissol::initializer::AllocationPlace allocationPlace() = 0;

  virtual std::unique_ptr<FrictionSolver> clone() = 0;
};

class FrictionSolverImpl : public FrictionSolver {
  public:
  explicit FrictionSolverImpl(seissol::initializer::parameters::DRParameters* userDRParameters)
      : drParameters(userDRParameters) {
    std::copy(
        &init::quadweights<
            Cfg>::Values[init::quadweights<Cfg>::Start[seissol::multisim::BasisFunctionDimension]],
        &init::quadweights<
            Cfg>::Values[init::quadweights<Cfg>::Stop[seissol::multisim::BasisFunctionDimension]],
        &spaceWeights[0]);
  }

  void copyStorageToLocal(DynamicRupture::Layer& layerData) override;

  seissol::initializer::AllocationPlace allocationPlace() override;

  protected:
  /**
   * Adjust initial stress by adding nucleation stress * nucleation function
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf.
   */
  real deltaT[Cfg::ConvergenceOrder] = {};
  real sumDt{};

  seissol::initializer::parameters::DRParameters* __restrict drParameters;
  ImpedancesAndEta<Real<Cfg>>* __restrict impAndEta{};
  ImpedanceMatrices<Cfg>* __restrict impedanceMatrices{};
  real mFullUpdateTime{};
  // CS = coordinate system
  real (*__restrict initialStressInFaultCS)[6][misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict nucleationStressInFaultCS)[6][misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict cohesion)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict mu)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict accumulatedSlipMagnitude)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slip1)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slip2)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRateMagnitude)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRate1)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRate2)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict ruptureTime)[misc::NumPaddedPoints<Cfg>]{};
  bool (*__restrict ruptureTimePending)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict peakSlipRate)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict traction1)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict traction2)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict imposedStatePlus)[tensor::QInterpolated<Cfg>::size()]{};
  real (*__restrict imposedStateMinus)[tensor::QInterpolated<Cfg>::size()]{};
  real spaceWeights[misc::NumPaddedPoints<Cfg>]{};
  DREnergyOutput<Cfg>* __restrict energyData{};
  DRGodunovData<Cfg>* __restrict godunovData{};
  real (*__restrict initialPressure)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict nucleationPressure)[misc::NumPaddedPoints<Cfg>]{};

  // be careful only for some FLs initialized:
  real (*__restrict dynStressTime)[misc::NumPaddedPoints<Cfg>]{};
  bool (*__restrict dynStressTimePending)[misc::NumPaddedPoints<Cfg>]{};

  real (*__restrict qInterpolatedPlus)[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()]{};
  real (*__restrict qInterpolatedMinus)[Cfg::ConvergenceOrder]
                                       [tensor::QInterpolated<Cfg>::size()]{};
};

using FrictionSolverFactory = std::function<std::unique_ptr<FrictionSolver>(ConfigVariant)>;

template <typename T, typename... Args>
FrictionSolverFactory makeFrictionSolverFactory(Args... args) {
  return [=](ConfigVariant variant) {
    return std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          return std::make_unique<T>(args...);
        },
        variant);
  };
}

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
