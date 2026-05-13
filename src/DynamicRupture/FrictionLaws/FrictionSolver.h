// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Parallel/Runtime/Stream.h"

#include <vector>

namespace seissol::dr::friction_law {
/**
 * Abstract Base for friction solver class with the public interface
 * Only needed to be able to store a shared_ptr<FrictionSolver> in MemoryManager and TimeCluster.
 * BaseFrictionLaw has a template argument for CRTP, hence, we can't store a pointer to any
 * BaseFrictionLaw.
 *
 * Note: this class, FrictionSolver, must be trivially copyable. It is (or was?) important for GPU
 * offloading.
 */
class FrictionSolver {
  public:
  explicit FrictionSolver(const FrictionLawParameters& userDRParameters)
      : drParameters_(userDRParameters) {}
  virtual ~FrictionSolver() = default;

  struct FrictionTime {
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
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  virtual void allocateAuxiliaryMemory(GlobalData* globalData) {
    spaceWeights_ = globalData->spaceWeights;
  }

  virtual seissol::initializer::AllocationPlace allocationPlace();

  virtual std::unique_ptr<FrictionSolver> clone() = 0;

  protected:
  /**
   * Adjust initial stress by adding nucleation stress * nucleation function
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf.
   */
  real deltaT_[misc::TimeSteps] = {};

  FrictionLawParameters drParameters_;
  ImpedancesAndEta* __restrict impAndEta_{};
  ImpedanceMatrices* __restrict impedanceMatrices_{};
  real fullUpdateTime_{};
  // CS = coordinate system
  real (*__restrict initialStressInFaultCS_)[6][misc::NumPaddedPoints]{};
  real (*__restrict nucleationStressInFaultCS_)[6][misc::NumPaddedPoints]{};
  real (*__restrict cohesion_)[misc::NumPaddedPoints]{};
  real (*__restrict mu_)[misc::NumPaddedPoints]{};
  real (*__restrict accumulatedSlipMagnitude_)[misc::NumPaddedPoints]{};
  real (*__restrict slip1_)[misc::NumPaddedPoints]{};
  real (*__restrict slip2_)[misc::NumPaddedPoints]{};
  real (*__restrict slipRateMagnitude_)[misc::NumPaddedPoints]{};
  real (*__restrict slipRate1_)[misc::NumPaddedPoints]{};
  real (*__restrict slipRate2_)[misc::NumPaddedPoints]{};
  real (*__restrict ruptureTime_)[misc::NumPaddedPoints]{};
  bool (*__restrict ruptureTimePending_)[misc::NumPaddedPoints]{};
  real (*__restrict peakSlipRate_)[misc::NumPaddedPoints]{};
  real (*__restrict traction1_)[misc::NumPaddedPoints]{};
  real (*__restrict traction2_)[misc::NumPaddedPoints]{};
  real (*__restrict imposedStatePlus_)[tensor::QInterpolated::size()]{};
  real (*__restrict imposedStateMinus_)[tensor::QInterpolated::size()]{};
  real* __restrict spaceWeights_{};
  DREnergyOutput* __restrict energyData_{};
  DRGodunovData* __restrict godunovData_{};
  real (*__restrict initialPressure_)[misc::NumPaddedPoints]{};
  real (*__restrict nucleationPressure_)[misc::NumPaddedPoints]{};

  // be careful only for some FLs initialized:
  real (*__restrict dynStressTime_)[misc::NumPaddedPoints]{};
  bool (*__restrict dynStressTimePending_)[misc::NumPaddedPoints]{};

  real (*__restrict qInterpolatedPlus_)[misc::TimeSteps][tensor::QInterpolated::size()]{};
  real (*__restrict qInterpolatedMinus_)[misc::TimeSteps][tensor::QInterpolated::size()]{};
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
