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
      : drParameters(userDRParameters) {}
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
    spaceWeights = globalData->spaceWeights;
  }

  virtual seissol::initializer::AllocationPlace allocationPlace();

  virtual std::unique_ptr<FrictionSolver> clone() = 0;

  protected:
  /**
   * Adjust initial stress by adding nucleation stress * nucleation function
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf.
   */
  real deltaT[misc::TimeSteps] = {};

  FrictionLawParameters drParameters;
  ImpedancesAndEta* __restrict impAndEta{};
  ImpedanceMatrices* __restrict impedanceMatrices{};
  real mFullUpdateTime{};
  // CS = coordinate system
  real (*__restrict initialStressInFaultCS)[6][misc::NumPaddedPoints]{};
  real (*__restrict nucleationStressInFaultCS)[6][misc::NumPaddedPoints]{};
  real (*__restrict cohesion)[misc::NumPaddedPoints]{};
  real (*__restrict mu)[misc::NumPaddedPoints]{};
  real (*__restrict accumulatedSlipMagnitude)[misc::NumPaddedPoints]{};
  real (*__restrict slip1)[misc::NumPaddedPoints]{};
  real (*__restrict slip2)[misc::NumPaddedPoints]{};
  real (*__restrict slipRateMagnitude)[misc::NumPaddedPoints]{};
  real (*__restrict slipRate1)[misc::NumPaddedPoints]{};
  real (*__restrict slipRate2)[misc::NumPaddedPoints]{};
  real (*__restrict ruptureTime)[misc::NumPaddedPoints]{};
  bool (*__restrict ruptureTimePending)[misc::NumPaddedPoints]{};
  real (*__restrict peakSlipRate)[misc::NumPaddedPoints]{};
  real (*__restrict traction1)[misc::NumPaddedPoints]{};
  real (*__restrict traction2)[misc::NumPaddedPoints]{};
  real (*__restrict imposedStatePlus)[tensor::QInterpolated::size()]{};
  real (*__restrict imposedStateMinus)[tensor::QInterpolated::size()]{};
  real* __restrict spaceWeights{};
  DREnergyOutput* __restrict energyData{};
  DRGodunovData* __restrict godunovData{};
  real (*__restrict initialPressure)[misc::NumPaddedPoints]{};
  real (*__restrict nucleationPressure)[misc::NumPaddedPoints]{};

  // be careful only for some FLs initialized:
  real (*__restrict dynStressTime)[misc::NumPaddedPoints]{};
  bool (*__restrict dynStressTimePending)[misc::NumPaddedPoints]{};

  real (*__restrict qInterpolatedPlus)[misc::TimeSteps][tensor::QInterpolated::size()]{};
  real (*__restrict qInterpolatedMinus)[misc::TimeSteps][tensor::QInterpolated::size()]{};
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
