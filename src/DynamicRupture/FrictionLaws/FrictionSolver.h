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
  explicit FrictionSolver(seissol::initializer::parameters::DRParameters* userDRParameters)
      : drParameters(userDRParameters) {
    std::copy(&init::quadweights::Values[init::quadweights::Start[0]],
              &init::quadweights::Values[init::quadweights::Stop[0]],
              &spaceWeights[0]);
  }
  virtual ~FrictionSolver() = default;

  virtual void evaluate(seissol::initializer::Layer& layerData,
                        const seissol::initializer::DynamicRupture* dynRup,
                        real fullUpdateTime,
                        const double timeWeights[ConvergenceOrder],
                        seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  /**
   * compute the DeltaT from the current timePoints call this function before evaluate
   * to set the correct DeltaT
   */
  void computeDeltaT(const double timePoints[ConvergenceOrder]);

  /**
   * copies all common parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* dynRup,
                          real fullUpdateTime);

  virtual seissol::initializer::AllocationPlace allocationPlace();

  protected:
  /**
   * Adjust initial stress by adding nucleation stress * nucleation function
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf.
   */
  real deltaT[ConvergenceOrder] = {};
  real sumDt{};

  seissol::initializer::parameters::DRParameters* __restrict drParameters;
  ImpedancesAndEta* __restrict impAndEta{};
  ImpedanceMatrices* __restrict impedanceMatrices{};
  real mFullUpdateTime{};
  // CS = coordinate system
  real (*__restrict initialStressInFaultCS)[6][misc::NumPaddedPoints]{};
  real (*__restrict nucleationStressInFaultCS[seissol::initializer::parameters::MaxNucleactions])
      [6][misc::NumPaddedPoints]{};
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
  real spaceWeights[misc::NumPaddedPoints]{};
  DREnergyOutput* __restrict energyData{};
  DRGodunovData* __restrict godunovData{};
  real (*__restrict initialPressure)[misc::NumPaddedPoints]{};
  real (*__restrict nucleationPressure[initializer::parameters::MaxNucleactions])
      [misc::NumPaddedPoints]{};

  // be careful only for some FLs initialized:
  real (*__restrict dynStressTime)[misc::NumPaddedPoints]{};
  bool (*__restrict dynStressTimePending)[misc::NumPaddedPoints]{};

  real (*__restrict qInterpolatedPlus)[ConvergenceOrder][tensor::QInterpolated::size()]{};
  real (*__restrict qInterpolatedMinus)[ConvergenceOrder][tensor::QInterpolated::size()]{};
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
