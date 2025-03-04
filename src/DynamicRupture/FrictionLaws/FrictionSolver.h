// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_

#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/DynamicRupture.h"
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

  seissol::initializer::parameters::DRParameters* drParameters;
  ImpedancesAndEta* impAndEta{};
  ImpedanceMatrices* impedanceMatrices{};
  real mFullUpdateTime{};
  // CS = coordinate system
  real (*initialStressInFaultCS)[misc::NumPaddedPoints][6]{};
  real (*nucleationStressInFaultCS)[misc::NumPaddedPoints][6]{};
  real (*cohesion)[misc::NumPaddedPoints]{};
  real (*mu)[misc::NumPaddedPoints]{};
  real (*accumulatedSlipMagnitude)[misc::NumPaddedPoints]{};
  real (*slip1)[misc::NumPaddedPoints]{};
  real (*slip2)[misc::NumPaddedPoints]{};
  real (*slipRateMagnitude)[misc::NumPaddedPoints]{};
  real (*slipRate1)[misc::NumPaddedPoints]{};
  real (*slipRate2)[misc::NumPaddedPoints]{};
  real (*ruptureTime)[misc::NumPaddedPoints]{};
  bool (*ruptureTimePending)[misc::NumPaddedPoints]{};
  real (*peakSlipRate)[misc::NumPaddedPoints]{};
  real (*traction1)[misc::NumPaddedPoints]{};
  real (*traction2)[misc::NumPaddedPoints]{};
  real (*imposedStatePlus)[tensor::QInterpolated::size()]{};
  real (*imposedStateMinus)[tensor::QInterpolated::size()]{};
  real spaceWeights[misc::NumPaddedPoints]{};
  DREnergyOutput* energyData{};
  DRGodunovData* godunovData{};
  real (*initialPressure)[misc::NumPaddedPoints]{};
  real (*nucleationPressure)[misc::NumPaddedPoints]{};

  // be careful only for some FLs initialized:
  real (*dynStressTime)[misc::NumPaddedPoints]{};
  bool (*dynStressTimePending)[misc::NumPaddedPoints]{};

  real (*qInterpolatedPlus)[ConvergenceOrder][tensor::QInterpolated::size()]{};
  real (*qInterpolatedMinus)[ConvergenceOrder][tensor::QInterpolated::size()]{};
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
