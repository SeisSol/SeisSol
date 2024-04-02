#ifndef SEISSOL_FRICTIONSOLVER_H
#define SEISSOL_FRICTIONSOLVER_H

#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/DynamicRupture.h"

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
                        seissol::initializer::DynamicRupture const* const dynRup,
                        real fullUpdateTime,
                        const double timeWeights[ConvergenceOrder]) = 0;

  /**
   * compute the DeltaT from the current timePoints call this function before evaluate
   * to set the correct DeltaT
   */
  void computeDeltaT(const double timePoints[ConvergenceOrder]);

  /**
   * copies all common parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          seissol::initializer::DynamicRupture const* const dynRup,
                          real fullUpdateTime);

  protected:
  /**
   * Adjust initial stress by adding nucleation stress * nucleation function
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf.
   */
  real deltaT[ConvergenceOrder] = {};
  real sumDt;

  seissol::initializer::parameters::DRParameters* drParameters;
  ImpedancesAndEta* impAndEta;
  ImpedanceMatrices* impedanceMatrices;
  real mFullUpdateTime;
  // CS = coordinate system
  real (*initialStressInFaultCS)[misc::numPaddedPoints][6];
  real (*nucleationStressInFaultCS)[misc::numPaddedPoints][6];
  real (*cohesion)[misc::numPaddedPoints];
  real (*mu)[misc::numPaddedPoints];
  real (*accumulatedSlipMagnitude)[misc::numPaddedPoints];
  real (*slip1)[misc::numPaddedPoints];
  real (*slip2)[misc::numPaddedPoints];
  real (*slipRateMagnitude)[misc::numPaddedPoints];
  real (*slipRate1)[misc::numPaddedPoints];
  real (*slipRate2)[misc::numPaddedPoints];
  real (*ruptureTime)[misc::numPaddedPoints];
  bool (*ruptureTimePending)[misc::numPaddedPoints];
  real (*peakSlipRate)[misc::numPaddedPoints];
  real (*traction1)[misc::numPaddedPoints];
  real (*traction2)[misc::numPaddedPoints];
  real (*imposedStatePlus)[tensor::QInterpolated::size()];
  real (*imposedStateMinus)[tensor::QInterpolated::size()];
  real spaceWeights[misc::numPaddedPoints];
  DREnergyOutput* energyData{};
  DRGodunovData* godunovData{};
  real (*initialPressure)[misc::numPaddedPoints];
  real (*nucleationPressure)[misc::numPaddedPoints];

  // be careful only for some FLs initialized:
  real (*dynStressTime)[misc::numPaddedPoints];
  bool (*dynStressTimePending)[misc::numPaddedPoints];

  real (*qInterpolatedPlus)[ConvergenceOrder][tensor::QInterpolated::size()];
  real (*qInterpolatedMinus)[ConvergenceOrder][tensor::QInterpolated::size()];
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_FRICTIONSOLVER_H
