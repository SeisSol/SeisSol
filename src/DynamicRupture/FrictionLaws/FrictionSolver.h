#ifndef SEISSOL_FRICTIONSOLVER_H
#define SEISSOL_FRICTIONSOLVER_H

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/DynamicRupture.h"
#include "Common/constants.hpp"

namespace seissol::dr::friction_law {
/**
 * Abstract Base for friction solver class with the public interface
 * Only needed to be able to store a shared_ptr<FrictionSolver> in MemoryManager and TimeCluster.
 * BaseFrictionLaw has a template argument for CRTP, hence, we can't store a pointer to any
 * BaseFrictionLaw.
 */
template <typename Config>
class FrictionSolver {
  public:
  using RealT = typename Config::RealT;
  // Note: FrictionSolver must be trivially copyable. It is important for GPU offloading
  explicit FrictionSolver(dr::DRParameters* userDrParameters) : drParameters(userDrParameters) {
    std::copy(
        &Yateto<Config>::Init::quadweights::Values[Yateto<Config>::Init::quadweights::Start[0]],
        &Yateto<Config>::Init::quadweights::Values[Yateto<Config>::Init::quadweights::Stop[0]],
        &spaceWeights[0]);
  }
  virtual ~FrictionSolver() = default;

  virtual void evaluate(seissol::initializers::Layer& layerData,
                        seissol::initializers::DynamicRupture<Config> const* const dynRup,
                        RealT fullUpdateTime,
                        const double timeWeights[Config::ConvergenceOrder]) = 0;

  /**
   * compute the DeltaT from the current timePoints call this function before evaluate
   * to set the correct DeltaT
   */
  void computeDeltaT(const double timePoints[Config::ConvergenceOrder]);

  /**
   * copies all common parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime);

  protected:
  /**
   * Adjust initial stress by adding nucleation stress * nucleation function
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf.
   */
  RealT deltaT[Config::ConvergenceOrder] = {};

  dr::DRParameters* drParameters;
  ImpedancesAndEta<Config>* impAndEta;
  RealT mFullUpdateTime;
  // CS = coordinate system
  RealT (*initialStressInFaultCS)[misc::numPaddedPoints<Config>][6];
  RealT (*nucleationStressInFaultCS)[misc::numPaddedPoints<Config>][6];
  RealT (*cohesion)[misc::numPaddedPoints<Config>];
  RealT (*mu)[misc::numPaddedPoints<Config>];
  RealT (*accumulatedSlipMagnitude)[misc::numPaddedPoints<Config>];
  RealT (*slip1)[misc::numPaddedPoints<Config>];
  RealT (*slip2)[misc::numPaddedPoints<Config>];
  RealT (*slipRateMagnitude)[misc::numPaddedPoints<Config>];
  RealT (*slipRate1)[misc::numPaddedPoints<Config>];
  RealT (*slipRate2)[misc::numPaddedPoints<Config>];
  RealT (*ruptureTime)[misc::numPaddedPoints<Config>];
  bool (*ruptureTimePending)[misc::numPaddedPoints<Config>];
  RealT (*peakSlipRate)[misc::numPaddedPoints<Config>];
  RealT (*traction1)[misc::numPaddedPoints<Config>];
  RealT (*traction2)[misc::numPaddedPoints<Config>];
  RealT (*imposedStatePlus)[Yateto<Config>::Tensor::QInterpolated::size()];
  RealT (*imposedStateMinus)[Yateto<Config>::Tensor::QInterpolated::size()];
  RealT spaceWeights[misc::numPaddedPoints<Config>];
  DREnergyOutput<Config>* energyData{};
  DRGodunovData<Config>* godunovData{};

  // be careful only for some FLs initialized:
  RealT (*dynStressTime)[misc::numPaddedPoints<Config>];
  bool (*dynStressTimePending)[misc::numPaddedPoints<Config>];

  RealT (
      *qInterpolatedPlus)[Config::ConvergenceOrder][Yateto<Config>::Tensor::QInterpolated::size()];
  RealT (
      *qInterpolatedMinus)[Config::ConvergenceOrder][Yateto<Config>::Tensor::QInterpolated::size()];
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_FRICTIONSOLVER_H
