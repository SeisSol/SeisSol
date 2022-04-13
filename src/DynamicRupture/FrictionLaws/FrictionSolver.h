#ifndef SEISSOL_FRICTIONSOLVER_H
#define SEISSOL_FRICTIONSOLVER_H

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/DynamicRupture.h"
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
  FrictionSolver(dr::DRParameters& userDrParameters) : drParameters(userDrParameters){};
  virtual ~FrictionSolver(){};

  virtual void evaluate(seissol::initializers::Layer& layerData,
                        seissol::initializers::DynamicRupture* dynRup,
                        real fullUpdateTime,
                        double timeWeights[CONVERGENCE_ORDER]) = 0;

  /**
   * compute the DeltaT from the current timePoints call this function before evaluate
   * to set the correct DeltaT
   */
  void computeDeltaT(double timePoints[CONVERGENCE_ORDER]);

  /**
   * copies all common parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  /**
   * Calculate traction and normal stress at the interface.
   * Using equations (A2) from Pelties et al. 2014
   * Definiton of eta and impedance Z are found in dissertation of Carsten Uphoff
   * @param ltsFace: current fault face to be evaluated
   * @returns
   * NormalStress XYStress, XZStress at the 2d face quadrature nodes evaluated at the time
   * quadrature points
   */
  FaultStresses precomputeStressFromQInterpolated(unsigned int ltsFace);

  /**
   * Integrate over all Time points with the time weights and calculate the traction for each side
   * according to Carsten Uphoff Thesis: EQ.: 4.60
   * IN: NormalStressGP, XYTractionResultGP, * XZTractionResultGP
   * OUT: imposedStatePlus, imposedStateMinus
   */
  void postcomputeImposedStateFromNewStress(const FaultStresses& faultStresses,
                                            const TractionResults& tractionResults,
                                            double timeWeights[CONVERGENCE_ORDER],
                                            unsigned int ltsFace);

  /**
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
   */
  real calcSmoothStepIncrement(real currentTime, real dt);

  /**
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
   */
  real calcSmoothStep(real currentTime);

  /**
   * output rupture front, saves update time of the rupture front
   * rupture front is the first registered change in slip rates that exceeds 0.001
   */
  void saveRuptureFrontOutput(unsigned int ltsFace);

  /**
   * Save the maximal computed slip rate magnitude in peakSlipRate
   */
  void savePeakSlipRateOutput(unsigned int ltsFace);

  /**
   * Compute and store element-averaged slip to determine the magnitude of an earthquake.
   * In calc_seissol.f90 this value will be multiplied by the element surface
   * and the seismic moment is outputted once at the end of the simulation.
   * @param tmpSlip
   * @param ltsFace
   */
  void saveAverageSlipOutput(std::array<real, misc::numPaddedPoints>& tmpSlip,
                             unsigned int ltsFace);

  real deltaT[CONVERGENCE_ORDER] = {};

  dr::DRParameters& drParameters;
  ImpedancesAndEta* impAndEta;
  real mFullUpdateTime;
  // CS = coordinate system
  real (*initialStressInFaultCS)[misc::numPaddedPoints][6];
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
  real (*tractionXY)[misc::numPaddedPoints];
  real (*tractionXZ)[misc::numPaddedPoints];
  real (*imposedStatePlus)[tensor::QInterpolated::size()];
  real (*imposedStateMinus)[tensor::QInterpolated::size()];

  // be careful only for some FLs initialized:
  real* averagedSlip;
  real (*dynStressTime)[misc::numPaddedPoints];
  bool (*dynStressTimePending)[misc::numPaddedPoints];

  real (*qInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  real (*qInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_FRICTIONSOLVER_H
