#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/DR_Parameters.h"
#include "DynamicRupture/DR_math.h"
#include "Initializer/DynamicRupture.h"
#include "Kernels/DynamicRupture.h"

namespace seissol::dr::friction_law {
class BaseFrictionLaw;
}

// Base class, has implementations of methods that are used by each friction law
class seissol::dr::friction_law::BaseFrictionLaw {

  public:
  /*
   * Destructor, if memory is allocated in this class, deallocate it here.
   */
  virtual ~BaseFrictionLaw() {}

  /*
   * Set the parameters from .par file with yaml to this class attributes.
   * This function is called at initialisation time. Could be extended to initialize more parameters
   * if needed.
   */
  void setInputParam(dr::DRParameters* DynRupParameter) { m_Params = DynRupParameter; }

  protected:
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0]; // DISC%Galerkin%nBndGP
  static constexpr int numOfPointsPadded =
      init::QInterpolated::Stop[0]; // number of points padded to next dividable number by four
  // YAML::Node m_InputParam;
  dr::DRParameters* m_Params;
  ImpedancesAndEta* impAndEta;
  real m_fullUpdateTime;
  real deltaT[CONVERGENCE_ORDER] = {};
  real (*initialStressInFaultCS)[numOfPointsPadded][6];
  real (*cohesion)[numOfPointsPadded];
  real (*mu)[numOfPointsPadded];
  real (*slip)[numOfPointsPadded];
  real (*slipStrike)[numOfPointsPadded];
  real (*slipDip)[numOfPointsPadded];
  real (*SlipRateMagnitude)[numOfPointsPadded];
  real (*slipRateStrike)[numOfPointsPadded];
  real (*slipRateDip)[numOfPointsPadded];
  real (*rupture_time)[numOfPointsPadded];
  bool (*RF)[numOfPointsPadded];
  real (*peakSR)[numOfPointsPadded];
  real (*tractionXY)[numOfPointsPadded];
  real (*tractionXZ)[numOfPointsPadded];
  real (*imposedStatePlus)[tensor::QInterpolated::size()];
  real (*imposedStateMinus)[tensor::QInterpolated::size()];

  // be careful only for some FLs initialized:
  real* averaged_Slip;

  /*
   * Struct that contains all input stresses and output stresses
   * IN: NormalStressGP, XYStressGP, XZStressGP (Godunov stresses computed by
   * precomputeStressFromQInterpolated) OUT: XYTractionResultGP, XZTractionResultGP and
   * NormalStressGP (used to compute resulting +/- sided stress results by
   * postcomputeImposedStateFromNewStress)
   */
  struct FaultStresses {
    real XYTractionResultGP[CONVERGENCE_ORDER][numOfPointsPadded] = {
        {}}; // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
    real XZTractionResultGP[CONVERGENCE_ORDER][numOfPointsPadded] = {
        {}}; // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
    real NormalStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
    real XYStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
    real XZStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
  };

  /*
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  virtual void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                  seissol::initializers::DynamicRupture* dynRup,
                                  real fullUpdateTime);
  /*
   * output:
   * NorStressGP, XYStressGP, XZStressGP
   *
   * input:
   * QInterpolatedPlus, QInterpolatedMinus, eta_p, Zp, Zp_neig, eta_s, Zs, Zs_neig
   *
   * Calculate godunov state from jump of plus and minus side
   * using equations (A2) from Pelites et al. 2014
   * Definiton of eta and impedance Z are found in dissertation of Carsten Uphoff
   */
  virtual void precomputeStressFromQInterpolated(
      FaultStresses& faultStresses,
      real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      unsigned int ltsFace);
  /*
   * Output: imposedStatePlus, imposedStateMinus
   *
   * Integrate over all Time points with the time weights and calculate the traction vor each side
   * according to Carsten Uphoff Thesis: EQ.: 4.60 IN: NormalStressGP, XYTractionResultGP,
   * XZTractionResultGP OUT: imposedStatePlus, imposedStateMinus
   */
  void postcomputeImposedStateFromNewStress(
      real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      const FaultStresses& faultStresses,
      double timeWeights[CONVERGENCE_ORDER],
      unsigned int ltsFace);

  /*
   * Function from NucleationFunctions_mod.f90
   */
  real Calc_SmoothStepIncrement(real current_time, real dt);

  /*
   * Function from NucleationFunctions_mod.f90
   */
  real Calc_SmoothStep(real current_time);

  /*
   * output rupture front, saves update time of the rupture front
   * rupture front is the first registered change in slip rates that exceeds 0.001
   */
  void saveRuptureFrontOutput(unsigned int ltsFace);

  /*
   * save the maximal computed slip rate magnitude in peakSR
   */
  void savePeakSlipRateOutput(unsigned int ltsFace);

  //---compute and store slip to determine the magnitude of an earthquake ---
  //    to this end, here the slip is computed and averaged per element
  //    in calc_seissol.f90 this value will be multiplied by the element surface
  //    and an output happened once at the end of the simulation
  void saveAverageSlipOutput(std::array<real, numOfPointsPadded>& tmpSlip, unsigned int ltsFace);

  public:
  /*
   * evaluates the current friction model
   * Friction laws (child classes) implement this function
   */
  virtual void
      evaluate(seissol::initializers::Layer& layerData,
               seissol::initializers::DynamicRupture* dynRup,
               real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real fullUpdateTime,
               double timeWeights[CONVERGENCE_ORDER]) = 0;

  /*
   * compute the DeltaT from the current timePoints
   * call this function before evaluate to set the correct DeltaT
   */
  void computeDeltaT(double timePoints[CONVERGENCE_ORDER]);
};

#endif // SEISSOL_BASEFRICTIONLAW_H
