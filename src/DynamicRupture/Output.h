#ifndef SEISSOL_OUTPUT_H
#define SEISSOL_OUTPUT_H

#include <yaml-cpp/yaml.h>

#include "Parameters.h"
#include "Initializer/DynamicRupture.h"
#include "Solver/Interoperability.h"

/*
 * Currently this is only an interface to the old Fortran output writer.
 * should be replaced fully implemented C++ output writer
 */
namespace seissol::dr::output {

/**
 * abstract class, implements output that all FLs have in common
 * Currently this class only copies values computed in C++ dynamic rupture back to Fortran to have
 * it available for the output writer
 */
class OutputBase {
  protected:
  // currently, not required, but later for fully C++ implementation important
  dr::DRParameters& drParameters;

  public:
  OutputBase(dr::DRParameters& drParameters) : drParameters(drParameters){};

  virtual ~OutputBase() {}

  /**
   * Copies values from LTS tree back to Fortran data structure
   * unnecessary performance loss due to this copy back
   */
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) {
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*slip)[size] = layerData.var(dynRup->slip);
    real(*slipStrike)[size] = layerData.var(dynRup->slipStrike);
    real(*slipDip)[size] = layerData.var(dynRup->slipDip);
    real(*rupture_time)[size] = layerData.var(dynRup->ruptureTime);
    real(*peakSR)[size] = layerData.var(dynRup->peakSlipRate);
    real(*tractionXY)[size] = layerData.var(dynRup->tractionXY);
    real(*tractionXZ)[size] = layerData.var(dynRup->tractionXZ);

    DRFaceInformation* faceInformation = layerData.var(dynRup->faceInformation);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortran(ltsFace,
                                                     meshFace,
                                                     slip,
                                                     slipStrike,
                                                     slipDip,
                                                     rupture_time,
                                                     peakSR,
                                                     tractionXY,
                                                     tractionXZ);
    }
  }

  /**
   * this function has no use yet.
   */
  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) = 0;
};

/**
 * output for the no fault friction law
 */
class OutputNoFault : public OutputBase {
  public:
  using OutputBase::OutputBase;

  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};

/**
 * Derived class for the linear slip weakening friction law (also for the one with forced rupture
 * time)
 */
class OutputLinearSlipWeakening : public OutputBase {
  public:
  using OutputBase::OutputBase;

  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);

    auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*dynStressTime)[size] = layerData.var(concreteLts->dynStressTime);
    real(*slipRateStrike)[size] = layerData.var(concreteLts->slipRateStrike);
    real(*slipRateDip)[size] = layerData.var(concreteLts->slipRateDip);
    real(*mu)[size] = layerData.var(concreteLts->mu);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averagedSlip, dynStressTime, slipRateStrike, slipRateDip, mu);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};

class OutputImposedSlipRates : public OutputBase {
  public:
  using OutputBase::OutputBase;

  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};

/**
 *  output for rate and state friction laws
 */
class OutputRateAndState : public OutputBase {
  public:
  using OutputBase::OutputBase;

  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);

    auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*dynStressTime)[size] = layerData.var(concreteLts->dynStressTime);
    real(*slipRateStrike)[size] = layerData.var(concreteLts->slipRateStrike);
    real(*slipRateDip)[size] = layerData.var(concreteLts->slipRateDip);
    real(*mu)[size] = layerData.var(concreteLts->mu);
    real(*stateVar)[size] = layerData.var(concreteLts->stateVariable);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averagedSlip, dynStressTime, slipRateStrike, slipRateDip, mu);
      e_interoperability.copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVar);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};

/**
 * extends output for regularized strength
 */
class OutputLinearSlipWeakeningBimaterial : public OutputLinearSlipWeakening {
  public:
  using OutputLinearSlipWeakening::OutputLinearSlipWeakening;

  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputLinearSlipWeakening::tiePointers(layerData, dynRup, e_interoperability);

    auto concreteLts =
        dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*regularisedStrength)[size] = layerData.var(concreteLts->regularisedStrength);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranStrength(
          ltsFace, meshFace, regularisedStrength);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};

/**
 * output for rate and state Friction laws with additional thermal pressurisation
 */
class OutputRateAndStateThermalPressurisation : public OutputRateAndState {
  public:
  using OutputRateAndState::OutputRateAndState;
  using OutputRateAndState::postCompute;
  using OutputRateAndState::tiePointers;
};
} // namespace seissol::dr::output

#endif // SEISSOL_OUTPUT_H
