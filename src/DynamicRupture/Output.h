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
class OutputBase;                // abstract class, implements output that all FLs have in common
class OutputNoFault;             // No SolverNoFaultFL0
class OutputLinearSlipWeakening; // for linear slip laws FL2, FL16
class OutputLinearSlipWeakeningBimaterial;     // extends output for Strength output
class OutputRateAndState;                      // dummy class, not implemented, could be replaced by
                                               // OutputRateAndStateFL103
class OutputRateAndStateFastVelocityWeakening; // output for rate and state Friction laws
class OutputImposedSlipRates;                  // SolverImposedSlipRatesFL33
} // namespace seissol::dr::output

/*
 * Currently this class only copies values computed in C++ dynamic rupture back to Fortran to have
 * it available for the output writer
 */
class seissol::dr::output::OutputBase {
  protected:
  // currently, not required, but later for fully C++ implementation important
  dr::DRParameters& drParameters;

  public:
  OutputBase(dr::DRParameters& drParameters) : drParameters(drParameters){};
  virtual ~OutputBase() {}

  /*
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

    real(*iniBulkXX)[size] = layerData.var(dynRup->iniBulkXX);
    real(*iniBulkYY)[size] = layerData.var(dynRup->iniBulkYY);
    real(*iniBulkZZ)[size] = layerData.var(dynRup->iniBulkZZ);
    real(*iniShearXY)[size] = layerData.var(dynRup->iniShearXY);
    real(*iniShearYZ)[size] = layerData.var(dynRup->iniShearYZ);
    real(*iniShearXZ)[size] = layerData.var(dynRup->iniShearXZ);

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

  /*
   * this function has no use yet.
   */
  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) = 0;
};

class seissol::dr::output::OutputNoFault : public seissol::dr::output::OutputBase {
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

class seissol::dr::output::OutputLinearSlipWeakening : public seissol::dr::output::OutputBase {
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
    std::cout << "output vars for Factory_FL_2\n";
  }
};

class seissol::dr::output::OutputImposedSlipRates : public seissol::dr::output::OutputBase {
  public:
  using OutputBase::OutputBase;
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // seissol::initializers::DR_FL_33 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33
    // &>(DynRup);
    std::cout << "output vars for Init_FL_33\n";
  }
};

class seissol::dr::output::OutputRateAndStateFastVelocityWeakening
    : public seissol::dr::output::OutputBase {
  public:
  using OutputBase::OutputBase;
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
    auto concreteLts =
        dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);
    // std::cout << "tie ptr for Init_FL_103\n";

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*dynStressTime)[size] = layerData.var(concreteLts->dynStressTime);
    real(*slipRateStrike)[size] = layerData.var(concreteLts->slipRateStrike);
    real(*slipRateDip)[size] = layerData.var(concreteLts->slipRateDip);
    real(*mu)[size] = layerData.var(concreteLts->mu);
    real(*stateVar)[size] = layerData.var(concreteLts->stateVariable);
    real(*initialStressInFaultCS)[size][6] = layerData.var(dynRup->initialStressInFaultCS);

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
    // seissol::initializers::DR_FL_103 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103
    // &>(DynRup);
    std::cout << "output vars for Init_FL_103\n";
  }
};

// output_Strength

class seissol::dr::output::OutputLinearSlipWeakeningBimaterial
    : public seissol::dr::output::OutputBase {
  public:
  using OutputBase::OutputBase;
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
    auto concreteLts =
        dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*dynStressTime)[size] = layerData.var(concreteLts->dynStressTime);
    real(*slipRateStrike)[size] = layerData.var(concreteLts->slipRateStrike);
    real(*slipRateDip)[size] = layerData.var(concreteLts->slipRateDip);
    real(*mu)[size] = layerData.var(concreteLts->mu);
    real(*regularisedStrength)[size] = layerData.var(concreteLts->regularisedStrength);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averagedSlip, dynStressTime, slipRateStrike, slipRateDip, mu);
      e_interoperability.copyFrictionOutputToFortranStrength(
          ltsFace, meshFace, regularisedStrength);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // seissol::initializers::DR_FL_103 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103
    // &>(DynRup);
    std::cout << "output vars for Init_FL_6\n";
  }
};

/*
 * Not implemented
 */
class seissol::dr::output::OutputRateAndState : public OutputBase {
  public:
  using OutputBase::OutputBase;
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
    // seissol::initializers::DR_FL_3 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3
    // *>(dynRup);
    std::cout << "tie ptr for Init_FL_3 (not implemented)\n";
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // seissol::initializers::DR_FL_3 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3
    // &>(DynRup);
    std::cout << "output vars for Init_FL_3\n";
  }
};

#endif // SEISSOL_OUTPUT_H
