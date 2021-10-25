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
class OutputBase;                   // abstract class, implements output that all FLs have in common
class OutputNoFaultFL0;             // No SolverNoFaultFL0
class OutputLinearSlipWeakeningFL2; // for linear slip laws FL2, FL16 and FL6
class OutputRateAndStateFL3;        // dummy class, not implemented, could be replaced by
                                    // OutputRateAndStateFL103
class OutputLinearBimaterialFL6;    // extends output for Strength output
class OutputImposedSlipRatesFL33;   // SolverImposedSlipRatesFL33
class OutputRateAndStateFL103;      // output for rate and state Friction laws
} // namespace seissol::dr::output

/*
 * Currently this class only copies values computed in C++ dynamic rupture back to Fortran to have
 * it available for the output writer
 */
class seissol::dr::output::OutputBase {
  protected:
  seissol::dr::DRParameters*
      m_Params; // currently, not required, but later for fully C++ implementation important
  public:
  virtual ~OutputBase() {}

  // set the parameters from .par file with yaml to this class attributes.
  void setInputParam(dr::DRParameters* DynRupParameter) { m_Params = DynRupParameter; }

  /*
   * Copies values from LTS tree back to Fortran data structure
   * unnecessary performance loss due to this copy back
   */
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) {
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*mu)[size] = layerData.var(dynRup->mu);
    real(*slip)[size] = layerData.var(dynRup->slip);
    real(*slipStrike)[size] = layerData.var(dynRup->slipStrike);
    real(*slipDip)[size] = layerData.var(dynRup->slipDip);
    real(*slipRateStrike)[size] = layerData.var(dynRup->slipRateStrike);
    real(*slipRateDip)[size] = layerData.var(dynRup->slipRateDip);
    real(*rupture_time)[size] = layerData.var(dynRup->rupture_time);
    real(*peakSR)[size] = layerData.var(dynRup->peakSR);
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
                                                     mu,
                                                     slip,
                                                     slipStrike,
                                                     slipDip,
                                                     slipRateStrike,
                                                     slipRateDip,
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

class seissol::dr::output::OutputNoFaultFL0 : public seissol::dr::output::OutputBase {
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
  }
  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};

class seissol::dr::output::OutputLinearSlipWeakeningFL2 : public seissol::dr::output::OutputBase {
  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {

    OutputBase::tiePointers(layerData, dynRup, e_interoperability);

    seissol::initializers::LTS_LinearSlipWeakeningFL2* ConcreteLts =
        dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningFL2*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(ConcreteLts->faceInformation);
    real* averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*dynStress_time)[size] = layerData.var(ConcreteLts->dynStress_time);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averaged_Slip, dynStress_time);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    std::cout << "output vars for Factory_FL_2\n";
  }
};

class seissol::dr::output::OutputImposedSlipRatesFL33 : public seissol::dr::output::OutputBase {
  public:
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

class seissol::dr::output::OutputRateAndStateFL103 : public seissol::dr::output::OutputBase {
  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
    seissol::initializers::LTS_RateAndStateFL103* ConcreteLts =
        dynamic_cast<seissol::initializers::LTS_RateAndStateFL103*>(dynRup);
    // std::cout << "tie ptr for Init_FL_103\n";

    DRFaceInformation* faceInformation = layerData.var(ConcreteLts->faceInformation);
    real* averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*dynStress_time)[size] = layerData.var(ConcreteLts->dynStress_time);
    real(*stateVar)[size] = layerData.var(ConcreteLts->stateVar);
    real(*initialStressInFaultCS)[size][6] = layerData.var(dynRup->initialStressInFaultCS);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averaged_Slip, dynStress_time);
      e_interoperability.copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVar);
      e_interoperability.copyFrictionOutputToFortranInitialStressInFaultCS(
          ltsFace, meshFace, initialStressInFaultCS);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // seissol::initializers::DR_FL_103 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103
    // &>(DynRup);
    std::cout << "output vars for Init_FL_103\n";
  }
};

// output_Strength

class seissol::dr::output::OutputLinearBimaterialFL6 : public seissol::dr::output::OutputBase {
  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    OutputBase::tiePointers(layerData, dynRup, e_interoperability);
    seissol::initializers::LTS_LinearBimaterialFL6* ConcreteLts =
        dynamic_cast<seissol::initializers::LTS_LinearBimaterialFL6*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(ConcreteLts->faceInformation);
    real* averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*dynStress_time)[size] = layerData.var(ConcreteLts->dynStress_time);
    real(*strength)[size] = layerData.var(ConcreteLts->strengthData);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averaged_Slip, dynStress_time);
      e_interoperability.copyFrictionOutputToFortranStrength(ltsFace, meshFace, strength);
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
class seissol::dr::output::OutputRateAndStateFL3 : public seissol::dr::output::OutputBase {
  public:
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
