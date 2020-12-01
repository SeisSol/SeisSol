//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_OUTPUT_H
#define SEISSOL_DR_OUTPUT_H

#include <yaml-cpp/yaml.h>



/*
 * Currently this is only an interface to the old Fortran output writer.
 * should be replaced fully implemented C++ output writer
 */
namespace seissol {
    namespace dr {
        namespace output {
            class Output_Base;      //abstrac class, implements output that all FLs have in common
            class Output_NoFaultFL0;  //No SolverNoFaultFL0
            class Output_LinearSlipWeakeningFL2;  //for linear slip laws FL2, FL16 and FL6
            class Output_RateAndStateFL3; //dummy class, not implemented, could be replaced by Output_RateAndStateFL103
            class Output_LinearBimaterialFL6; //extends output for Strength output
            class Output_ImposedSlipRatesFL33; //SolverImposedSlipRatesFL33
            class Output_RateAndStateFL103;   //output for rate and state Friction laws
        }
    }
}

/*
 * Currently this class only copies values computed in C++ dynamic rupture back to Fortran to have it availabe for the output writer
 */
class seissol::dr::output::Output_Base{
protected:
  dr::DrParameterT *m_Params; //currently not required, but later for fully C++ implementation important
public:
    virtual ~Output_Base() {}

  //set the parameters from .par file with yaml to this class attributes.
  void setInputParam(dr::DrParameterT *DynRupParameter) {
    m_Params = DynRupParameter;
  }

  /*
   * Copies values from LTS tree back to Fortran data structure
   * unnecessary performance loss due to this copy back
   */
    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) {

        real  (*mu)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->mu);
        real  (*slip)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip);
        real  (*slipStrike)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipStrike);
        real  (*slipDip)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipDip);
        real  (*slipRateStrike)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRateStrike);
        real  (*slipRateDip)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRateDip);
        real  (*rupture_time)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->rupture_time);
        real  (*peakSR)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->peakSR);
        real  (*tractionXY)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->tractionXY);
        real  (*tractionXZ)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->tractionXZ);

        DRFaceInformation*                    faceInformation = layerData.var(dynRup->faceInformation);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
            unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
            e_interoperability.copyFrictionOutputToFortran(ltsFace, meshFace, mu, slip, slipStrike, slipDip, slipRateStrike, slipRateDip, rupture_time, peakSR, tractionXY, tractionXZ);
        }
    }

    /*
     * this function has no use yet.
     */
    virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) = 0;



};

class seissol::dr::output::Output_NoFaultFL0 : public seissol::dr::output::Output_Base {
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
  }
  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    //do nothing
  }
};

class seissol::dr::output::Output_LinearSlipWeakeningFL2 : public seissol::dr::output::Output_Base {
public:
    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) override {

        Output_Base::tiePointers(layerData, dynRup, e_interoperability);

        seissol::initializers::LTS_LinearSlipWeakeningFL2 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningFL2 *>(dynRup);

        DRFaceInformation*                    faceInformation = layerData.var(ConcreteLts->faceInformation);
        real  *averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
        real  (*dynStress_time)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->dynStress_time);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
            unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
            e_interoperability.copyFrictionOutputToFortranFL2(ltsFace,  meshFace,
                    averaged_Slip,  dynStress_time);
        }
    }

    virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
        std::cout << "output vars for Factory_FL_2\n";
    }
};




class seissol::dr::output::Output_ImposedSlipRatesFL33 : public seissol::dr::output::Output_Base {
public:
    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) override {
        Output_Base::tiePointers(layerData, dynRup, e_interoperability);
    }

    virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
        //seissol::initializers::DR_FL_33 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 &>(DynRup);
        std::cout << "output vars for Init_FL_33\n";
    }
};

class seissol::dr::output::Output_RateAndStateFL103 : public seissol::dr::output::Output_Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
    seissol::initializers::LTS_RateAndStateFL103 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_RateAndStateFL103 *>(dynRup);
    //std::cout << "tie ptr for Init_FL_103\n";


    DRFaceInformation*                    faceInformation = layerData.var(ConcreteLts->faceInformation);
    real  *averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
    real  (*dynStress_time)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->dynStress_time);
    real  (*stateVar)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->stateVar);
    real  (*initialStressInFaultCS)[init::QInterpolated::Stop[0]][6] = layerData.var(dynRup->initialStressInFaultCS);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(ltsFace,  meshFace, averaged_Slip,  dynStress_time);
      e_interoperability.copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVar);
      e_interoperability.copyFrictionOutputToFortranInitialStressInFaultCS(ltsFace, meshFace, initialStressInFaultCS);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    //seissol::initializers::DR_FL_103 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103 &>(DynRup);
    std::cout << "output vars for Init_FL_103\n";
  }
};


//output_Strength

class seissol::dr::output::Output_LinearBimaterialFL6 : public seissol::dr::output::Output_Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
    seissol::initializers::LTS_LinearBimaterialFL6 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_LinearBimaterialFL6 *>(dynRup);

    DRFaceInformation*                    faceInformation = layerData.var(ConcreteLts->faceInformation);
    real  *averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
    real  (*dynStress_time)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->dynStress_time);
    real  (*strength)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->strengthData);


#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(ltsFace,  meshFace, averaged_Slip,  dynStress_time);
      e_interoperability.copyFrictionOutputToFortranStrength(ltsFace, meshFace, strength);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    //seissol::initializers::DR_FL_103 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103 &>(DynRup);
    std::cout << "output vars for Init_FL_6\n";
  }
};

/*
 * Not implemented
 */
class seissol::dr::output::Output_RateAndStateFL3 : public seissol::dr::output::Output_Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
    //seissol::initializers::DR_FL_3 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3 *>(dynRup);
    std::cout << "tie ptr for Init_FL_3 (not implemented)\n";
  }

  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    //seissol::initializers::DR_FL_3 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3 &>(DynRup);
    std::cout << "output vars for Init_FL_3\n";
  }
};


#endif //SEISSOL_DR_OUTPUT_H
