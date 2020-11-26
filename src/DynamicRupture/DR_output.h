//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_OUTPUT_H
#define SEISSOL_DR_OUTPUT_H

#include <yaml-cpp/yaml.h>

namespace seissol {
    namespace dr {
        namespace output {
            class Output_Base;
            class Output_FL_0;
            class Output_FL_2;
            class Output_FL_3;
            class Output_FL_6;
            class Output_FL_33;
            class Output_FL_103;
        }
    }
}


//TODO: RF-files (e.g. tpv5-RF-00000-TID-00) are missing in the C++ output
class seissol::dr::output::Output_Base{
protected:
  dr::DrParameterT *m_Params;
public:
    virtual ~Output_Base() {}

  //set the parameters from .par file with yaml to this class attributes.
  void setInputParam(dr::DrParameterT *DynRupParameter) {
    m_Params = DynRupParameter;
  }

    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) {

        real  (*mu)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->mu);
        real  (*slip)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip);
        real  (*slip1)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip1);
        real  (*slip2)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip2);
        real  (*slipRateStrike)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRateStrike);
        real  (*slipRateDip)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRateDip);
        real  (*rupture_time)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->rupture_time);
        real  (*peakSR)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->peakSR);
        real  (*tracXY)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->tracXY);
        real  (*tracXZ)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->tracXZ);

        DRFaceInformation*                    faceInformation = layerData.var(dynRup->faceInformation);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
            unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
            e_interoperability.copyFrictionOutputToFortran(ltsFace, meshFace, mu, slip, slip1, slip2, slipRateStrike, slipRateDip, rupture_time, peakSR, tracXY, tracXZ);
        }
    }

    virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) = 0;



};

class seissol::dr::output::Output_FL_0 : public seissol::dr::output::Output_Base {
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
  }
  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    //do nothing
  }
};

class seissol::dr::output::Output_FL_2 : public seissol::dr::output::Output_Base {
public:
    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) override {

        Output_Base::tiePointers(layerData, dynRup, e_interoperability);

        seissol::initializers::DR_linear *ConcreteLts = dynamic_cast<seissol::initializers::DR_linear *>(dynRup);

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

class seissol::dr::output::Output_FL_3 : public seissol::dr::output::Output_Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
    //seissol::initializers::DR_FL_3 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3 *>(dynRup);
    std::cout << "tie ptr for Init_FL_3\n";
  }

  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    //seissol::initializers::DR_FL_3 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3 &>(DynRup);
    std::cout << "output vars for Init_FL_3\n";
  }
};



class seissol::dr::output::Output_FL_33 : public seissol::dr::output::Output_Base {
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

class seissol::dr::output::Output_FL_103 : public seissol::dr::output::Output_Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
    seissol::initializers::DR_FL_103 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103 *>(dynRup);
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

class seissol::dr::output::Output_FL_6 : public seissol::dr::output::Output_Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    Output_Base::tiePointers(layerData, dynRup, e_interoperability);
    seissol::initializers::DR_FL_6 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_6 *>(dynRup);

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

#endif //SEISSOL_DR_OUTPUT_H
