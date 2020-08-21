//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_OUTPUT_H
#define SEISSOL_DR_OUTPUT_H

#include <yaml-cpp/yaml.h>

namespace seissol {
    namespace dr {
        namespace output {
            class Base;
            class FL_2;
            class FL_3;
            class FL_33;
        }
    }
}



class seissol::dr::output::Base{
public:
    virtual ~Base() {}
    void setInputParam(const YAML::Node& Param) {m_InputParam = Param;}

    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) {

        real  (*mu)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->mu);
        real  (*slip)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip);
        real  (*slip1)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip1);
        real  (*slip2)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip2);
        real  (*slipRate1)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRate1);
        real  (*slipRate2)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRate2);
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
            e_interoperability.copyFrictionOutputToFortran(ltsFace,  meshFace, mu, slip, slip1, slip2, slipRate1, slipRate2, rupture_time, peakSR, tracXY, tracXZ);
        }
    }

    virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) = 0;

protected:
  YAML::Node m_InputParam;
};

class seissol::dr::output::FL_2 : public seissol::dr::output::Base {
public:
    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) override {

        Base::tiePointers(layerData, dynRup, e_interoperability);

        seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);

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
        std::cout << "output vars for FL_2\n";
    }
};

class seissol::dr::output::FL_3 : public seissol::dr::output::Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);
    std::cout << "tie ptr for FL_33\n";
  }

  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    seissol::initializers::DR_FL_33 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 &>(DynRup);
    std::cout << "output vars for FL_33\n";
  }
};



class seissol::dr::output::FL_33 : public seissol::dr::output::Base {
public:
    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) override {
        seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);
        std::cout << "tie ptr for FL_33\n";
    }

    virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
        seissol::initializers::DR_FL_33 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 &>(DynRup);
        std::cout << "output vars for FL_33\n";
    }
};

#endif //SEISSOL_DR_OUTPUT_H
