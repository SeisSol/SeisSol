//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_OUTPUT_H
#define SEISSOL_DR_OUTPUT_H

#include "DR_LTS_Base.h"

namespace seissol {
    namespace dr {
        namespace output {
            class Base {
            public:
                virtual ~Base() {}

                virtual void tiePointers(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        seissol::Interoperability &e_interoperability) = 0;

                virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) = 0;
            };

            class FL_2 : public Base {
            public:
                virtual void tiePointers(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        seissol::Interoperability &e_interoperability) override {

                    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);
                    //std::cout << "tie ptr for FL_2\n";

                    DRFaceInformation*                    faceInformation                                                   = layerData.var(ConcreteLts->faceInformation);
                    real  (*mu)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->mu);
                    real  (*slip)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->slip);
                    real  (*slip1)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->slip1);
                    real  (*slip2)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->slip2);
                    real  (*rupture_time)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->rupture_time);
                    real  (*peakSR)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->peakSR);
                    real  (*dynStress_time)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->dynStress_time);
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
                        unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
                        e_interoperability.copyFrictionOutputToFortran(ltsFace,  meshFace,
                                 mu, slip, slip1, slip2, rupture_time, peakSR, dynStress_time);
                    }
                }

                virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
                    std::cout << "output vars for FL_2\n";
                }
            };

            class FL_16 : public Base {
            public:
                virtual void tiePointers(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        seissol::Interoperability &e_interoperability) override {
                    seissol::initializers::DR_FL_16 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_16 *>(dynRup);
                    std::cout << "tie ptr for FL_16\n";
                }

                virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
                    std::cout << "output vars for FL_16\n";
                }
            };

            class FL_33 : public Base {
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
        }
    }
}

#endif //SEISSOL_DR_OUTPUT_H
