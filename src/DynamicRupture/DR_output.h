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

                virtual void tiePointers(seissol::initializers::DynamicRupture &DynRup) = 0;

                virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) = 0;
            };

            class FL_2 : public Base {
            public:
                virtual void tiePointers(seissol::initializers::DynamicRupture &DynRup) override {
                    seissol::initializers::DR_FL_2 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 &>(DynRup);
                    //std::cout << "tie ptr for FL_2\n";
                }

                virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
                    std::cout << "output vars for FL_2\n";
                }
            };

            class FL_16 : public Base {
            public:
                virtual void tiePointers(seissol::initializers::DynamicRupture &DynRup) override {
                    seissol::initializers::DR_FL_16 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_16 &>(DynRup);
                    std::cout << "tie ptr for FL_16\n";
                }

                virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
                    std::cout << "output vars for FL_16\n";
                }
            };

            class FL_33 : public Base {
            public:
                virtual void tiePointers(seissol::initializers::DynamicRupture &DynRup) override {
                    seissol::initializers::DR_FL_33 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 &>(DynRup);
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
