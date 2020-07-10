//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_FRICTION_LAW_H
#define SEISSOL_DR_FRICTION_LAW_H

#include <c++/8.3.0/iostream>
#include "DR_LTS_Base.h"

namespace seissol {
    namespace dr {
        namespace fr_law {
            class Base {
            public:
                virtual ~Base() {}

                virtual void evaluate(seissol::initializers::DynamicRupture &dynRup) = 0;
            };

            class FL_2 : public Base {
            public:
                virtual void evaluate(seissol::initializers::DynamicRupture &dynRup) override {
                    seissol::initializers::DR_FL_2 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 &>(dynRup);
                    //std::cout << "computing DR for FL_2\n";
                }
            };

            class FL_17 : public Base {
            public:
                virtual void hook() {}

                virtual void evaluate(seissol::initializers::DynamicRupture &dynRup) override {
                    seissol::initializers::DR_FL_16 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_16 &>(dynRup);
                    //std::cout << "computing ";
                    hook();
                    //std::cout << " DR for FL_16\n";
                }
            };

            class FL_16 : public FL_17 {
            public:
                virtual void hook() override {
                    //std::cout << "(hook)";
                }
            };

            class FL_33 : public Base {
            public:
                virtual void evaluate(seissol::initializers::DynamicRupture &dynRup) override {
                    seissol::initializers::DR_FL_33 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 &>(dynRup);
                    std::cout << "computing DR for FL_33\n";
                }
            };
        }
    }
}

#endif //SEISSOL_DR_FRICTION_LAW_H
