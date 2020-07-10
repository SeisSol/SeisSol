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

                virtual void evaluate(lts::Base &Lts) = 0;
            };

            class FL_17 : public Base {
            public:
                virtual void hook() {}

                virtual void evaluate(lts::Base &Lts) override {
                    lts::FL_16 &ConcreteLts = dynamic_cast<lts::FL_16 &>(Lts);
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
                virtual void evaluate(lts::Base &Lts) override {
                    lts::FL_33 &ConcreteLts = dynamic_cast<lts::FL_33 &>(Lts);
                    std::cout << "computing DR for FL_33\n";
                }
            };
        }
    }
}

#endif //SEISSOL_DR_FRICTION_LAW_H
