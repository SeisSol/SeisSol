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

                virtual void tiePointers(lts::Base &Lts) = 0;

                virtual void postCompute(lts::Base &Lts) = 0;
            };

            class FL_16 : public Base {
            public:
                virtual void tiePointers(lts::Base &Lts) override {
                    lts::FL_16 &ConcreteLts = dynamic_cast<lts::FL_16 &>(Lts);
                    std::cout << "tie ptr for FL_16\n";
                }

                virtual void postCompute(lts::Base &Lts) override {
                    std::cout << "output vars for FL_16\n";
                }
            };

            class FL_33 : public Base {
            public:
                virtual void tiePointers(lts::Base &Lts) override {
                    lts::FL_33 &ConcreteLts = dynamic_cast<lts::FL_33 &>(Lts);
                    std::cout << "tie ptr for FL_33\n";
                }

                virtual void postCompute(lts::Base &Lts) override {
                    lts::FL_33 &ConcreteLts = dynamic_cast<lts::FL_33 &>(Lts);
                    std::cout << "output vars for FL_33\n";
                }
            };
        }
    }
}

#endif //SEISSOL_DR_OUTPUT_H
