//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_INITIALIZER_BASE_H
#define SEISSOL_DR_INITIALIZER_BASE_H

#include <c++/8.3.0/iostream>
#include "DR_LTS_Base.h"

namespace seissol {
    namespace dr {
        namespace initializer {
            class Base {
            public:
                virtual ~Base() {}

                virtual void initializeFrictionMatrices(lts::Base *Lts, initializers::LTSTree* dynRupTree) {
                    // Note: actually the vars. belong to DrLtsTree
                    Lts->Slip = 0;

                    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
                        real *mu = it->var(Lts->mu);
                        for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
                            mu[ltsFace] = 1;
                        }
                    }
                    std::cout << "init DR for Base\n";
                }
            };

            class FL_16 : public Base {
            public:
                virtual void initializeFrictionMatrices(lts::Base *Lts, initializers::LTSTree* dynRupTree) override {
                    Base::initializeFrictionMatrices(Lts, dynRupTree);
                    lts::FL_16 *ConcreteLts = dynamic_cast<lts::FL_16 *>(Lts);
                    ConcreteLts->InitStress = 2; // Note: actually the vars. belong to DrLtsTree
                    std::cout << "init DR for FL_16\n";
                }
            };

            class FL_33 : public Base {
            public:
                virtual void initializeFrictionMatrices(lts::Base *Lts, initializers::LTSTree* dynRupTree) override {
                    Base::initializeFrictionMatrices(Lts, dynRupTree);
                    lts::FL_33 *ConcreteLts = dynamic_cast<lts::FL_33 *>(Lts);
                    ConcreteLts->InitPressure = 3; // Note: actually the vars. belong to DrLtsTree
                    std::cout << "init DR for FL_33\n";
                }
            };
        }
    }
}
#endif //SEISSOL_DR_INITIALIZER_BASE_H
