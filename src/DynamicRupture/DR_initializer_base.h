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

                virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup, initializers::LTSTree* dynRupTree) {
                    // Note: actually the vars. belong to DrLtsTree

                    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
                        real *mu = it->var(dynRup->mu);
                        for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
                            mu[ltsFace] = 1;
                        }
                    }
                    std::cout << "init DR for Base\n";
                }
            };
            class FL_2 : public Base {
            public:
                virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup, initializers::LTSTree* dynRupTree) override {
                    Base::initializeFrictionMatrices(dynRup, dynRupTree);
                    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);
                    ConcreteLts->InitStress = 2; // Note: actually the vars. belong to DrLtsTree
                    std::cout << "init DR for FL_2\n";
                }
            };

            class FL_16 : public Base {
            public:
                virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup, initializers::LTSTree* dynRupTree) override {
                    Base::initializeFrictionMatrices(dynRup, dynRupTree);
                    seissol::initializers::DR_FL_16 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_16 *>(dynRup);
                    ConcreteLts->InitStress = 2; // Note: actually the vars. belong to DrLtsTree
                    std::cout << "init DR for FL_16\n";
                }
            };

            class FL_33 : public Base {
            public:
                virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup, initializers::LTSTree* dynRupTree) override {
                    Base::initializeFrictionMatrices(dynRup, dynRupTree);
                    seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);
                    ConcreteLts->InitPressure = 3; // Note: actually the vars. belong to DrLtsTree
                    std::cout << "init DR for FL_33\n";
                }
            };
        }
    }
}
#endif //SEISSOL_DR_INITIALIZER_BASE_H
