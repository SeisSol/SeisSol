//
// Created by adrian on 09.07.20.
//


#ifndef SEISSOL_DRBASE_H
#define SEISSOL_DRBASE_H

#include <c++/8.3.0/iostream>
#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <generated_code/tensor.h>

/*
 * file can be deleted
 *
namespace seissol {
    namespace dr {
        namespace lts {
            struct Base {
                virtual ~Base() {}
                int Slip{};
                initializers::Variable<real>                                                    mu;

                virtual void addVars(initializers::LTSTree& tree) {
                    initializers::LayerMask mask = initializers::LayerMask(Ghost);
                    tree.addVar(      mu,                             mask,                 1,      seissol::memory::Standard );

                    std::cout << "add the base variables to the tree\n";
                }
            };

            struct FL_16 : public Base {
                int InitStress{};

                virtual void addVars(initializers::LTSTree& tree) {
                    Base::addVars(tree);
                    std::cout << "add additional variables for FL_16\n";
                }
            };

            struct FL_33 : public Base {
                int InitPressure{};

                virtual void addVars(initializers::LTSTree& tree) {
                    Base::addVars(tree);
                    std::cout << "add additional variables for FL_33\n";
                }
            };
        }
    }
}
*/

#endif //SEISSOL_DRBASE_H
