//
// Created by adrian on 08.07.20.
//

#ifndef SEISSOL_DR_FACTORY_H
#define SEISSOL_DR_FACTORY_H

#include <iostream>
#include <tuple>
#include <stdexcept>
#include "DR_LTS_Base.h"
#include "DR_initializer_base.h"
#include "DR_friction_law.h"
#include "DR_output.h"


namespace seissol {
    namespace dr {
        namespace factory {
            using products = std::tuple<lts::Base*, initializer::Base*, fr_law::Base*, output::Base*>;
            class AbstractFactory {
            public:
                virtual ~AbstractFactory() {}
                virtual products produce() = 0;
            };
            class FL_16 : public AbstractFactory {
                virtual products produce() override {
                    return std::make_tuple(new lts::FL_16,
                                           new initializer::FL_16,
                                           new fr_law::FL_16,
                                           new output::FL_16);
                }
            };
            class FL_17 : public AbstractFactory {
                virtual products produce() override {
                    return std::make_tuple(new lts::FL_16,
                                           new initializer::FL_16,
                                           new fr_law::FL_17,
                                           new output::FL_16);
                }
            };
            class FL_33 : public AbstractFactory {
                virtual products produce() override {
                    return std::make_tuple(new lts::FL_33,
                                           new initializer::FL_33,
                                           new fr_law::FL_33,
                                           new output::FL_33);
                }
            };
            AbstractFactory* getFactory(Friction_law_type FrictionLawID);
        }
    }
}

int temporary_main();

#endif //SEISSOL_DR_FACTORY_H
