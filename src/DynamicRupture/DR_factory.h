//
// Created by adrian on 08.07.20.
//

#ifndef SEISSOL_DR_FACTORY_H
#define SEISSOL_DR_FACTORY_H

#include <iostream>
#include <tuple>
#include <stdexcept>
#include <Initializer/DynamicRupture.h>
#include "DR_initializer_base.h"
#include "DR_friction_law.h"
#include "DR_output.h"


namespace seissol {
  namespace dr {
    namespace factory {
      class AbstractFactory;
      struct FL_2;
      struct FL_3; //aging law
      struct FL_4; //slip law
      struct FL_16;
      struct FL_33;
      struct FL_103;
      seissol::dr::factory::AbstractFactory* getFactory(Friction_law_type FrictionLawID);
    }
  }
}


using products = std::tuple<seissol::initializers::DynamicRupture*,seissol::dr::initializer::Base*, seissol::dr::fr_law::Base*, seissol::dr::output::Base*>;
class seissol::dr::factory::AbstractFactory {
public:
    virtual ~AbstractFactory() {}
    virtual products produce() = 0;
};
class seissol::dr::factory::FL_2 : public seissol::dr::factory::AbstractFactory {
    virtual products produce() override {
        return std::make_tuple(new seissol::initializers::DR_FL_2,
                               new seissol::dr::initializer::FL_2,
                               new seissol::dr::fr_law::FL_2,
                               new seissol::dr::output::FL_2);
    }
};
class seissol::dr::factory::FL_3 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_3,
                           new seissol::dr::initializer::FL_3,
                           new seissol::dr::fr_law::FL_3,
                           new seissol::dr::output::FL_3);
  }
};

class seissol::dr::factory::FL_4 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_3,
                           new seissol::dr::initializer::FL_3,
                           new seissol::dr::fr_law::FL_4,
                           new seissol::dr::output::FL_3);
  }
};

class seissol::dr::factory::FL_16 : public seissol::dr::factory::AbstractFactory {
    virtual products produce() override {
        return std::make_tuple(new seissol::initializers::DR_FL_2,
                               new seissol::dr::initializer::FL_2,
                               new seissol::dr::fr_law::FL_16,
                               new seissol::dr::output::FL_2);
    }
};

class seissol::dr::factory::FL_33 : public seissol::dr::factory::AbstractFactory {
    virtual products produce() override {
        return std::make_tuple(new seissol::initializers::DR_FL_33,
                               new seissol::dr::initializer::FL_33,
                               new seissol::dr::fr_law::FL_33,
                               new seissol::dr::output::FL_33);
    }
};

class seissol::dr::factory::FL_103 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_103,
                           new seissol::dr::initializer::FL_103,
                           new seissol::dr::fr_law::FL_103,
                           new seissol::dr::output::FL_103);
  }
};

seissol::dr::factory::AbstractFactory* seissol::dr::factory::getFactory(Friction_law_type FrictionLawID);

int temporary_main();

#endif //SEISSOL_DR_FACTORY_H
