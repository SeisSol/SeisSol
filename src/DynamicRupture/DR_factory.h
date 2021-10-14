//
// Created by adrian on 08.07.20.
//

#ifndef SEISSOL_DR_FACTORY_H
#define SEISSOL_DR_FACTORY_H

#include <iostream>
#include <tuple>
#include <stdexcept>
#include <Initializer/DynamicRupture.h>
#include "DynamicRupture/Initializers/Initializers.h"
#include "DynamicRupture/FrictionLaws/FrictionLaws.h"
#include "DR_output.h"

namespace seissol::dr::factory {
  using products = std::tuple<seissol::initializers::DynamicRupture*,seissol::dr::initializers::BaseDRInitializer*, seissol::dr::friction_law::BaseFrictionLaw*, seissol::dr::output::Output_Base*>;
  class AbstractFactory;
  struct Factory_FL_0; // no fault
  struct Factory_FL_2; // linear slip weakening
  struct Factory_FL_3; // aging law
  struct Factory_FL_4; // slip law
  struct Factory_FL_6; // bimaterial lsw
  struct Factory_FL_7; // velocity weakening
  struct Factory_FL_16; // linears slip weakening
  struct Factory_FL_33; // imposes slip
  struct Factory_FL_103; // rate and state
  struct Factory_FL_103_Thermal; // rate and state with thermal pressurization
  seissol::dr::factory::AbstractFactory* getFactory(dr::DRParameters *DynRupParameter);
}

class seissol::dr::factory::AbstractFactory {
public:
  virtual ~AbstractFactory() {}
  virtual products produce() = 0;
};

class seissol::dr::factory::Factory_FL_0 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

class seissol::dr::factory::Factory_FL_2 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

class seissol::dr::factory::Factory_FL_3 : public seissol::dr::factory::AbstractFactory {
  virtual products produce();
};

class seissol::dr::factory::Factory_FL_4 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

class seissol::dr::factory::Factory_FL_6 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};


class seissol::dr::factory::Factory_FL_7 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

class seissol::dr::factory::Factory_FL_16 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

class seissol::dr::factory::Factory_FL_33 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

class seissol::dr::factory::Factory_FL_103 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

class seissol::dr::factory::Factory_FL_103_Thermal : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override;
};

seissol::dr::factory::AbstractFactory* seissol::dr::factory::getFactory(dr::DRParameters *DynRupParameter);

#endif //SEISSOL_DR_FACTORY_H
