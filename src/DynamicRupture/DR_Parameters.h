//
// Created by adrian on 03.09.20.
//

#ifndef SEISSOL_DR_PARAMETERS_H
#define SEISSOL_DR_PARAMETERS_H

#include <yaml-cpp/yaml.h>
#include "Initializer/InputAux.hpp"

namespace seissol {
  namespace dr {
    struct DrParameterT;
  }
}

struct seissol::dr::DrParameterT {
  int OutputPointType{3};
  int SlipRateOutputType{1};
  int FrictionLawType{0};
  int BackgroundType{0};
  bool IsRfOutputOn{false};
  bool IsDsOutputOn{false};
  bool IsMagnitudeOutputOn{false};
  bool IsEnergyRateOutputOn{false};
  bool IsGpWiseOutput{false};
  bool IsTermalPressureOn{false};
  int EnergyRatePrintTimeInterval{1};
  bool IsInstaHealingOn{false};
  real t_0{0.0};
  real rs_f0{0.0};
  real rs_b{0.0};
  real rs_sr0{0.0};
  real mu_w{0.0};


  void setAllInputParam(const YAML::Node& Params) {
    using namespace initializers;

    const YAML::Node& DrParams = Params["dynamicrupture"];
    OutputPointType = getParamIfExists(DrParams, "outputpointtype", 3);
    SlipRateOutputType = getParamIfExists(DrParams, "sliprateoutputtype", 1);
    FrictionLawType = getParamIfExists(DrParams, "fl", 0);
    BackgroundType = getParamIfExists(DrParams, "backgroundtype", 0);
    IsRfOutputOn = getParamIfExists(DrParams, "rf_output_on", false);
    IsDsOutputOn = getParamIfExists(DrParams, "ds_output_on", false);
    IsMagnitudeOutputOn = getParamIfExists(DrParams, "magnitude_output_on", false);
    IsEnergyRateOutputOn = getParamIfExists(DrParams, "energy_rate_output_on", false);
    IsGpWiseOutput = getParamIfExists(DrParams, "gpwise", false);
    IsTermalPressureOn = getParamIfExists(DrParams, "thermalpress", false);
    BackgroundType = getParamIfExists(DrParams, "energy_rate_printtimeinterval", 1);
    IsInstaHealingOn = getParamIfExists(DrParams, "inst_healing", false);
    t_0 = getParamIfExists(DrParams, "t_0", 0.0);
    rs_f0 = getParamIfExists(DrParams, "rs_f0", 0.0);
    rs_b = getParamIfExists(DrParams, "rs_b", 0.0);
    rs_sr0 = getParamIfExists(DrParams, "rs_sr0", 0.0);
    mu_w = getParamIfExists(DrParams, "mu_w", 0.0);
  }
};

#endif //SEISSOL_DR_PARAMETERS_H
