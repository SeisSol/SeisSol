#ifndef SEISSOL_PARAMETERS_H
#define SEISSOL_PARAMETERS_H

#include <yaml-cpp/yaml.h>

#include "Initializer/InputAux.hpp"
#include "Kernels/precision.hpp"
#include "Typedefs.hpp"

namespace seissol::dr {
struct DRParameters;
} // namespace seissol::dr

/*
 * Saves all dynamic rupture parameter read from parameter.par file
 * if values are not defined they are set to an initial value (mostly 0)
 */
struct seissol::dr::DRParameters {
  static constexpr unsigned int TP_grid_nz = 60;
  int outputPointType{3};
  int slipRateOutputType{1};
  FrictionLawType frictionLawType{0};
  int backgroundType{0};
  bool isRfOutputOn{false};
  bool isDsOutputOn{false};
  bool isMagnitudeOutputOn{false};
  bool isEnergyRateOutputOn{false};
  bool isGpWiseOutput{false};
  bool isThermalPressureOn{false};
  int energyRatePrintTimeInterval{1};
  bool isInstaHealingOn{false};
  real t_0{0.0};
  real rs_f0{0.0};
  real rs_a{0.0};
  real rs_b{0.0};
  real rs_sr0{0.0};
  real mu_w{0.0};
  real alpha_th{0.0};
  real rho_c{0.0};
  real tP_lambda{0.0};
  real iniTemp{0.0};
  real iniPressure{0.0};
  real v_star{0.0}; // Prakash-Clifton regularization parameter
  real prakash_length{0.0};

  void setAllInputParam(const YAML::Node& Params) {
    using namespace initializers;

    const YAML::Node& DrParams = Params["dynamicrupture"];
    outputPointType = getParamIfExists(DrParams, "outputpointtype", 3);
    slipRateOutputType = getParamIfExists(DrParams, "sliprateoutputtype", 1);
    frictionLawType = static_cast<FrictionLawType>(getParamIfExists(DrParams, "fl", 0));
    backgroundType = getParamIfExists(DrParams, "backgroundtype", 0);
    isRfOutputOn = getParamIfExists(DrParams, "rf_output_on", false);
    isDsOutputOn = getParamIfExists(DrParams, "ds_output_on", false);
    isMagnitudeOutputOn = getParamIfExists(DrParams, "magnitude_output_on", false);
    isEnergyRateOutputOn = getParamIfExists(DrParams, "energy_rate_output_on", false);
    isGpWiseOutput = getParamIfExists(DrParams, "gpwise", false);
    isThermalPressureOn = getParamIfExists(DrParams, "thermalpress", false);
    backgroundType = getParamIfExists(DrParams, "energy_rate_printtimeinterval", 1);
    isInstaHealingOn = getParamIfExists(DrParams, "inst_healing", false);
    t_0 = getParamIfExists(DrParams, "t_0", 0.0);
    rs_f0 = getParamIfExists(DrParams, "rs_f0", 0.0);
    rs_a = getParamIfExists(DrParams, "rs_a", 0.0);
    rs_b = getParamIfExists(DrParams, "rs_b", 0.0);
    rs_sr0 = getParamIfExists(DrParams, "rs_sr0", 0.0);
    mu_w = getParamIfExists(DrParams, "mu_w", 0.0);

    // if ThermalPress == true
    alpha_th = getParamIfExists(DrParams, "alpha_th", 0.0);
    rho_c = getParamIfExists(DrParams, "rho_c", 0.0);
    tP_lambda = getParamIfExists(DrParams, "tp_lambda", 0.0);
    iniTemp = getParamIfExists(DrParams, "initemp", 0.0);
    iniPressure = getParamIfExists(DrParams, "inipressure", 0.0);

    // Prakash-Clifton regularization parameter
    v_star = getParamIfExists(DrParams, "v_star", 0.0);
    prakash_length = getParamIfExists(DrParams, "L", 0.0);
  }
};

#endif // SEISSOL_PARAMETERS_H
