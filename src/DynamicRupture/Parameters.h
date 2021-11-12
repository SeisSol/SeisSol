#ifndef SEISSOL_PARAMETERS_H
#define SEISSOL_PARAMETERS_H

#include <yaml-cpp/yaml.h>

#include "Initializer/InputAux.hpp"
#include "DynamicRupture/Typedefs.hpp"
#include "Kernels/precision.hpp"
#include "Typedefs.hpp"

#include <Eigen/Dense>

namespace seissol::dr {
struct DRParameters;
inline DRParameters readParametersFromYaml(YAML::Node& params);
} // namespace seissol::dr

/*
 * Saves all dynamic rupture parameter read from parameter.par file
 * if values are not defined they are set to an initial value (mostly 0)
 */
struct seissol::dr::DRParameters {
  static constexpr unsigned int TP_grid_nz = 60;
  int outputPointType{3};
  Eigen::Vector3d referencePoint;
  int slipRateOutputType{1};
  FrictionLawType frictionLawType{0};
  int backgroundType{0};
  bool isRfOutputOn{false};
  bool isDsOutputOn{false};
  bool isMagnitudeOutputOn{false};
  bool isEnergyRateOutputOn{false};
  bool isGpWiseInitialization{true};
  bool isThermalPressureOn{false};
  int energyRatePrintTimeInterval{1};
  bool isInstaHealingOn{false};
  real t_0{0.0};
  real rs_f0{0.0};
  real rs_a{0.0};
  real rs_b{0.0};
  real rs_sr0{0.0};
  real rs_initialSlipRate1{0.0};
  real rs_initialSlipRate2{0.0};
  real mu_w{0.0};
  real alpha_th{0.0};
  real rho_c{0.0};
  real tP_lambda{0.0};
  real iniTemp{0.0};
  real iniPressure{0.0};
  real v_star{0.0}; // Prakash-Clifton regularization parameter
  real prakash_length{0.0};
  std::string faultFileName{""};
};

inline seissol::dr::DRParameters seissol::dr::readParametersFromYaml(YAML::Node& params) {
  DRParameters drParameters;
  const YAML::Node& yamlParams = params["dynamicrupture"];
  
  if (params["dynamicrupture"]) {
    double xref = 0.0;
    initializers::updateIfExists(yamlParams, "xref", xref);
    double yref = 0.0;
    initializers::updateIfExists(yamlParams, "yref", yref);
    double zref = 0.0;
    initializers::updateIfExists(yamlParams, "zref", zref);
    drParameters.referencePoint = {xref, yref, zref};

    initializers::updateIfExists(yamlParams, "outputpointtype", drParameters.outputPointType);
    initializers::updateIfExists(yamlParams, "sliprateoutputtype", drParameters.slipRateOutputType);
    initializers::updateIfExists(yamlParams, "fl", drParameters.frictionLawType);
    initializers::updateIfExists(yamlParams, "backgroundtype", drParameters.backgroundType);
    initializers::updateIfExists(yamlParams, "rf_output_on", drParameters.isRfOutputOn);
    initializers::updateIfExists(yamlParams, "ds_output_on", drParameters.isDsOutputOn);
    initializers::updateIfExists(yamlParams, "magnitude_output_on", drParameters.isMagnitudeOutputOn);
    initializers::updateIfExists(
        yamlParams, "energy_rate_output_on", drParameters.isEnergyRateOutputOn);
    initializers::updateIfExists(yamlParams, "gpwise", drParameters.isGpWiseInitialization);
    initializers::updateIfExists(yamlParams, "thermalpress", drParameters.isThermalPressureOn);
    initializers::updateIfExists(
        yamlParams, "energy_rate_printtimeinterval", drParameters.backgroundType);
    initializers::updateIfExists(yamlParams, "inst_healing", drParameters.isInstaHealingOn);
    initializers::updateIfExists(yamlParams, "t_0", drParameters.t_0);
    initializers::updateIfExists(yamlParams, "rs_f0", drParameters.rs_f0);
    initializers::updateIfExists(yamlParams, "rs_a", drParameters.rs_a);
    initializers::updateIfExists(yamlParams, "rs_b", drParameters.rs_b);
    initializers::updateIfExists(yamlParams, "rs_sr0", drParameters.rs_sr0);
    initializers::updateIfExists(yamlParams, "rs_inisliprate1", drParameters.rs_initialSlipRate1);
    initializers::updateIfExists(yamlParams, "rs_inisliprate2", drParameters.rs_initialSlipRate2);
    initializers::updateIfExists(yamlParams, "mu_w", drParameters.mu_w);

    // Thermal Pressurisation parameters
    initializers::updateIfExists(yamlParams, "alpha_th", drParameters.alpha_th);
    initializers::updateIfExists(yamlParams, "rho_c", drParameters.rho_c);
    initializers::updateIfExists(yamlParams, "tp_lambda", drParameters.tP_lambda);
    initializers::updateIfExists(yamlParams, "initemp", drParameters.iniTemp);
    initializers::updateIfExists(yamlParams, "inipressure", drParameters.iniPressure);

    // Prakash-Clifton regularization parameters
    initializers::updateIfExists(yamlParams, "v_star", drParameters.v_star);
    initializers::updateIfExists(yamlParams, "L", drParameters.prakash_length);

    // filename of the yaml file describing the fault parameters
    initializers::updateIfExists(yamlParams, "modelfilename", drParameters.faultFileName);
  }

  return drParameters;
}
#endif // SEISSOL_PARAMETERS_H
