#ifndef SEISSOL_PARAMETERS_H
#define SEISSOL_PARAMETERS_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Typedefs.hpp"
#include "Initializer/InputAux.hpp"
#include "Kernels/precision.hpp"
#include "Typedefs.hpp"

#include <Eigen/Dense>

namespace seissol::dr {
/**
 * Saves all dynamic rupture parameter read from parameter.par file
 * if values are not defined they are set to an initial value (mostly 0)
 */
struct DRParameters {
  bool isDynamicRuptureEnabled{true};
  int outputPointType{3};
  Eigen::Vector3d referencePoint;
  int slipRateOutputType{1};
  FrictionLawType frictionLawType{0};
  int backgroundType{0};
  bool isRfOutputOn{false};
  bool isDsOutputOn{false};
  bool isThermalPressureOn{false};
  real t0{0.0};
  real rsF0{0.0};
  real rsA{0.0};
  real rsB{0.0};
  real rsSr0{0.0};
  real rsInitialSlipRate1{0.0};
  real rsInitialSlipRate2{0.0};
  real muW{0.0};
  real thermalDiffusivity{0.0};
  real heatCapacity{0.0};
  real undrainedTPResponse{0.0};
  real initialTemperature{0.0};
  real initialPressure{0.0};
  real vStar{0.0}; // Prakash-Clifton regularization parameter
  real prakashLength{0.0};
  std::string faultFileName{""};
};

inline DRParameters readParametersFromYaml(YAML::Node& params) {
  DRParameters drParameters;
  const YAML::Node& yamlDrParams = params["dynamicrupture"];

  if (params["dynamicrupture"]) {
    double xref = 0.0;
    initializers::updateIfExists(yamlDrParams, "xref", xref);
    double yref = 0.0;
    initializers::updateIfExists(yamlDrParams, "yref", yref);
    double zref = 0.0;
    initializers::updateIfExists(yamlDrParams, "zref", zref);
    drParameters.referencePoint = {xref, yref, zref};

    initializers::updateIfExists(yamlDrParams, "outputpointtype", drParameters.outputPointType);
    initializers::updateIfExists(
        yamlDrParams, "sliprateoutputtype", drParameters.slipRateOutputType);
    initializers::updateIfExists(yamlDrParams, "fl", drParameters.frictionLawType);
    initializers::updateIfExists(yamlDrParams, "backgroundtype", drParameters.backgroundType);
    initializers::updateIfExists(yamlDrParams, "thermalpress", drParameters.isThermalPressureOn);
    initializers::updateIfExists(yamlDrParams, "t_0", drParameters.t0);
    initializers::updateIfExists(yamlDrParams, "rs_f0", drParameters.rsF0);
    initializers::updateIfExists(yamlDrParams, "rs_a", drParameters.rsA);
    initializers::updateIfExists(yamlDrParams, "rs_b", drParameters.rsB);
    initializers::updateIfExists(yamlDrParams, "rs_sr0", drParameters.rsSr0);
    initializers::updateIfExists(yamlDrParams, "rs_inisliprate1", drParameters.rsInitialSlipRate1);
    initializers::updateIfExists(yamlDrParams, "rs_inisliprate2", drParameters.rsInitialSlipRate2);
    initializers::updateIfExists(yamlDrParams, "rs_muw", drParameters.muW);

    // Thermal Pressurization parameters
    initializers::updateIfExists(
        yamlDrParams, "tp_thermaldiffusivity", drParameters.thermalDiffusivity);
    initializers::updateIfExists(yamlDrParams, "tp_heatcapacity", drParameters.heatCapacity);
    initializers::updateIfExists(
        yamlDrParams, "tp_undrainedtpresponse", drParameters.undrainedTPResponse);
    initializers::updateIfExists(yamlDrParams, "tp_initemp", drParameters.initialTemperature);
    initializers::updateIfExists(yamlDrParams, "tp_inipressure", drParameters.initialPressure);

    // Prakash-Clifton regularization parameters
    initializers::updateIfExists(yamlDrParams, "pc_vstar", drParameters.vStar);
    initializers::updateIfExists(yamlDrParams, "pc_prakashlength", drParameters.prakashLength);

    // filename of the yaml file describing the fault parameters
    initializers::updateIfExists(yamlDrParams, "modelfilename", drParameters.faultFileName);
  }

  const YAML::Node& yamlElementwiseParams = params["elementwise"];
  if (params["elementwise"]) {
    // check whether we need rupture time and dynamic stress time outputs
    std::array<bool, 12> mask;
    initializers::convertStringToMask(yamlElementwiseParams["outputmask"].as<std::string>(), mask);
    drParameters.isRfOutputOn = drParameters.isRfOutputOn || mask[9];
    drParameters.isDsOutputOn = drParameters.isDsOutputOn || mask[10];
  }
  const YAML::Node& yamlPickpointParams = params["pickpoint"];
  if (params["pickpoint"]) {
    std::array<bool, 12> mask;
    initializers::convertStringToMask(yamlPickpointParams["outputmask"].as<std::string>(), mask);
    drParameters.isRfOutputOn = drParameters.isRfOutputOn || mask[9];
    drParameters.isDsOutputOn = drParameters.isDsOutputOn || mask[10];
  }
  // if there is no filename given for the fault, assume that we do not use dynamic rupture
  if (drParameters.faultFileName == "") {
    drParameters.isDynamicRuptureEnabled = false;
  }

  return drParameters;
}
} // namespace seissol::dr
#endif // SEISSOL_PARAMETERS_H
