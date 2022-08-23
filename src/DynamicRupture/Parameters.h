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

inline std::unique_ptr<DRParameters> readParametersFromYaml(std::shared_ptr<YAML::Node>& params) {
  std::unique_ptr<DRParameters> drParameters = std::make_unique<DRParameters>();
  const YAML::Node& yamlDrParams = (*params)["dynamicrupture"];

  using namespace seissol::initializers;

  if ((*params)["dynamicrupture"]) {
    double xref = getWithDefault(yamlDrParams, "xref", 0.0);
    double yref = getWithDefault(yamlDrParams, "yref", 0.0);
    double zref = getWithDefault(yamlDrParams, "zref", 0.0);
    drParameters->referencePoint = {xref, yref, zref};

    drParameters->outputPointType = getWithDefault(yamlDrParams, "outputpointtype", 3);
    drParameters->slipRateOutputType = getWithDefault(yamlDrParams, "sliprateoutputtype", 1);
    drParameters->frictionLawType =
        static_cast<FrictionLawType>(getWithDefault(yamlDrParams, "fl", 0));
    drParameters->backgroundType = getWithDefault(yamlDrParams, "backgroundtype", 0);
    drParameters->isThermalPressureOn = getWithDefault(yamlDrParams, "thermalpress", false);
    drParameters->t0 = getWithDefault(yamlDrParams, "t_0", 0.0);
    drParameters->rsF0 = getWithDefault(yamlDrParams, "rs_f0", 0.0);
    drParameters->rsA = getWithDefault(yamlDrParams, "rs_a", 0.0);
    drParameters->rsB = getWithDefault(yamlDrParams, "rs_b", 0.0);
    drParameters->rsSr0 = getWithDefault(yamlDrParams, "rs_sr0", 0.0);
    drParameters->rsInitialSlipRate1 = getWithDefault(yamlDrParams, "rs_inisliprate1", 0.0);
    drParameters->rsInitialSlipRate2 = getWithDefault(yamlDrParams, "rs_inisliprate2", 0.0);
    drParameters->muW = getWithDefault(yamlDrParams, "mu_w", 0.0);

    // Thermal Pressurization parameters
    drParameters->thermalDiffusivity = getWithDefault(yamlDrParams, "alpha_th", 0.0);
    drParameters->heatCapacity = getWithDefault(yamlDrParams, "rho_c", 0.0);
    drParameters->undrainedTPResponse = getWithDefault(yamlDrParams, "tp_lambda", 0.0);
    drParameters->initialTemperature = getWithDefault(yamlDrParams, "initemp", 0.0);
    drParameters->initialPressure = getWithDefault(yamlDrParams, "inipressure", 0.0);

    // Prakash-Clifton regularization parameters
    drParameters->vStar = getWithDefault(yamlDrParams, "vstar", 0.0);
    drParameters->prakashLength = getWithDefault(yamlDrParams, "prakashlength", 0.0);

    // filename of the yaml file describing the fault parameters
    drParameters->faultFileName = getWithDefault(yamlDrParams, "modelfilename", std::string(""));
  }

  const YAML::Node& yamlElementwiseParams = (*params)["elementwise"];
  if ((*params)["elementwise"]) {
    // check whether we need rupture time and dynamic stress time outputs
    std::array<bool, 12> mask;
    initializers::convertStringToMask(yamlElementwiseParams["outputmask"].as<std::string>(), mask);
    drParameters->isRfOutputOn = drParameters->isRfOutputOn || mask[9];
    drParameters->isDsOutputOn = drParameters->isDsOutputOn || mask[10];
  }
  const YAML::Node& yamlPickpointParams = (*params)["pickpoint"];
  if ((*params)["pickpoint"]) {
    std::array<bool, 12> mask;
    initializers::convertStringToMask(yamlPickpointParams["outputmask"].as<std::string>(), mask);
    drParameters->isRfOutputOn = drParameters->isRfOutputOn || mask[9];
    drParameters->isDsOutputOn = drParameters->isDsOutputOn || mask[10];
  }
  // if there is no filename given for the fault, assume that we do not use dynamic rupture
  if (drParameters->faultFileName == "") {
    drParameters->isDynamicRuptureEnabled = false;
  }

  return drParameters;
}
} // namespace seissol::dr
#endif // SEISSOL_PARAMETERS_H
