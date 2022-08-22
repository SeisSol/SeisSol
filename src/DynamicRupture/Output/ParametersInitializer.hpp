#ifndef SEISSOL_DR_OUTPUT_PARAMETERS_INITIALIZER_HPP
#define SEISSOL_DR_OUTPUT_PARAMETERS_INITIALIZER_HPP

#include "DataTypes.hpp"
#include "FaultRefiner/FaultRefiners.hpp"
#include "Initializer/InputAux.hpp"
#include <utils/logger.h>
#include <yaml-cpp/yaml.h>

namespace seissol::dr::output {
class ParametersInitializer {
  public:
  explicit ParametersInitializer(const YAML::Node& userData) : data(userData) {}

  GeneralParams getDrGeneralParams() {
    using namespace seissol::initializers;
    GeneralParams params{};

    if (!data["dynamicrupture"]) {
      logError() << "dynamic rupture params. is not provided in the namelist";
    }

    const YAML::Node& drSettings = data["dynamicrupture"];
    auto outputPointID = static_cast<int>(OutputType::None);
    updateIfExists(drSettings, "outputpointtype", outputPointID);
    params.outputPointType = static_cast<OutputType>(outputPointID);

    int slipRateOutputType{1};
    updateIfExists(drSettings, "sliprateoutputtype", slipRateOutputType);

    switch (slipRateOutputType) {
    case 0: {
      params.slipRateOutputType = SlipRateOutputType::VelocityDifference;
      break;
    }
    case 1: {
      params.slipRateOutputType = SlipRateOutputType::TractionsAndFailure;
      break;
    }
    default: {
      logError() << "given slip rate output type (" << slipRateOutputType << ") is not supported";
    }
    }

    updateIfExists(drSettings, "thermalpress", params.isThermalPressurizationOn);

    const YAML::Node& outputParams = data["output"];
    updateIfExists(outputParams, "faultoutputflag", params.faultOutputFlag);
    updateIfExists(outputParams, "outputfile", params.outputFilePrefix);
    updateIfExists(outputParams, "checkpointbackend", params.checkPointBackend);

    const YAML::Node& abortCriteriaParams = data["abortcriteria"];
    updateIfExists(abortCriteriaParams, "endtime", params.endTime);
    updateIfExists(abortCriteriaParams, "maxiteration", params.maxIteration);

#ifdef USE_HDF
    updateIfExists(outputParams, "xdmfwriterbackend", params.xdmfWriterBackend);
#else
    params.xdmfWriterBackend = std::string("posix");
    updateIfExists(outputParams, "xdmfwriterbackend", params.xdmfWriterBackend);
#endif

    return params;
  }

  PickpointParams getPickPointParams() {
    using namespace seissol::initializers;
    PickpointParams ppParams{};

    if (!data["pickpoint"]) {
      logError() << "pickpoint output parameters for dynamic rupture is not provided";
    }

    const YAML::Node& ppData = data["pickpoint"];
    updateIfExists(ppData, "printtimeinterval", ppParams.printTimeInterval);
    updateIfExists(ppData, "ppfilename", ppParams.ppFileName);

    if (ppData["outputmask"]) {
      convertStringToMask(ppData["outputmask"].as<std::string>(), ppParams.outputMask);
    }

    return ppParams;
  }

  ElementwiseFaultParams getElementwiseFaultParams() {
    using namespace seissol::initializers;

    ElementwiseFaultParams ewParams{};

    if (!data["elementwise"]) {
      logError() << "elementwise fault output parameters for dynamic rupture is not provided";
    }

    const YAML::Node& ewData = data["elementwise"];
    updateIfExists(ewData, "printtimeinterval_sec", ewParams.printTimeIntervalSec);
    updateIfExists(ewData, "refinement", ewParams.refinement);

    int refinementStrategy{2};
    updateIfExists(ewData, "refinement_strategy", refinementStrategy);
    ewParams.refinementStrategy = refiner::castToRefinerType(refinementStrategy);

    if (ewData["outputmask"]) {
      convertStringToMask(ewData["outputmask"].as<std::string>(), ewParams.outputMask);
    }

    return ewParams;
  }

  private:
  const YAML::Node& data;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_PARAMETERS_INITIALIZER_HPP
