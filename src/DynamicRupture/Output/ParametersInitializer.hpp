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
    const auto outputPointID =
        getWithDefault(drSettings, "outputpointtype", static_cast<int>(OutputType::None));
    params.outputPointType = static_cast<OutputType>(outputPointID);

    auto slipRateOutputType = getWithDefault(drSettings, "sliprateoutputtype", 1);

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

    params.isThermalPressurizationOn = getWithDefault(drSettings, "thermalpress", false);

    const YAML::Node& outputParams = data["output"];
    params.faultOutputFlag = getWithDefault(outputParams, "faultoutputflag", false);
    params.outputFilePrefix = getWithDefault(outputParams, "outputfile", std::string("data"));
    params.checkPointBackend =
        getWithDefault(outputParams, "checkpointbackend", std::string("none"));

    const YAML::Node& abortCriteriaParams = data["abortcriteria"];
    params.endTime = getWithDefault(abortCriteriaParams, "endtime", 15.0);
    params.maxIteration =
        getWithDefault(abortCriteriaParams, "maxiteration", static_cast<size_t>(1000000000));

#ifdef USE_HDF
    params.xdmfWriterBackend =
        getWithDefault(outputParams, "xdmfwriterbackend", std::string("hdf5"));
#else
    params.xdmfWriterBackend =
        getWithDefault(outputParams, "xdmfwriterbackend", std::string("posix"));
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
    ppParams.printTimeInterval = getWithDefault(ppData, "printtimeinterval", 1);
    ppParams.ppFileName = getWithDefault(ppData, "ppfilename", std::string(""));

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
    ewParams.printTimeIntervalSec = getWithDefault(ewData, "printtimeinterval_sec", 1.0);
    ewParams.refinement = getWithDefault(ewData, "refinement", 2);

    auto refinementStrategy = getWithDefault(ewData, "refinement_strategy", 2);
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
