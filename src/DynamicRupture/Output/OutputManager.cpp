#include "OutputManager.hpp"
#include "DynamicRupture/Output/Base.hpp"
#include "SeisSol.h"
#include "ResultWriter/common.hpp"
#include <filesystem>
#include <fstream>
#include <type_traits>
#include <unordered_map>

struct NativeFormat {};
struct WideFormat {};
template <typename T, typename U = NativeFormat>
struct FormattedBuildInType {
  T value;
};

template <typename T, typename U = NativeFormat>
auto makeFormatted(T value) {
  return FormattedBuildInType<T, U>{value};
}

template <typename T, typename U = NativeFormat>
std::ostream& operator<<(std::ostream& stream, FormattedBuildInType<T, U> obj) {
  if constexpr (std::is_floating_point_v<T>) {
    stream << std::setprecision(16) << std::scientific << obj.value;
  } else if constexpr (std::is_integral_v<T> && std::is_same_v<U, WideFormat>) {
    stream << std::setw(5) << std::setfill('0') << obj.value;
  } else {
    stream << obj.value;
  }
  return stream;
}

namespace seissol::dr::output {
OutputManager::~OutputManager() {
  auto deallocateVars = [](auto& var, int) { var.releaseData(); };
  misc::forEach(ppOutputData.vars, deallocateVars);
  misc::forEach(ewOutputData.vars, deallocateVars);
}

void OutputManager::setInputParam(const YAML::Node& inputData, MeshReader& userMesher) {
  using namespace initializers;

  ParametersInitializer reader(inputData);
  impl->setMeshReader(&userMesher);
  generalParams = reader.getDrGeneralParams();

  // adjust general output parameters
  generalParams.isRfTimeOn = generalParams.isRfOutputOn;
  if (generalParams.isDsOutputOn && !generalParams.isRfOutputOn) {
    generalParams.isRfOutputOn = true;
    generalParams.isRfTimeOn = true;
  }

  bool bothEnabled = generalParams.outputPointType == OutputType::AtPickpointAndElementwise;
  bool pointEnabled = generalParams.outputPointType == OutputType::AtPickpoint || bothEnabled;
  bool elementwiseEnabled = generalParams.outputPointType == OutputType::Elementwise || bothEnabled;
  if (pointEnabled) {
    logInfo() << "Enabling on-fault receiver output";
    ppOutputBuilder = std::make_unique<PickPointBuilder>();
    ppOutputBuilder->setMeshReader(&userMesher);
    pickpointParams = reader.getPickPointParams();
    ppOutputBuilder->setParams(pickpointParams);
  }
  if (elementwiseEnabled) {
    logInfo() << "Enabling 2D fault output";
    ewOutputBuilder = std::make_unique<ElementWiseBuilder>();
    ewOutputBuilder->setMeshReader(&userMesher);
    elementwiseParams = reader.getElementwiseFaultParams();
    ewOutputBuilder->setParams(elementwiseParams);
  }
  if (!elementwiseEnabled && !pointEnabled) {
    logInfo() << "No dynamic rupture output enabled";
  }
}

void OutputManager::setLtsData(seissol::initializers::LTSTree* wpTree,
                               seissol::initializers::LTS* wpDescr,
                               seissol::initializers::Lut* wpLut,
                               seissol::initializers::LTSTree* drTree,
                               seissol::initializers::DynamicRupture* drDescr) {
  impl->setLtsData(wpTree, wpDescr, wpLut, drTree, drDescr);
}

std::string OutputManager::constructPickpointReceiverFileName(const int receiverGlobalIndex) const {
  std::stringstream fileName;
  fileName << generalParams.outputFilePrefix << "-new-faultreceiver-"
           << makeFormatted<int, WideFormat>(receiverGlobalIndex);
#ifdef PARALLEL
  fileName << '-' << makeFormatted<int, WideFormat>(MPI::mpi.rank());
#endif
  fileName << ".dat";
  return fileName.str();
}

void OutputManager::initElementwiseOutput() {
  ewOutputBuilder->build(&ewOutputData);

  const auto& receiverPoints = ewOutputData.receiverPoints;
  auto cellConnectivity = getCellConnectivity(receiverPoints);
  auto vertices = getAllVertices(receiverPoints);
  constexpr auto maxNumVars = std::tuple_size<DrVarsT>::value;
  auto intMask = convertMaskFromBoolToInt<maxNumVars>(this->elementwiseParams.outputMask);

  std::string newFaultFilePrefix = generalParams.outputFilePrefix + std::string("-new");
  double printTime = this->elementwiseParams.printTimeIntervalSec;
  auto backendType = seissol::writer::backendType(generalParams.xdmfWriterBackend.c_str());

  std::vector<real*> dataPointers;
  auto recordPointers = [&dataPointers](auto& var, int) {
    if (var.isActive) {
      for (int dim = 0; dim < var.dim(); ++dim)
        dataPointers.push_back(var.data[dim]);
    }
  };
  misc::forEach(ewOutputData.vars, recordPointers);

  seissol::SeisSol::main.secondFaultWriter().init(
      cellConnectivity.data(),
      vertices.data(),
      static_cast<unsigned int>(receiverPoints.size()),
      static_cast<unsigned int>(3 * receiverPoints.size()),
      &intMask[0],
      const_cast<const real**>(dataPointers.data()),
      newFaultFilePrefix.data(),
      printTime,
      backendType);

  seissol::SeisSol::main.secondFaultWriter().setupCallbackObject(this);
}

void OutputManager::initPickpointOutput() {
  ppOutputBuilder->build(&ppOutputData);

  std::stringstream baseHeader;
  baseHeader << "VARIABLES = \"Time\"";
  size_t labelCounter = 0;
  auto collectVariableNames = [&baseHeader, &labelCounter](auto& var, int i) {
    if (var.isActive) {
      for (int dim = 0; dim < var.dim(); ++dim) {
        baseHeader << " ,\"" << writer::FaultWriterExecutor::getLabelName(labelCounter) << '\"';
        ++labelCounter;
      }
    } else {
      labelCounter += var.dim();
    }
  };
  misc::forEach(ppOutputData.vars, collectVariableNames);

  auto& outputData = ppOutputData;
  for (const auto& receiver : outputData.receiverPoints) {
    const size_t globalIndex = receiver.globalReceiverIndex;
    auto fileName = constructPickpointReceiverFileName(globalIndex);
    if (!std::filesystem::exists(fileName)) {

      std::ofstream file(fileName, std::ios_base::out);
      if (file.is_open()) {
        std::stringstream title;
        title << "TITLE = \"Temporal Signal for fault receiver number " << globalIndex << "\"";

        file << title.str() << '\n';
        file << baseHeader.str() << '\n';

        auto& point = const_cast<ExtVrtxCoords&>(receiver.global);
        file << "# x1\t" << makeFormatted(point[0]) << '\n';
        file << "# x2\t" << makeFormatted(point[1]) << '\n';
        file << "# x3\t" << makeFormatted(point[2]) << '\n';

      } else {
        logError() << "cannot open " << fileName;
      }
      file.close();
    }
  }
}

void OutputManager::init() {
  if (ewOutputBuilder) {
    initElementwiseOutput();
  }
  if (ppOutputBuilder) {
    initPickpointOutput();
  }
}

void OutputManager::initFaceToLtsMap() { impl->initFaceToLtsMap(); }

bool OutputManager::isAtPickpoint(double time, double dt) {
  bool isFirstStep = iterationStep == 0;

  const double abortTime = std::min(generalParams.endTime, generalParams.maxIteration * dt);
  const bool isCloseToTimeOut = (abortTime - time) < (dt * timeMargin);

  const int printTimeInterval = this->pickpointParams.printTimeInterval;
  const bool isOutputIteration = iterationStep % printTimeInterval == 0;

  return ppOutputBuilder && (isFirstStep || isOutputIteration || isCloseToTimeOut);
}

void OutputManager::writePickpointOutput(double time, double dt) {
  if (this->ppOutputBuilder) {
    if (this->isAtPickpoint(time, dt)) {

      auto& outputData = ppOutputData;
      impl->calcFaultOutput(OutputType::AtPickpoint, ppOutputData, generalParams, time);

      const bool isMaxCacheLevel =
          outputData.currentCacheLevel >= static_cast<size_t>(this->pickpointParams.maxPickStore);
      const bool isCloseToEnd = (generalParams.endTime - time) < dt * timeMargin;
      if (isMaxCacheLevel || isCloseToEnd) {
        for (size_t pointId = 0; pointId < outputData.receiverPoints.size(); ++pointId) {

          std::stringstream data;
          for (size_t level = 0; level < outputData.currentCacheLevel; ++level) {
            data << makeFormatted(outputData.cachedTime[level]) << '\t';
            auto recordResults = [pointId, level, &data](auto& var, int) {
              if (var.isActive) {
                for (int dim = 0; dim < var.dim(); ++dim) {
                  data << makeFormatted(var(dim, level, pointId)) << '\t';
                }
              }
            };
            misc::forEach(outputData.vars, recordResults);
            data << '\n';
          }

          auto fileName = constructPickpointReceiverFileName(
              outputData.receiverPoints[pointId].globalReceiverIndex);
          std::ofstream file(fileName, std::ios_base::app);
          if (file.is_open()) {
            file << data.str();
          } else {
            logError() << "cannot open " << fileName;
          }
          file.close();
        }
        outputData.currentCacheLevel = 0;
      }
    }
    iterationStep += 1;
  }
}

void OutputManager::updateElementwiseOutput() {
  if (this->ewOutputBuilder) {
    impl->calcFaultOutput(OutputType::Elementwise, ewOutputData, generalParams);
  }
}

void OutputManager::tiePointers(seissol::initializers::Layer& layerData,
                                seissol::initializers::DynamicRupture* description,
                                seissol::Interoperability& eInteroperability) {
  impl->tiePointers(layerData, description, eInteroperability);
}
} // namespace seissol::dr::output