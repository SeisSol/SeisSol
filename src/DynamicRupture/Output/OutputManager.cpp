#include "Common/filesystem.h"
#include "DynamicRupture/Output/OutputManager.hpp"
#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"
#include "SeisSol.h"
#include <Initializer/Parameters/OutputParameters.h>
#include <Initializer/Parameters/SeisSolParameters.h>
#include <fstream>
#include <type_traits>
#include <unordered_map>

struct NativeFormat {};
struct WideFormat {};
template <typename T, typename U = NativeFormat>
struct FormattedBuildinType {
  T value;
};

template <typename T, typename U = NativeFormat>
auto makeFormatted(T value) {
  return FormattedBuildinType<T, U>{value};
}

template <typename T, typename U = NativeFormat>
std::ostream& operator<<(std::ostream& stream, FormattedBuildinType<T, U> obj) {
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
std::string buildFileName(std::string namePrefix,
                          std::string nameSuffix,
                          std::string fileExtension = std::string()) {
  std::stringstream fileName;
  fileName << namePrefix << '-' << nameSuffix;
  if (fileExtension.empty()) {
    return fileName.str();
  } else {
    fileName << '.' << fileExtension;
    return fileName.str();
  }
}

std::string buildMPIFileName(std::string namePrefix,
                             std::string nameSuffix,
                             std::string fileExtension = std::string()) {
#ifdef PARALLEL
  std::stringstream suffix;
  suffix << nameSuffix << '-' << makeFormatted<int, WideFormat>(MPI::mpi.rank());
  return buildFileName(namePrefix, suffix.str(), fileExtension);
#else
  return buildFileName(namePrefix, nameSuffix, fileExtension);
#endif
}

std::string buildIndexedMPIFileName(std::string namePrefix,
                                    int index,
                                    std::string nameSuffix,
                                    std::string fileExtension = std::string()) {
  std::stringstream suffix;
#ifdef PARALLEL
  suffix << nameSuffix << '-' << makeFormatted<int, WideFormat>(index) << '-'
         << makeFormatted<int, WideFormat>(MPI::mpi.rank());
#else
  suffix << nameSuffix << '-' << makeFormatted<int, WideFormat>(index);
#endif
  return buildFileName(namePrefix, suffix.str(), fileExtension);
}

OutputManager::OutputManager(std::unique_ptr<ReceiverOutput> concreteImpl,
                             seissol::SeisSol& seissolInstance)
    : seissolInstance(seissolInstance), ewOutputData(std::make_shared<ReceiverOutputData>()),
      ppOutputData(std::make_shared<ReceiverOutputData>()), impl(std::move(concreteImpl)) {
  backupTimeStamp = utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(0L));
}

OutputManager::~OutputManager() { flushPickpointDataToFile(); }

void OutputManager::setInputParam(seissol::geometry::MeshReader& userMesher) {
  using namespace initializer;
  meshReader = &userMesher;

  impl->setMeshReader(&userMesher);

  const auto& seissolParameters = seissolInstance.getSeisSolParameters();
  const bool bothEnabled = seissolParameters.drParameters.outputPointType ==
                           seissol::initializer::parameters::OutputType::AtPickpointAndElementwise;
  const bool pointEnabled = seissolParameters.drParameters.outputPointType ==
                                seissol::initializer::parameters::OutputType::AtPickpoint ||
                            bothEnabled;
  const bool elementwiseEnabled = seissolParameters.drParameters.outputPointType ==
                                      seissol::initializer::parameters::OutputType::Elementwise ||
                                  bothEnabled;
  const int rank = seissol::MPI::mpi.rank();
  if (pointEnabled) {
    logInfo(rank) << "Enabling on-fault receiver output";
    ppOutputBuilder = std::make_unique<PickPointBuilder>();
    ppOutputBuilder->setMeshReader(&userMesher);
    ppOutputBuilder->setParams(seissolParameters.output.pickpointParameters);
  }
  if (elementwiseEnabled) {
    logInfo(rank) << "Enabling 2D fault output";
    ewOutputBuilder = std::make_unique<ElementWiseBuilder>();
    ewOutputBuilder->setMeshReader(&userMesher);
    ewOutputBuilder->setParams(seissolParameters.output.elementwiseParameters);
  }
  if (!elementwiseEnabled && !pointEnabled) {
    logInfo(rank) << "No dynamic rupture output enabled";
  }
}

void OutputManager::setLtsData(seissol::initializer::LTSTree* userWpTree,
                               seissol::initializer::LTS* userWpDescr,
                               seissol::initializer::Lut* userWpLut,
                               seissol::initializer::LTSTree* userDrTree,
                               seissol::initializer::DynamicRupture* userDrDescr) {
  wpDescr = userWpDescr;
  wpTree = userWpTree;
  wpLut = userWpLut;
  drTree = userDrTree;
  drDescr = userDrDescr;
  impl->setLtsData(wpTree, wpDescr, wpLut, drTree, drDescr);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();
  const bool bothEnabled = seissolParameters.drParameters.outputPointType ==
                           seissol::initializer::parameters::OutputType::AtPickpointAndElementwise;
  const bool pointEnabled = seissolParameters.drParameters.outputPointType ==
                                seissol::initializer::parameters::OutputType::AtPickpoint ||
                            bothEnabled;
  const bool elementwiseEnabled = seissolParameters.drParameters.outputPointType ==
                                      seissol::initializer::parameters::OutputType::Elementwise ||
                                  bothEnabled;
  if (pointEnabled) {
    ppOutputBuilder->setLtsData(userWpTree, userWpDescr, userWpLut);
  }
  if (elementwiseEnabled) {
    ewOutputBuilder->setLtsData(userWpTree, userWpDescr, userWpLut);
  }
}

void OutputManager::initElementwiseOutput() {
  ewOutputBuilder->build(ewOutputData);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  const auto& receiverPoints = ewOutputData->receiverPoints;
  const auto cellConnectivity = getCellConnectivity(receiverPoints);
  const auto faultTags = getFaultTags(receiverPoints);
  const auto vertices = getAllVertices(receiverPoints);
  constexpr auto maxNumVars = std::tuple_size<DrVarsT>::value;
  const auto outputMask = seissolParameters.output.elementwiseParameters.outputMask;
  const auto intMask = convertMaskFromBoolToInt<maxNumVars>(outputMask);

  const double printTime = seissolParameters.output.elementwiseParameters.printTimeIntervalSec;
  const auto backendType = seissolParameters.output.xdmfWriterBackend;

  std::vector<real*> dataPointers;
  auto recordPointers = [&dataPointers](auto& var, int) {
    if (var.isActive) {
      for (int dim = 0; dim < var.dim(); ++dim)
        dataPointers.push_back(var.data[dim]);
    }
  };
  misc::forEach(ewOutputData->vars, recordPointers);

  seissolInstance.faultWriter().init(cellConnectivity.data(),
                                     vertices.data(),
                                     faultTags.data(),
                                     static_cast<unsigned int>(receiverPoints.size()),
                                     static_cast<unsigned int>(3 * receiverPoints.size()),
                                     &intMask[0],
                                     const_cast<const real**>(dataPointers.data()),
                                     seissolParameters.output.prefix.data(),
                                     printTime,
                                     backendType,
                                     backupTimeStamp);

  seissolInstance.faultWriter().setupCallbackObject(this);
}

void OutputManager::initPickpointOutput() {
  ppOutputBuilder->build(ppOutputData);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  std::stringstream baseHeader;
  baseHeader << "VARIABLES = \"Time\"";
  size_t labelCounter = 0;
  auto collectVariableNames = [&baseHeader, &labelCounter](auto& var, int) {
    if (var.isActive) {
      for (int dim = 0; dim < var.dim(); ++dim) {
        baseHeader << " ,\"" << writer::FaultWriterExecutor::getLabelName(labelCounter) << '\"';
        ++labelCounter;
      }
    } else {
      labelCounter += var.dim();
    }
  };
  misc::forEach(ppOutputData->vars, collectVariableNames);

  auto& outputData = ppOutputData;
  for (const auto& receiver : outputData->receiverPoints) {
    const size_t globalIndex = receiver.globalReceiverIndex + 1;

    auto fileName =
        buildIndexedMPIFileName(seissolParameters.output.prefix, globalIndex, "faultreceiver");
    seissol::generateBackupFileIfNecessary(fileName, "dat", {backupTimeStamp});
    fileName += ".dat";

    if (!seissol::filesystem::exists(fileName)) {
      std::ofstream file(fileName, std::ios_base::out);
      if (file.is_open()) {
        std::stringstream title;
        title << "TITLE = \"Temporal Signal for fault receiver number " << globalIndex << "\"";

        file << title.str() << '\n';
        file << baseHeader.str() << '\n';

        const auto& point = const_cast<ExtVrtxCoords&>(receiver.global);
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
  impl->allocateMemory({ppOutputData, ewOutputData});
}

void OutputManager::initFaceToLtsMap() {
  if (drTree) {
    const size_t readerFaultSize = meshReader->getFault().size();
    const size_t ltsFaultSize = drTree->getNumberOfCells(Ghost);

    faceToLtsMap.resize(std::max(readerFaultSize, ltsFaultSize));
    for (auto it = drTree->beginLeaf(seissol::initializer::LayerMask(Ghost));
         it != drTree->endLeaf();
         ++it) {

      DRFaceInformation* faceInformation = it->var(drDescr->faceInformation);
      for (size_t ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        faceToLtsMap[faceInformation[ltsFace].meshFace] = std::make_pair(&(*it), ltsFace);
      }
    }
  }
  impl->setFaceToLtsMap(&faceToLtsMap);
}

bool OutputManager::isAtPickpoint(double time, double dt) {
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();
  const bool isFirstStep = iterationStep == 0;
  const double abortTime = seissolParameters.timeStepping.endTime;
  const bool isCloseToTimeOut = (abortTime - time) < (dt * timeMargin);

  const int printTimeInterval = seissolParameters.output.pickpointParameters.printTimeInterval;
  const bool isOutputIteration = iterationStep % printTimeInterval == 0;

  return (isFirstStep || isOutputIteration || isCloseToTimeOut);
}

void OutputManager::writePickpointOutput(double time, double dt) {
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();
  if (this->ppOutputBuilder) {
    if (this->isAtPickpoint(time, dt)) {

      const auto& outputData = ppOutputData;
      impl->calcFaultOutput(seissol::initializer::parameters::OutputType::AtPickpoint,
                            seissolParameters.drParameters.slipRateOutputType,
                            ppOutputData,
                            time);

      const bool isMaxCacheLevel =
          outputData->currentCacheLevel >=
          static_cast<size_t>(seissolParameters.output.pickpointParameters.maxPickStore);
      const bool isCloseToEnd = (seissolParameters.timeStepping.endTime - time) < dt * timeMargin;

      if (isMaxCacheLevel || isCloseToEnd) {
        this->flushPickpointDataToFile();
      }
    }
    ++iterationStep;
  }
}

void OutputManager::flushPickpointDataToFile() {
  auto& outputData = ppOutputData;
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  for (size_t pointId = 0; pointId < outputData->receiverPoints.size(); ++pointId) {
    std::stringstream data;
    for (size_t level = 0; level < outputData->currentCacheLevel; ++level) {
      data << makeFormatted(outputData->cachedTime[level]) << '\t';
      auto recordResults = [pointId, level, &data](auto& var, int) {
        if (var.isActive) {
          for (int dim = 0; dim < var.dim(); ++dim) {
            data << makeFormatted(var(dim, level, pointId)) << '\t';
          }
        }
      };
      misc::forEach(outputData->vars, recordResults);
      data << '\n';
    }

    const auto globalIndex = outputData->receiverPoints[pointId].globalReceiverIndex + 1;
    const auto fileName = buildIndexedMPIFileName(
        seissolParameters.output.prefix, globalIndex, "faultreceiver", "dat");

    std::ofstream file(fileName, std::ios_base::app);
    if (file.is_open()) {
      file << data.str();
    } else {
      logError() << "cannot open " << fileName;
    }
    file.close();
  }
  outputData->currentCacheLevel = 0;
}

void OutputManager::updateElementwiseOutput() {
  if (this->ewOutputBuilder) {
    const auto& seissolParameters = seissolInstance.getSeisSolParameters();
    impl->calcFaultOutput(seissol::initializer::parameters::OutputType::Elementwise,
                          seissolParameters.drParameters.slipRateOutputType,
                          ewOutputData);
  }
}
} // namespace seissol::dr::output
