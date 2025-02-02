// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/Output/OutputManager.h"
#include "Common/Filesystem.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/Builders/ElementWiseBuilder.h"
#include "DynamicRupture/Output/Builders/PickPointBuilder.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/Geometry.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Writer/Writer.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Memory/Tree/Lut.h"
#include "ResultWriter/FaultWriterExecutor.h"
#include "SeisSol.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <fstream>
#include <init.h>
#include <iomanip>
#include <ios>
#include <kernel.h>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <utils/logger.h>
#include <utils/timeutils.h>
#include <vector>

namespace {

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

std::string buildFileName(const std::string& namePrefix,
                          const std::string& nameSuffix,
                          const std::string& fileExtension = std::string()) {
  std::stringstream fileName;
  fileName << namePrefix << '-' << nameSuffix;
  if (fileExtension.empty()) {
    return fileName.str();
  } else {
    fileName << '.' << fileExtension;
    return fileName.str();
  }
}

std::string buildIndexedMPIFileName(const std::string& namePrefix,
                                    int index,
                                    const std::string& nameSuffix,
                                    const std::string& fileExtension = std::string()) {
  std::stringstream suffix;
#ifdef PARALLEL
  suffix << nameSuffix << '-' << makeFormatted<int, WideFormat>(index) << '-'
         << makeFormatted<int, WideFormat>(seissol::MPI::mpi.rank());
#else
  suffix << nameSuffix << '-' << makeFormatted<int, WideFormat>(index);
#endif
  return buildFileName(namePrefix, suffix.str(), fileExtension);
}

} // namespace

namespace seissol::dr::output {

OutputManager::OutputManager(std::unique_ptr<ReceiverOutput> concreteImpl,
                             seissol::SeisSol& seissolInstance)
    : seissolInstance(seissolInstance), ewOutputData(std::make_shared<ReceiverOutputData>()),
      ppOutputData(std::make_shared<ReceiverOutputData>()), impl(std::move(concreteImpl)) {
  backupTimeStamp = utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(nullptr));
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
  if (pointEnabled) {
    logInfo() << "Enabling on-fault receiver output";
    ppOutputBuilder = std::make_unique<PickPointBuilder>();
    ppOutputBuilder->setMeshReader(&userMesher);
    ppOutputBuilder->setParams(seissolParameters.output.pickpointParameters);
  }
  if (elementwiseEnabled) {
    logInfo() << "Enabling 2D fault output";
    ewOutputBuilder = std::make_unique<ElementWiseBuilder>();
    ewOutputBuilder->setMeshReader(&userMesher);
    ewOutputBuilder->setParams(seissolParameters.output.elementwiseParameters);
  }
  if (!elementwiseEnabled && !pointEnabled) {
    logInfo() << "No dynamic rupture output enabled";
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
  initFaceToLtsMap();
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
    ppOutputBuilder->setLtsData(userWpTree, userWpDescr, userWpLut, userDrTree, userDrDescr);
    ppOutputBuilder->setVariableList(impl->getOutputVariables());
    ppOutputBuilder->setFaceToLtsMap(&globalFaceToLtsMap);
  }
  if (elementwiseEnabled) {
    ewOutputBuilder->setLtsData(userWpTree, userWpDescr, userWpLut, userDrTree, userDrDescr);
    ewOutputBuilder->setFaceToLtsMap(&globalFaceToLtsMap);
  }
}

void OutputManager::initElementwiseOutput() {
  ewOutputBuilder->build(ewOutputData);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  const auto& receiverPoints = ewOutputData->receiverPoints;
  const auto cellConnectivity = getCellConnectivity(receiverPoints);
  const auto faultTags = getFaultTags(receiverPoints);
  const auto vertices = getAllVertices(receiverPoints);
  constexpr auto MaxNumVars = std::tuple_size<DrVarsT>::value;
  const auto outputMask = seissolParameters.output.elementwiseParameters.outputMask;
  const auto intMask = convertMaskFromBoolToInt<MaxNumVars>(outputMask);

  const double printTime = seissolParameters.output.elementwiseParameters.printTimeIntervalSec;
  const auto backendType = seissolParameters.output.xdmfWriterBackend;

  if (seissolParameters.output.elementwiseParameters.vtkorder < 0) {
    std::vector<real*> dataPointers;
    auto recordPointers = [&dataPointers](auto& var, int) {
      if (var.isActive) {
        for (int dim = 0; dim < var.dim(); ++dim) {
          dataPointers.push_back(var.data[dim]);
        }
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
  } else {
    // Code to be refactored.

    auto order = seissolParameters.output.elementwiseParameters.vtkorder;
    if (order == 0) {
      logError() << "VTK order 0 is currently not supported for the elementwise fault output.";
    }

    io::instance::mesh::VtkHdfWriter writer("fault-elementwise",
                                            receiverPoints.size() /
                                                seissol::init::vtk2d::Shape[order][1],
                                            2,
                                            order);

    writer.addPointProjector([=](double* target, std::size_t index) {
      for (std::size_t i = 0; i < seissol::init::vtk2d::Shape[order][1]; ++i) {
        for (int j = 0; j < 3; ++j) {
          target[i * 3 + j] =
              receiverPoints[seissol::init::vtk2d::Shape[order][1] * index + i].global.coords[j];
        }
      }
    });

    misc::forEach(ewOutputData->vars, [&](auto& var, int i) {
      if (var.isActive) {
        for (int d = 0; d < var.dim(); ++d) {
          auto* data = var.data[d];
          writer.addPointData<real>(VariableLabels[i][d],
                                    std::vector<std::size_t>(),
                                    [=](real* target, std::size_t index) {
                                      std::memcpy(
                                          target,
                                          data + seissol::init::vtk2d::Shape[order][1] * index,
                                          sizeof(real) * seissol::init::vtk2d::Shape[order][1]);
                                    });
        }
      }
    });

    auto& self = *this;
    writer.addHook([&](std::size_t, double) { self.updateElementwiseOutput(); });

    io::writer::ScheduledWriter schedWriter;
    schedWriter.interval = printTime;
    schedWriter.name = "fault-elementwise";
    schedWriter.planWrite = writer.makeWriter();

    seissolInstance.getOutputManager().addOutput(schedWriter);
  }
}

void OutputManager::initPickpointOutput() {
  ppOutputBuilder->build(ppOutputData);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  if (seissolParameters.output.pickpointParameters.collectiveio) {
    logError() << "Collective IO for the Fault Pickpoint output is still under construction.";
  }

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
  for (size_t i = 0; i < outputData->receiverPoints.size(); ++i) {
    const auto& receiver = outputData->receiverPoints[i];
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

        // output coordinates
        file << "# x1\t" << makeFormatted(point[0]) << '\n';
        file << "# x2\t" << makeFormatted(point[1]) << '\n';
        file << "# x3\t" << makeFormatted(point[2]) << '\n';

        // stress info
        std::array<real, 6> rotatedInitialStress{};

        {
          auto [layer, face] = faceToLtsMap.at(receiver.faultFaceIndex);

          const auto* initialStressVar = layer->var(drDescr->initialStressInFaultCS);
          const auto* initialStress = reinterpret_cast<const real*>(initialStressVar[face]);

          seissol::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
          alignAlongDipAndStrikeKernel.stressRotationMatrix =
              outputData->stressGlbToDipStrikeAligned[i].data();
          alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
              outputData->stressFaceAlignedToGlb[i].data();

          alignAlongDipAndStrikeKernel.initialStress = initialStress;
          alignAlongDipAndStrikeKernel.rotatedStress = rotatedInitialStress.data();
          alignAlongDipAndStrikeKernel.execute();
        }

        file << "# P_0\t" << makeFormatted(rotatedInitialStress[0]) << '\n';
        file << "# T_s\t" << makeFormatted(rotatedInitialStress[3]) << '\n';
        file << "# T_d\t" << makeFormatted(rotatedInitialStress[5]) << '\n';

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

void OutputManager::initFaceToLtsMap() {
  if (drTree != nullptr) {
    const size_t readerFaultSize = meshReader->getFault().size();
    const size_t ltsFaultSize = drTree->getNumberOfCells(Ghost);

    faceToLtsMap.resize(std::max(readerFaultSize, ltsFaultSize));
    globalFaceToLtsMap.resize(faceToLtsMap.size());
    for (auto& layer : drTree->leaves(Ghost)) {

      DRFaceInformation* faceInformation = layer.var(drDescr->faceInformation);
      for (size_t ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
        faceToLtsMap[faceInformation[ltsFace].meshFace] = std::make_pair(&layer, ltsFace);
      }
    }

    DRFaceInformation* faceInformation = drTree->var(drDescr->faceInformation);
    for (size_t ltsFace = 0; ltsFace < ltsFaultSize; ++ltsFace) {
      globalFaceToLtsMap[faceInformation[ltsFace].meshFace] = ltsFace;
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
