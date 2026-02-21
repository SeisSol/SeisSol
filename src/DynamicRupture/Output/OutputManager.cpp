// SPDX-FileCopyrightText: 2022 SeisSol Group
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
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Writer/Writer.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Runtime/Stream.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
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
  if (index >= 0) {
    suffix << nameSuffix << '-' << makeFormatted<int, WideFormat>(index);
  } else {
    suffix << nameSuffix << "-r" << makeFormatted<int, WideFormat>(seissol::Mpi::mpi.rank());
  }
  return buildFileName(namePrefix, suffix.str(), fileExtension);
}

} // namespace

namespace seissol::dr::output {

OutputManager::OutputManager(std::unique_ptr<ReceiverOutput> concreteImpl,
                             seissol::SeisSol& seissolInstance)
    : seissolInstance(seissolInstance), ewOutputData(std::make_shared<ReceiverOutputData>()),
      impl(std::move(concreteImpl)) {
  backupTimeStamp = utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(nullptr));
}

OutputManager::~OutputManager() = default;

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
    ppOutputBuilder->setTimestep(seissolInstance.getMemoryManager().clusterLayout().minimumTimestep,
                                 seissolParameters.timeStepping.endTime);
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

void OutputManager::setLtsData(LTS::Storage& userWpStorage,
                               LTS::Backmap& userWpBackmap,
                               DynamicRupture::Storage& userDrStorage) {
  wpStorage = &userWpStorage;
  wpBackmap = &userWpBackmap;
  drStorage = &userDrStorage;
  impl->setLtsData(userWpStorage, userWpBackmap, userDrStorage);
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
    ppOutputBuilder->setLtsData(userWpStorage, userWpBackmap, userDrStorage);
    ppOutputBuilder->setVariableList(impl->getOutputVariables());
    ppOutputBuilder->setFaceToLtsMap(&globalFaceToLtsMap);
  }
  if (elementwiseEnabled) {
    ewOutputBuilder->setLtsData(userWpStorage, userWpBackmap, userDrStorage);
    ewOutputBuilder->setFaceToLtsMap(&globalFaceToLtsMap);
  }
}

void OutputManager::initElementwiseOutput() {
  logInfo() << "Setting up the fault output.";
  ewOutputBuilder->build(ewOutputData);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  const auto& receiverPoints = ewOutputData->receiverPoints;
  const auto cellConnectivity = getCellConnectivity(receiverPoints);
  const auto faultTags = getFaultTags(receiverPoints);
  const auto vertices = getAllVertices(receiverPoints);
  constexpr auto MaxNumVars = std::tuple_size_v<DrVarsT>;
  const auto outputMask = seissolParameters.output.elementwiseParameters.outputMask;
  const auto intMask = convertMaskFromBoolToInt<MaxNumVars>(outputMask);

  const double printTime = seissolParameters.output.elementwiseParameters.printTimeIntervalSec;

  const auto orderPre = seissolParameters.output.elementwiseParameters.vtkorder;

  const auto order = static_cast<uint32_t>(std::max(0, orderPre));

  const auto pointCount =
      io::instance::geometry::numPoints(order, io::instance::geometry::Shape::Triangle);
  const auto dataCount = std::max(pointCount, static_cast<std::size_t>(1));

  const auto config = io::instance::geometry::WriterConfig{
      order,
      orderPre < 0 ? io::instance::geometry::WriterFormat::Xdmf
                   : io::instance::geometry::WriterFormat::Vtk,
      seissolParameters.output.xdmfWriterBackend ==
              seissol::initializer::parameters::XdmfBackend::Posix
          ? io::instance::geometry::WriterBackend::Binary
          : io::instance::geometry::WriterBackend::Hdf5,
      io::instance::geometry::WriterGroup::FullSnapshot,
      seissolParameters.output.hdfcompress};

  auto writer = io::instance::geometry::GeometryWriter("fault-elementwise",
                                                       receiverPoints.size() / pointCount /
                                                           multisim::NumSimulations,
                                                       io::instance::geometry::Shape::Triangle,
                                                       config,
                                                       1);

  writer.addPointProjector([=](double* target, std::size_t index, std::size_t) {
    for (std::size_t i = 0; i < pointCount; ++i) {
      for (std::size_t j = 0; j < Cell::Dim; ++j) {
        target[i * Cell::Dim + j] =
            receiverPoints[(pointCount * index + i) * multisim::NumSimulations].global.coords[j];
      }
    }
  });

  const auto rank = seissol::Mpi::mpi.rank();
  writer.addCellData<int>(
      "partition", {}, true, [=](int* target, std::size_t, std::size_t) { target[0] = rank; });

  writer.addCellData<int>(
      "fault-tag", {}, true, [=, &receiverPoints](int* target, std::size_t index, std::size_t) {
        *target = receiverPoints[index * multisim::NumSimulations].faultTag;
      });

  writer.addCellData<std::size_t>(
      "global-id",
      {},
      true,
      [=, &receiverPoints](std::size_t* target, std::size_t index, std::size_t) {
        *target = receiverPoints[index * multisim::NumSimulations].elementGlobalIndex * 4 +
                  receiverPoints[index * multisim::NumSimulations].localFaceSideId;
      });

  misc::forEach(ewOutputData->vars, [&](auto& var, int i) {
    if (var.isActive) {
      for (std::size_t d = 0; d < var.dim(); ++d) {
        auto* data = var.data[d];
        const auto variableName = [&](std::size_t d, std::size_t s) {
          if constexpr (multisim::MultisimEnabled) {
            return VariableLabels[i][d] + "_" + std::to_string(s);
          } else {
            return VariableLabels[i][d];
          }
        };
        for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
          writer.addGeometryOutput<real>(
              variableName(d, s),
              std::vector<std::size_t>(),
              false,
              [=](real* target, std::size_t index, std::size_t) {
                for (std::size_t i = 0; i < dataCount; ++i) {
                  target[index + i] = data[(dataCount * index + i) * multisim::NumSimulations + s];
                }
              });
        }
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

void OutputManager::initPickpointOutput() {
  logInfo() << "Setting up on-fault receivers.";
  ppOutputBuilder->build(ppOutputData);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  seissolInstance.pickpointWriter().enable(
      seissolParameters.output.pickpointParameters.writeInterval);
  seissolInstance.pickpointWriter().setupWriter([&]() { flushPickpointDataToFile(); });

  if (seissolParameters.output.pickpointParameters.collectiveio) {
    logError() << "Collective IO for the on-fault receiver output is still under construction.";
  }

  for (auto& [id, outputData] : ppOutputData) {
    const bool allReceiversInOneFilePerRank =
        seissolParameters.output.pickpointParameters.aggregate;
    auto& files = ppFiles[id];

    if (allReceiversInOneFilePerRank) {
      // aggregate all receivers per rank

      files.resize(1);
      auto fileName = buildIndexedMPIFileName(seissolParameters.output.prefix, -1, "faultreceiver");
      fileName += ".dat";
      std::vector<std::size_t> receivers(outputData->receiverPoints.size());
      std::iota(receivers.begin(), receivers.end(), 0);
      files[0] = PickpointFile{fileName, receivers};
    } else {
      // aggregate at least all fused simulations

      std::unordered_map<std::size_t, std::vector<std::size_t>> globalIndexMap;
      for (size_t i = 0; i < outputData->receiverPoints.size(); ++i) {
        globalIndexMap[outputData->receiverPoints[i].globalReceiverIndex].push_back(i);
      }

      files.resize(globalIndexMap.size());
      std::size_t counter = 0;
      for (const auto& [index, receivers] : globalIndexMap) {
        auto fileName =
            buildIndexedMPIFileName(seissolParameters.output.prefix, index + 1, "faultreceiver");
        seissol::generateBackupFileIfNecessary(fileName, "dat", {backupTimeStamp});
        fileName += ".dat";

        files[counter] = PickpointFile{fileName, receivers};
        ++counter;
      }
    }

    std::stringstream baseHeader;

    auto suffix = [&allReceiversInOneFilePerRank](auto pointIndex, auto simIndex) {
      std::string suffix;

      if (allReceiversInOneFilePerRank) {
        suffix += "-" + std::to_string(pointIndex);
      }

      if constexpr (seissol::multisim::MultisimEnabled) {
        suffix += "-" + std::to_string(simIndex);
      }

      return suffix;
    };

    const size_t actualPointCount =
        allReceiversInOneFilePerRank ? outputData->receiverPoints.size() / multisim::NumSimulations
                                     : 1;

    for (std::size_t pointIndex = 0; pointIndex < actualPointCount; ++pointIndex) {
      for (std::size_t simIndex = 0; simIndex < multisim::NumSimulations; ++simIndex) {
        size_t labelCounter = 0;
        auto collectVariableNames =
            [&baseHeader, &labelCounter, &simIndex, &pointIndex, suffix](auto& var, int i) {
              if (var.isActive) {
                for (std::size_t dim = 0; dim < var.dim(); ++dim) {
                  baseHeader << " ,\"" << VariableLabels[i][dim]
                             << suffix(pointIndex + 1, simIndex + 1) << '\"';
                  ++labelCounter;
                }
              } else {
                labelCounter += var.dim();
              }
            };
        misc::forEach(outputData->vars, collectVariableNames);
      }
    }

    for (size_t i = 0; i < files.size(); ++i) {
      const auto& ppfile = files[i];

      if (!seissol::filesystem::exists(ppfile.fileName)) {
        std::ofstream file(ppfile.fileName, std::ios_base::out);
        if (file.is_open()) {
          std::stringstream title;

          title << "TITLE = \"Temporal Signal for fault receiver number(s) and simulation(s)";
          for (const auto& gIdx : ppfile.indices) {
            const auto& receiver = outputData->receiverPoints[gIdx];
            const size_t globalIndex = receiver.globalReceiverIndex + 1;
            const size_t simIndex = receiver.simIndex + 1;
            title << " " << globalIndex << "," << simIndex << ";";
          }
          title << "\"";

          file << title.str() << '\n';
          file << "VARIABLES = \"Time\"";

          file << baseHeader.str();

          file << '\n';

          for (const auto& gIdx : ppfile.indices) {
            const auto& receiver = outputData->receiverPoints[gIdx];
            const size_t globalIndex = receiver.globalReceiverIndex + 1;
            const size_t simIndex = receiver.simIndex;
            const auto& point = receiver.global;

            // output coordinates
            if (simIndex == 0) {
              file << "# Receiver number " << globalIndex << '\n';
              file << "# x1\t" << makeFormatted(point[0]) << '\n';
              file << "# x2\t" << makeFormatted(point[1]) << '\n';
              file << "# x3\t" << makeFormatted(point[2]) << '\n';
            }

            // stress info
            std::array<real, 6> rotatedInitialStress{};
            {
              auto [layer, face] = faceToLtsMap.at(receiver.faultFaceIndex);

              const auto* initialStressVar = layer->var<DynamicRupture::InitialStressInFaultCS>();
              const auto* initialStress = initialStressVar[face];
              std::array<real, 6> unrotatedInitialStress{};
              for (std::size_t stressVar = 0; stressVar < unrotatedInitialStress.size();
                   ++stressVar) {
                unrotatedInitialStress[stressVar] = initialStress[stressVar][receiver.gpIndex];
              }

              seissol::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
              alignAlongDipAndStrikeKernel.stressRotationMatrix =
                  outputData->stressGlbToDipStrikeAligned[i].data();
              alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
                  outputData->stressFaceAlignedToGlb[i].data();

              alignAlongDipAndStrikeKernel.initialStress = unrotatedInitialStress.data();
              alignAlongDipAndStrikeKernel.rotatedStress = rotatedInitialStress.data();
              alignAlongDipAndStrikeKernel.execute();
            }

            file << "# P_0" << simIndex + 1 << "\t" << makeFormatted(rotatedInitialStress[0])
                 << '\n';
            file << "# T_s" << simIndex + 1 << "\t" << makeFormatted(rotatedInitialStress[3])
                 << '\n';
            file << "# T_d" << simIndex + 1 << "\t" << makeFormatted(rotatedInitialStress[5])
                 << '\n';
          }
        } else {
          logError() << "cannot open " << ppfile.fileName;
        }
        file.close();
      }
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
  if (drStorage != nullptr) {
    ::seissol::initializer::StorageBackmap<1> backmap;
    backmap.setSize(meshReader->getFault().size());

    faceToLtsMap.resize(meshReader->getFault().size());
    globalFaceToLtsMap.resize(faceToLtsMap.size());

    const auto* globalFaceInformation = drStorage->var<DynamicRupture::FaceInformation>();
    for (auto& layer : drStorage->leaves()) {
      const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
      for (size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        faceToLtsMap[faceInformation[ltsFace].meshFace] = std::make_pair(&layer, ltsFace);
        backmap.addElement(layer.id(),
                           globalFaceInformation,
                           faceInformation,
                           faceInformation[ltsFace].meshFace,
                           ltsFace);
      }
    }

    for (size_t i = 0; i < meshReader->getFault().size(); ++i) {
      globalFaceToLtsMap[i] = backmap.get(i);
    }
  }
  impl->setFaceToLtsMap(&faceToLtsMap);
}

bool OutputManager::isAtPickpoint(double time, double dt) {
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();
  const bool isFirstStep = iterationStep == 0;
  const double abortTime = seissolParameters.timeStepping.endTime;
  const bool isCloseToTimeOut = (abortTime - time) < (dt * TimeMargin);

  const int printTimeInterval = seissolParameters.output.pickpointParameters.printTimeInterval;
  const bool isOutputIteration = iterationStep % printTimeInterval == 0;

  return (isFirstStep || isOutputIteration || isCloseToTimeOut);
}

void OutputManager::writePickpointOutput(std::size_t layerId,
                                         double time,
                                         double dt,
                                         double meshDt,
                                         double meshInDt,
                                         parallel::runtime::StreamRuntime& runtime) {
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();
  if (this->ppOutputBuilder) {
    if (this->isAtPickpoint(time, dt)) {
      const auto findResult = ppOutputData.find(layerId);
      if (findResult != ppOutputData.end()) {
        const auto& outputData = findResult->second;

        if (outputData->currentCacheLevel >= outputData->maxCacheLevel) {
          logError() << "Error: not enough space for on-fault receiver allocated."
                     << outputData->maxCacheLevel;
        }

        impl->calcFaultOutput(seissol::initializer::parameters::OutputType::AtPickpoint,
                              seissolParameters.drParameters.slipRateOutputType,
                              outputData,
                              runtime,
                              time,
                              meshDt,
                              meshInDt);
      }
    }
    ++iterationStep;
  }
}

void OutputManager::writePickpointOutput(double time, double dt) {
  for (const auto& [id, _] : ppOutputData) {
    writePickpointOutput(id, time, dt, 0, 1, runtime);
  }
}

void OutputManager::flushPickpointDataToFile() {
  for (auto& [layerId, outputData] : ppOutputData) {
    for (const auto& ppfile : ppFiles.at(layerId)) {
      std::stringstream data;
      for (size_t level = 0; level < outputData->currentCacheLevel; ++level) {
        data << makeFormatted(outputData->cachedTime[level]) << '\t';
        for (std::size_t pointId : ppfile.indices) {
          auto recordResults = [pointId, level, &data](auto& var, int) {
            if (var.isActive) {
              for (std::size_t dim = 0; dim < var.dim(); ++dim) {
                data << makeFormatted(var(dim, level, pointId)) << '\t';
              }
            }
          };
          misc::forEach(outputData->vars, recordResults);
        }
        data << '\n';
      }

      std::ofstream file(ppfile.fileName, std::ios_base::app);
      if (file.is_open()) {
        file << data.str();
      } else {
        logError() << "cannot open" << ppfile.fileName;
      }
      file.close();
    }
    outputData->currentCacheLevel = 0;
  }
}

void OutputManager::updateElementwiseOutput() {
  if (this->ewOutputBuilder) {
    const auto& seissolParameters = seissolInstance.getSeisSolParameters();
    impl->calcFaultOutput(seissol::initializer::parameters::OutputType::Elementwise,
                          seissolParameters.drParameters.slipRateOutputType,
                          ewOutputData,
                          runtime);
    runtime.wait();
  }
}
} // namespace seissol::dr::output
