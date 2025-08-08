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
#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Writer/Writer.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"
#include <GeneratedCode/init.h>
#include <GeneratedCode/kernel.h>
#include <IO/Instance/Geometry/Typedefs.h>
#include <Memory/Descriptor/LTS.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/MultipleSimulations.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
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
    suffix << nameSuffix << "-r" << makeFormatted<int, WideFormat>(seissol::MPI::mpi.rank());
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

OutputManager::~OutputManager() {
  for (const auto& [id, _] : ppOutputData) {
    flushPickpointDataToFile(id);
  }
}

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
  constexpr auto MaxNumVars = std::tuple_size<DrVarsT>::value;
  const auto outputMask = seissolParameters.output.elementwiseParameters.outputMask;
  const auto intMask = convertMaskFromBoolToInt<MaxNumVars>(outputMask);

  const double printTime = seissolParameters.output.elementwiseParameters.printTimeIntervalSec;
  const auto backendType = seissolParameters.output.xdmfWriterBackend;

  auto order = seissolParameters.output.elementwiseParameters.vtkorder;

  // explicitly take the first config here
  using Cfg0 = Config0;
  const auto pointCount = order > 1 ? seissol::init::vtk2d<Cfg0>::Shape[order][1]
                                    : seissol::init::vtk2d<Cfg0>::Shape[1][1];
  const auto dataCount = order > 1 ? seissol::init::vtk2d<Cfg0>::Shape[order][1] : 1;

  io::instance::geometry::GeometryWriter writer("fault-elementwise",
                                                receiverPoints.size() / pointCount,
                                                io::instance::geometry::Shape::Triangle,
                                                order);

  writer.addPointProjector([=](double* target, std::size_t index) {
    for (std::size_t i = 0; i < pointCount; ++i) {
      for (int j = 0; j < 3; ++j) {
        target[i * 3 + j] = receiverPoints[pointCount * index + i].global.coords[j];
      }
    }
  });

  writer.addCellData<int>(
      "fault-tag", {}, true, [=, &receiverPoints](int* target, std::size_t index) {
        *target = receiverPoints[index].faultTag;
      });

  writer.addCellData<std::size_t>(
      "global-id", {}, true, [=, &receiverPoints](std::size_t* target, std::size_t index) {
        *target =
            receiverPoints[index].elementGlobalIndex * 4 + receiverPoints[index].localFaceSideId;
      });

  misc::forEach(ewOutputData->vars, [&](auto& var, int i) {
    if (var.isActive) {
      for (int d = 0; d < var.dim(); ++d) {
        auto* data = var.data[d];
        writer.addGeometryOutput<float>(VariableLabels[i][d],
                                        std::vector<std::size_t>(),
                                        false,
                                        [=](float* target, std::size_t index) {
                                          std::copy_n(data + dataCount * index, dataCount, target);
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

void OutputManager::initPickpointOutput() {
  logInfo() << "Setting up on-fault receivers.";
  ppOutputBuilder->build(ppOutputData);
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();

  if (seissolParameters.output.pickpointParameters.collectiveio) {
    logError() << "Collective IO for the on-fault receiver output is still under construction.";
  }

  std::size_t maxSims = 1;
  for (const auto& element : seissolInstance.meshReader().getElements()) {
    std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          maxSims = std::max(maxSims, Cfg::NumSimulations);
        },
        ConfigVariantList[element.configId]);
  }
  MPI_Allreduce(MPI_IN_PLACE,
                &maxSims,
                1,
                seissol::MPI::castToMpiType<std::size_t>(),
                MPI_MAX,
                seissol::MPI::mpi.comm());

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

    auto suffix = [&allReceiversInOneFilePerRank, maxSims](auto pointIndex, auto simIndex) {
      std::string suffix;

      if (allReceiversInOneFilePerRank) {
        suffix += "-" + std::to_string(pointIndex);
      }

      if (maxSims > 1) {
        suffix += "-" + std::to_string(simIndex);
      }

      return suffix;
    };

    const size_t actualPointCount =
        allReceiversInOneFilePerRank ? outputData->receiverPoints.size() / maxSims : 1;

    for (std::size_t pointIndex = 0; pointIndex < actualPointCount; ++pointIndex) {
      for (std::size_t simIndex = 0; simIndex < maxSims; ++simIndex) {
        auto collectVariableNames = [&baseHeader, &simIndex, &pointIndex, suffix](auto& var,
                                                                                  int index) {
          if (var.isActive) {
            for (int dim = 0; dim < var.dim(); ++dim) {
              baseHeader << " ,\"" << VariableLabels[index][dim]
                         << suffix(pointIndex + 1, simIndex + 1) << '\"';
            }
          }
        };
        misc::forEach(outputData->vars, collectVariableNames);
      }
    }

    for (size_t i = 0; i < files.size(); ++i) {
      const auto& receiver = outputData->receiverPoints[i];

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
            const auto& point = const_cast<ExtVrtxCoords&>(receiver.global);

            // output coordinates
            if (simIndex == 0) {
              file << "# Receiver number " << globalIndex << '\n';
              file << "# x1\t" << makeFormatted(point[0]) << '\n';
              file << "# x2\t" << makeFormatted(point[1]) << '\n';
              file << "# x3\t" << makeFormatted(point[2]) << '\n';
            }

            // stress info
            auto [layer, face] = faceToLtsMap.at(receiver.faultFaceIndex);
            layer->wrap([&](auto cfg) {
              using Cfg = decltype(cfg);
              using real = Real<Cfg>;

              std::array<real, 6> rotatedInitialStress{};
              {
                auto [layer, face] = faceToLtsMap.at(receiver.faultFaceIndex);

                const auto* initialStressVar =
                    layer->var<DynamicRupture::InitialStressInFaultCS>(Cfg());
                const auto* initialStress = initialStressVar[face];
                std::array<real, 6> unrotatedInitialStress{};
                for (std::size_t stressVar = 0; stressVar < unrotatedInitialStress.size();
                     ++stressVar) {
                  unrotatedInitialStress[stressVar] = initialStress[stressVar][receiver.gpIndex];
                }

                const auto& transformData =
                    std::get<TransformData<Cfg>>(outputData->transformData[i]);

                seissol::dynamicRupture::kernel::rotateInitStress<Cfg> alignAlongDipAndStrikeKernel;
                alignAlongDipAndStrikeKernel.stressRotationMatrix =
                    transformData.stressGlbToDipStrikeAligned.data();
                alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
                    transformData.stressFaceAlignedToGlb.data();

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
            });
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
    const size_t readerFaultSize = meshReader->getFault().size();
    const size_t ltsFaultSize = drStorage->size(Ghost);

    faceToLtsMap.resize(std::max(readerFaultSize, ltsFaultSize));
    globalFaceToLtsMap.resize(faceToLtsMap.size());
    for (auto& layer : drStorage->leaves(Ghost)) {

      const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
      for (size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        faceToLtsMap[faceInformation[ltsFace].meshFace] = std::make_pair(&layer, ltsFace);
      }
    }

    const auto* faceInformation = drStorage->var<DynamicRupture::FaceInformation>();
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

void OutputManager::writePickpointOutput(int64_t clusterId,
                                         double time,
                                         double dt,
                                         parallel::runtime::StreamRuntime& runtime) {
  const auto& seissolParameters = seissolInstance.getSeisSolParameters();
  if (this->ppOutputBuilder) {
    if (this->isAtPickpoint(time, dt)) {
      const auto findResult = ppOutputData.find(clusterId);
      if (findResult != ppOutputData.end()) {
        const auto& outputData = findResult->second;

        impl->calcFaultOutput(seissol::initializer::parameters::OutputType::AtPickpoint,
                              seissolParameters.drParameters.slipRateOutputType,
                              outputData,
                              runtime,
                              time);

        const bool isMaxCacheLevel =
            outputData->currentCacheLevel >=
            static_cast<size_t>(seissolParameters.output.pickpointParameters.maxPickStore);
        const bool isCloseToEnd = (seissolParameters.timeStepping.endTime - time) < dt * timeMargin;

        if (isMaxCacheLevel || isCloseToEnd) {
          // we need to wait for all data to be (internally) written to write it out
          auto& callRuntime =
              outputData->extraRuntime.has_value() ? outputData->extraRuntime.value() : runtime;
          callRuntime.wait();

          this->flushPickpointDataToFile(clusterId);
        }
      }
    }
    ++iterationStep;
  }
}

void OutputManager::writePickpointOutput(double time, double dt) {
  for (const auto& [id, _] : ppOutputData) {
    writePickpointOutput(id, time, dt, runtime);
  }
  runtime.wait();
}

void OutputManager::flushPickpointDataToFile(int64_t clusterId) {
  auto& outputData = ppOutputData.at(clusterId);

  for (const auto& ppfile : ppFiles.at(clusterId)) {
    std::stringstream data;
    for (size_t level = 0; level < outputData->currentCacheLevel; ++level) {
      data << makeFormatted(outputData->cachedTime[level]) << '\t';
      for (std::size_t pointId : ppfile.indices) {
        auto recordResults = [pointId, level, &data](auto& var, int) {
          if (var.isActive) {
            for (int dim = 0; dim < var.dim(); ++dim) {
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
      logError() << "cannot open " << ppfile.fileName;
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
                          ewOutputData,
                          runtime);
    runtime.wait();
  }
}
} // namespace seissol::dr::output
