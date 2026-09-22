// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/Output/OutputManager.h"

#include "Common/Constants.h"
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
#include "IO/Datatype/Inference.h"
#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Point/Grouping.h"
#include "IO/Instance/Point/Hdf5Table.h"
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
#include <limits>
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
    : seissolInstance_(seissolInstance), ewOutputData_(std::make_shared<ReceiverOutputData>()),
      impl_(std::move(concreteImpl)) {
  backupTimeStamp_ = utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(nullptr));
}

OutputManager::~OutputManager() = default;

void OutputManager::setInputParam(seissol::geometry::MeshReader& userMesher) {
  using namespace initializer;
  meshReader_ = &userMesher;

  impl_->setMeshReader(&userMesher);

  const auto& seissolParameters = seissolInstance_.parameters();
  impl_->setDrParameters(&seissolParameters.drParameters);

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
    ppOutputBuilder_ = std::make_unique<PickPointBuilder>();
    ppOutputBuilder_->setMeshReader(&userMesher);
    ppOutputBuilder_->setParams(seissolParameters.output.pickpointParameters);
    ppOutputBuilder_->setTimestep(seissolInstance_.memoryManager().clusterLayout().minimumTimestep,
                                  seissolParameters.timeStepping.endTime);
  }
  if (elementwiseEnabled) {
    logInfo() << "Enabling 2D fault output";
    ewOutputBuilder_ = std::make_unique<ElementWiseBuilder>();
    ewOutputBuilder_->setMeshReader(&userMesher);
    ewOutputBuilder_->setParams(seissolParameters.output.elementwiseParameters);
  }
  if (!elementwiseEnabled && !pointEnabled) {
    logInfo() << "No dynamic rupture output enabled";
  }
}

void OutputManager::setLtsData(LTS::Storage& userWpStorage,
                               LTS::Backmap& userWpBackmap,
                               DynamicRupture::Storage& userDrStorage) {
  wpStorage_ = &userWpStorage;
  wpBackmap_ = &userWpBackmap;
  drStorage_ = &userDrStorage;
  impl_->setLtsData(userWpStorage, userWpBackmap, userDrStorage);
  initFaceToLtsMap();
  const auto& seissolParameters = seissolInstance_.parameters();
  const bool bothEnabled = seissolParameters.drParameters.outputPointType ==
                           seissol::initializer::parameters::OutputType::AtPickpointAndElementwise;
  const bool pointEnabled = seissolParameters.drParameters.outputPointType ==
                                seissol::initializer::parameters::OutputType::AtPickpoint ||
                            bothEnabled;
  const bool elementwiseEnabled = seissolParameters.drParameters.outputPointType ==
                                      seissol::initializer::parameters::OutputType::Elementwise ||
                                  bothEnabled;
  if (pointEnabled) {
    ppOutputBuilder_->setLtsData(userWpStorage, userWpBackmap, userDrStorage);
    ppOutputBuilder_->setVariableList(impl_->getOutputVariables());
    ppOutputBuilder_->setFaceToLtsMap(&faceToLtsMap_);
  }
  if (elementwiseEnabled) {
    ewOutputBuilder_->setLtsData(userWpStorage, userWpBackmap, userDrStorage);
    ewOutputBuilder_->setFaceToLtsMap(&faceToLtsMap_);
  }
}

namespace {
//! @brief The file grouping a time series mode asks the writer for.
io::instance::geometry::WriterGroup
    writerGroupOf(seissol::initializer::parameters::TimeSeriesMode mode) {
  switch (mode) {
  case seissol::initializer::parameters::TimeSeriesMode::Incremental:
    return io::instance::geometry::WriterGroup::IncrementalSnapshot;
  case seissol::initializer::parameters::TimeSeriesMode::Monolith:
    return io::instance::geometry::WriterGroup::Monolith;
  case seissol::initializer::parameters::TimeSeriesMode::Snapshot:
    return io::instance::geometry::WriterGroup::FullSnapshot;
  }
  return io::instance::geometry::WriterGroup::FullSnapshot;
}
} // namespace

void OutputManager::initElementwiseOutput() {
  logInfo() << "Setting up the fault output.";
  ewOutputBuilder_->build(ewOutputData_);
  const auto& seissolParameters = seissolInstance_.parameters();

  const auto& receiverPoints = ewOutputData_->receiverPoints;

  const double writeInterval = seissolParameters.output.elementwiseParameters.printTimeIntervalSec;

  const auto orderPre = seissolParameters.output.elementwiseParameters.vtkorder;

  const auto order = static_cast<uint32_t>(std::max(0, orderPre));

  const auto dataCount =
      io::instance::geometry::numPoints(order, io::instance::geometry::Shape::Triangle);
  const auto pointCount = io::instance::geometry::numPoints(
      std::max(order, 1U), io::instance::geometry::Shape::Triangle);

  const auto format = orderPre < 0 ? io::instance::geometry::WriterFormat::Xdmf
                                   : io::instance::geometry::WriterFormat::Vtk;

  const auto config = io::instance::geometry::WriterConfig{
      order,
      format,
      seissolParameters.output.xdmfWriterBackend ==
              seissol::initializer::parameters::XdmfBackend::Posix
          ? io::instance::geometry::WriterBackend::Binary
          : io::instance::geometry::WriterBackend::Hdf5,
      io::instance::geometry::supportedWriterGroup(
          writerGroupOf(seissolParameters.output.elementwiseParameters.timeSeries),
          format,
          "elementwise fault"),
      seissolParameters.output.hdfcompress};

  auto writer = io::instance::geometry::GeometryWriter(
      "fault",
      receiverPoints.size() / dataCount / multisim::NumSimulations,
      io::instance::geometry::Shape::Triangle,
      config,
      1,

      [=](double* target, std::size_t index, std::size_t) {
        if (order > 0) {
          for (std::size_t i = 0; i < pointCount; ++i) {
            for (std::size_t j = 0; j < Cell::Dim; ++j) {
              target[i * Cell::Dim + j] =
                  receiverPoints[(pointCount * index + i) * multisim::NumSimulations]
                      .global.coords[j];
            }
          }
        } else {
          const auto& triangle = receiverPoints[index * multisim::NumSimulations].globalTriangle;
          for (std::size_t i = 0; i < pointCount; ++i) {
            for (std::size_t j = 0; j < Cell::Dim; ++j) {
              target[i * Cell::Dim + j] = triangle.point(i).coords[j];
            }
          }
        }
      });

  const auto rank = seissol::Mpi::mpi.rank();
  writer.addCellData<int>(
      "partition", {}, true, [=](int* target, std::size_t, std::size_t) { target[0] = rank; });

  writer.addCellData<int>(
      "fault-tag", {}, true, [=, &receiverPoints](int* target, std::size_t index, std::size_t) {
        *target = faultTagOfCell(receiverPoints, index, dataCount, multisim::NumSimulations);
      });

  writer.addCellData<std::size_t>(
      "global-id",
      {},
      true,
      [=, &receiverPoints](std::size_t* target, std::size_t index, std::size_t) {
        *target = globalFaceIdOfCell(receiverPoints, index, dataCount, multisim::NumSimulations);
      });

  misc::forEach(ewOutputData_->vars, [&](const auto& var, int i) {
    if (var.isActive) {
      for (std::size_t d = 0; d < var.dim(); ++d) {
        const auto* data = var[d];
        const auto variableName = [&](std::size_t d, std::size_t s) {
          if constexpr (multisim::MultisimEnabled) {
            return VariableLabels[i][d] + "-" + std::to_string(s);
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
                  target[i] = data[(dataCount * index + i) * multisim::NumSimulations + s];
                }
              });
        }
      }
    }
  });

  auto& self = *this;
  writer.addHook([&](std::size_t, double time) {
    self.seissolInstance_.dofSync().syncDofs(time);
    self.updateElementwiseOutput();
  });

  io::writer::ScheduledWriter schedWriter;
  schedWriter.interval = writeInterval;
  schedWriter.name = "fault";
  schedWriter.planWrite = writer.makeWriter();

  seissolInstance_.outputManager().addOutput(schedWriter);
}

void OutputManager::initPickpointOutput() {
  logInfo() << "Setting up on-fault receivers.";
  ppOutputBuilder_->build(ppOutputData_);
  const auto& seissolParameters = seissolInstance_.parameters();

  seissolInstance_.pickpointWriter().enable(
      seissolParameters.output.pickpointParameters.writeInterval);
  seissolInstance_.pickpointWriter().setupWriter([&]() { flushPickpointDataToFile(); });

  if (seissolParameters.output.pickpointParameters.format ==
      seissol::initializer::parameters::ReceiverOutputFormat::Hdf5) {
    initPickpointTable();
    return;
  }

  if (seissolParameters.output.pickpointParameters.collectiveio) {
    logError() << "Collective IO for the on-fault receiver output is still under construction.";
  }

  for (auto& [id, outputData] : ppOutputData_) {
    const bool allReceiversInOneFilePerRank =
        seissolParameters.output.pickpointParameters.aggregate;
    auto& files = ppFiles_[id];

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
        seissol::generateBackupFileIfNecessary(fileName, "dat", {backupTimeStamp_});
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
            [&baseHeader, &labelCounter, &simIndex, &pointIndex, suffix](const auto& var, int i) {
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
            const size_t simIndex = receiver.simIndex + 1;
            const auto& point = receiver.global;

            // output coordinates
            if (simIndex == 1) {
              file << "# Receiver number " << globalIndex << '\n';
              file << "# x1\t" << makeFormatted(point[0]) << '\n';
              file << "# x2\t" << makeFormatted(point[1]) << '\n';
              file << "# x3\t" << makeFormatted(point[2]) << '\n';
              file << "# face-global-id\t" << receiver.globalFaultFaceId() << '\n';
              file << "# plus-cell-global-id\t" << receiver.elementGlobalIndex << '\n';
              file << "# plus-face-side\t" << receiver.localFaceSideId << '\n';
              file << "# minus-cell-global-id\t" << receiver.elementNeighborGlobalIndex << '\n';
              file << "# minus-face-side\t" << receiver.localNeighborFaceSideId << '\n';
            }

            // stress info
            std::array<real, 6> rotatedInitialStress{};
            {
              const auto position = faceToLtsMap_.get(receiver.faultFaceIndex);

              const auto* initialStress =
                  drStorage_->lookup<DynamicRupture::InitialStressInFaultCS>(position);
              std::array<real, 6> unrotatedInitialStress{};
              for (std::size_t stressVar = 0; stressVar < unrotatedInitialStress.size();
                   ++stressVar) {
                unrotatedInitialStress[stressVar] = initialStress[stressVar][receiver.gpIndex];
              }

              seissol::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
              alignAlongDipAndStrikeKernel.stressRotationMatrix =
                  outputData->stressGlbToDipStrikeAligned[gIdx].data();
              alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
                  outputData->stressFaceAlignedToGlb[gIdx].data();

              alignAlongDipAndStrikeKernel.initialStress = unrotatedInitialStress.data();
              alignAlongDipAndStrikeKernel.rotatedStress = rotatedInitialStress.data();
              alignAlongDipAndStrikeKernel.execute();
            }

            {
              using namespace misc::quantity_indices;
              file << "# P_0" << simIndex << "\t" << makeFormatted(rotatedInitialStress[XX])
                   << '\n';
              file << "# T_s" << simIndex << "\t" << makeFormatted(rotatedInitialStress[XY])
                   << '\n';
              file << "# T_d" << simIndex << "\t" << makeFormatted(rotatedInitialStress[XZ])
                   << '\n';
            }
          }
        } else {
          logError() << "Cannot open fault receiver file" << ppfile.fileName;
        }
        file.close();
      }
    }
  }
}

void OutputManager::init() {
  if (ewOutputBuilder_) {
    initElementwiseOutput();
  }
  if (ppOutputBuilder_) {
    initPickpointOutput();
  }
}

void OutputManager::initFaceToLtsMap() {
  if (drStorage_ != nullptr) {
    faceToLtsMap_.setSize(meshReader_->getFault().size());

    const auto* globalFaceInformation = drStorage_->var<DynamicRupture::FaceInformation>();
    for (auto& layer : drStorage_->leaves()) {
      const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
      for (size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {

        faceToLtsMap_.addElement(layer.id(),
                                 globalFaceInformation,
                                 faceInformation,
                                 faceInformation[ltsFace].meshFace,
                                 ltsFace);
      }
    }
  }
  impl_->setFaceToLtsMap(&faceToLtsMap_);
}

bool OutputManager::isAtPickpoint(double time, double dt) {
  const auto& seissolParameters = seissolInstance_.parameters();
  const bool isFirstStep = iterationStep_ == 0;
  const double abortTime = seissolParameters.timeStepping.endTime;
  const bool isCloseToTimeOut = (abortTime - time) < (dt * TimeMargin);

  const int printTimeInterval = seissolParameters.output.pickpointParameters.printTimeInterval;
  const bool isOutputIteration = iterationStep_ % printTimeInterval == 0;

  return (isFirstStep || isOutputIteration || isCloseToTimeOut);
}

void OutputManager::writePickpointOutput(std::size_t layerId,
                                         double time,
                                         double dt,
                                         double meshDt,
                                         double meshInDt,
                                         parallel::runtime::StreamRuntime& runtime) {
  const auto& seissolParameters = seissolInstance_.parameters();
  if (this->ppOutputBuilder_) {
    if (this->isAtPickpoint(time, dt)) {
      const auto findResult = ppOutputData_.find(layerId);
      if (findResult != ppOutputData_.end()) {
        const auto& outputData = findResult->second;

        if (outputData->currentCacheLevel >= outputData->maxCacheLevel) {
          // our calculation was off (maybe due to many intermediate sync points), so resize

          outputData->maxCacheLevel = outputData->currentCacheLevel + 1;
          const auto newCacheLevel = outputData->maxCacheLevel;
          outputData->cachedTime.resize(newCacheLevel);
          misc::forEach(outputData->vars,
                        [newCacheLevel](auto& var, int) { var.resizeCache(newCacheLevel); });
        }

        impl_->calcFaultOutput(seissol::initializer::parameters::OutputType::AtPickpoint,
                               seissolParameters.drParameters.slipRateOutputType,
                               outputData,
                               runtime,
                               time,
                               meshDt,
                               meshInDt);
      }
    }
    ++iterationStep_;
  }
}

void OutputManager::writePickpointOutput(double time, double dt) {
  for (const auto& [id, _] : ppOutputData_) {
    writePickpointOutput(id, time, dt, 0, 1, runtime_);
  }
}

void OutputManager::initPickpointTable() {
  const auto& seissolParameters = seissolInstance_.parameters();

  // A receiver of the fault output is one point of one simulation, so what it records is the
  // time plus the active components -- the per-simulation column names the text files carry are
  // a row of their own here.
  std::vector<io::instance::point::TableQuantity> quantitySet;
  quantitySet.push_back(
      io::instance::point::TableQuantity{"Time", io::datatype::inferDatatype<real>()});
  // A rank without on-fault receivers has no point to describe; the table learns the quantity
  // sets of the other ranks when it groups the points.
  if (!ppOutputData_.empty()) {
    misc::forEach(ppOutputData_.begin()->second->vars, [&](const auto& var, int i) {
      if (var.isActive) {
        for (std::size_t dim = 0; dim < var.dim(); ++dim) {
          quantitySet.push_back(io::instance::point::TableQuantity{
              VariableLabels[i][dim], io::datatype::inferDatatype<real>()});
        }
      }
    });
  }

  // the rows of a rank are its receivers by the number they have in the receiver file, and the
  // simulations of one of them next to each other
  ppTableRows_.clear();
  for (const auto& [layerId, outputData] : ppOutputData_) {
    for (std::size_t point = 0; point < outputData->receiverPoints.size(); ++point) {
      ppTableRows_.emplace_back(layerId, point);
    }
  }
  std::sort(ppTableRows_.begin(), ppTableRows_.end(), [this](const auto& a, const auto& b) {
    const auto& left = ppOutputData_.at(a.first)->receiverPoints[a.second];
    const auto& right = ppOutputData_.at(b.first)->receiverPoints[b.second];
    return std::tie(left.globalReceiverIndex, left.simIndex) <
           std::tie(right.globalReceiverIndex, right.simIndex);
  });

  const std::vector<std::vector<io::instance::point::TableQuantity>> pointQuantities(
      ppTableRows_.size(), quantitySet);
  ppTable_ = std::make_unique<io::instance::point::Hdf5Table>(
      "faultreceivers",
      pointQuantities,
      seissol::Mpi::mpi.comm(),
      seissolParameters.output.pickpointParameters.samplechunk);

  // what the header of a text file states about a receiver, as a column each
  std::vector<std::uint64_t> receiverIds;
  std::vector<std::uint64_t> simulations;
  std::vector<std::uint64_t> faceIds;
  std::vector<std::int64_t> plusCells;
  std::vector<std::int64_t> plusSides;
  std::vector<std::int64_t> minusCells;
  std::vector<std::int64_t> minusSides;
  std::vector<double> coordinates;
  for (const auto& [layerId, point] : ppTableRows_) {
    const auto& receiver = ppOutputData_.at(layerId)->receiverPoints[point];
    receiverIds.push_back(static_cast<std::uint64_t>(receiver.globalReceiverIndex));
    simulations.push_back(static_cast<std::uint64_t>(receiver.simIndex));
    faceIds.push_back(static_cast<std::uint64_t>(receiver.globalFaultFaceId()));
    plusCells.push_back(static_cast<std::int64_t>(receiver.elementGlobalIndex));
    plusSides.push_back(receiver.localFaceSideId);
    minusCells.push_back(static_cast<std::int64_t>(receiver.elementNeighborGlobalIndex));
    minusSides.push_back(receiver.localNeighborFaceSideId);
    for (std::size_t dimension = 0; dimension < Cell::Dim; ++dimension) {
      coordinates.push_back(receiver.global.coords[dimension]);
    }
  }
  ppTable_->addPointData("ReceiverId", {}, receiverIds);
  ppTable_->addPointData("SimulationIndex", {}, simulations);
  ppTable_->addPointData("FaceId", {}, faceIds);
  ppTable_->addPointData("PlusCellId", {}, plusCells);
  ppTable_->addPointData("PlusFaceSide", {}, plusSides);
  ppTable_->addPointData("MinusCellId", {}, minusCells);
  ppTable_->addPointData("MinusFaceSide", {}, minusSides);
  ppTable_->addPointData("Coordinates", {Cell::Dim}, coordinates);

  io::writer::ScheduledWriter scheduled;
  scheduled.name = "faultreceivers";
  scheduled.interval = seissolParameters.output.pickpointParameters.writeInterval;
  scheduled.planWrite = [this, plan = ppTable_->makeWriter()](
                            const std::string& prefix, std::size_t counter, double time) {
    collectPickpointSamples();
    return plan(prefix, counter, time);
  };
  seissolInstance_.outputManager().addOutput(scheduled);
}

void OutputManager::collectPickpointSamples() {
  const auto& grouping = ppTable_->grouping();

  // How far a table grows has to be the same everywhere, and the clusters a rank holds do not all
  // cache the same number of samples between two writes.
  std::vector<std::size_t> samples(grouping.groupCount(), 0);
  for (std::size_t row = 0; row < ppTableRows_.size(); ++row) {
    const auto group = grouping.group[row];
    samples[group] =
        std::max(samples[group], ppOutputData_.at(ppTableRows_[row].first)->currentCacheLevel);
  }
  if (!samples.empty()) {
    MPI_Allreduce(MPI_IN_PLACE,
                  samples.data(),
                  static_cast<int>(samples.size()),
                  seissol::Mpi::castToMpiType<std::size_t>(),
                  MPI_MAX,
                  seissol::Mpi::mpi.comm());
  }

  std::vector<real*> storage(grouping.groupCount(), nullptr);
  for (std::size_t group = 0; group < grouping.groupCount(); ++group) {
    auto* prepared = ppTable_->prepare(group, samples[group]);
    storage[group] = reinterpret_cast<real*>(prepared);
    // A receiver that cached fewer samples than the longest one of its table leaves the rest of
    // its column unset, and a zero there is a value a reader cannot tell from a measurement.
    const auto values = samples[group] * ppTable_->localPointCount(group) *
                        ppTable_->sampleSize(group) / sizeof(real);
    std::fill_n(storage[group], values, std::numeric_limits<real>::quiet_NaN());
  }

  for (std::size_t row = 0; row < ppTableRows_.size(); ++row) {
    const auto [layerId, point] = ppTableRows_[row];
    auto& outputData = *ppOutputData_.at(layerId);
    const auto group = grouping.group[row];
    const auto column = ppTable_->localRow(row);
    const auto points = ppTable_->localPointCount(group);
    const auto components = ppTable_->sampleSize(group) / sizeof(real);

    for (std::size_t level = 0; level < std::min(outputData.currentCacheLevel, samples[group]);
         ++level) {
      auto* target = storage[group] + (level * points + column) * components;
      std::size_t position = 0;
      target[position++] = static_cast<real>(outputData.cachedTime[level]);
      misc::forEach(outputData.vars, [&](const auto& var, int) {
        if (var.isActive) {
          for (std::size_t dim = 0; dim < var.dim(); ++dim) {
            target[position++] = var(dim, level, point);
          }
        }
      });
    }
  }

  for (auto& [layerId, outputData] : ppOutputData_) {
    outputData->currentCacheLevel = 0;
  }
}

void OutputManager::flushPickpointDataToFile() {
  if (ppTable_ != nullptr) {
    // the tables are filled and handed over by the writer they are registered with
    return;
  }

  for (auto& [layerId, outputData] : ppOutputData_) {
    for (const auto& ppfile : ppFiles_.at(layerId)) {
      std::stringstream data;
      for (size_t level = 0; level < outputData->currentCacheLevel; ++level) {
        data << makeFormatted(outputData->cachedTime[level]) << '\t';
        for (std::size_t pointId : ppfile.indices) {
          auto recordResults = [pointId, level, &data](const auto& var, int) {
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
  if (this->ewOutputBuilder_) {
    const auto& seissolParameters = seissolInstance_.parameters();
    impl_->calcFaultOutput(seissol::initializer::parameters::OutputType::Elementwise,
                           seissolParameters.drParameters.slipRateOutputType,
                           ewOutputData_,
                           runtime_);
    runtime_.wait();
  }
}
} // namespace seissol::dr::output
