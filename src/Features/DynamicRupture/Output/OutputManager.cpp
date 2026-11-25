// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Features/DynamicRupture/Output/OutputManager.h"

#include "Alignment.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Filesystem.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/IO/Instance/Geometry/Geometry.h"
#include "Core/IO/Instance/Geometry/Typedefs.h"
#include "Core/IO/Writer/Writer.h"
#include "Core/IO/WriterParameters.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/Helper.h"
#include "Core/Parallel/MPI.h"
#include "Core/Parallel/Runtime/HostActor.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Run/Context.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Output/Builders/ElementWiseBuilder.h"
#include "Features/DynamicRupture/Output/DataTypes.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Features/DynamicRupture/Output/Geometry.h"
#include "Features/DynamicRupture/Output/OutputAux.h"
#include "Features/DynamicRupture/Output/PickpointOutput.h"
#include "Features/DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "WavePropagation/Output/DofSync.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <utils/timeutils.h>
#include <vector>

namespace seissol::dr::output {

OutputManager::OutputManager(std::unique_ptr<ReceiverOutput> concreteImpl,
                             RunContext& run,
                             const seissol::initializer::parameters::DRParameters& drParameters,
                             FaultWriterSettings writer,
                             writer::DofSync& dofSync)
    : run_(run), drParameters_(drParameters), writer_(writer), dofSync_(dofSync),
      ewOutputData_(std::make_shared<ReceiverOutputData>()), impl_(std::move(concreteImpl)) {
  backupTimeStamp_ = utils::TimeUtils::timeAsString(BackupTimeStampFormat, time(nullptr));
}

OutputManager::~OutputManager() = default;

void OutputManager::setOutputs(
    const seissol::initializer::parameters::ElementwiseFaultParameters* elementwise,
    const std::vector<
        std::pair<std::string, seissol::initializer::parameters::PickpointParameters>>&
        pickpoints) {
  elementwise_.reset();
  if (elementwise != nullptr) {
    elementwise_ = *elementwise;
  }
  pickpointOutputs_.clear();
  for (const auto& [name, parameters] : pickpoints) {
    pickpointOutputs_.push_back(
        std::make_unique<PickpointOutput>(name, parameters, *impl_, run_, drParameters_));
  }
}

std::vector<PickpointOutput*> OutputManager::pickpointOutputs() {
  std::vector<PickpointOutput*> outputs;
  outputs.reserve(pickpointOutputs_.size());
  for (const auto& output : pickpointOutputs_) {
    outputs.push_back(output.get());
  }
  return outputs;
}

void OutputManager::setInputParam(seissol::geometry::MeshReader& userMesher,
                                  const FaultOutputSetup& setup) {
  using namespace initializer;
  meshReader_ = &userMesher;

  impl_->setMeshReader(&userMesher);

  impl_->setDrParameters(&drParameters_);
  setup_ = setup;
  for (auto& output : pickpointOutputs_) {
    output->setInputParam(userMesher, setup);
  }
  if (elementwise_.has_value()) {
    logInfo() << "Enabling 2D fault output";
    ewOutputBuilder_ = std::make_unique<ElementWiseBuilder>();
    ewOutputBuilder_->setMeshReader(&userMesher);
    ewOutputBuilder_->setParams(elementwise_.value());
    ewOutputBuilder_->setSimulationCount(setup.simulationCount, setup.fusedCount);
  }
  if (!elementwise_.has_value() && pickpointOutputs_.empty()) {
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
  for (auto& output : pickpointOutputs_) {
    output->setLtsData(userWpStorage, userWpBackmap, userDrStorage, &faceToLtsMap_);
  }
  if (elementwise_.has_value()) {
    ewOutputBuilder_->setLtsData(userWpStorage, userWpBackmap, userDrStorage);
    ewOutputBuilder_->setFaceToLtsMap(&faceToLtsMap_);
  }
}

namespace {
//! @brief The file grouping a time series mode asks the writer for; none if there is no mode.
std::optional<io::instance::geometry::WriterGroup>
    writerGroupOf(std::optional<seissol::initializer::parameters::TimeSeriesMode> mode) {
  if (!mode.has_value()) {
    return std::nullopt;
  }
  switch (mode.value()) {
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
  // the receivers of the simulations of the run follow each other
  const auto numSimulations = setup_.simulationCount;

  const auto& receivers = ewOutputData_->receivers;
  // the receivers in the order of the output: cell by cell, then point by point, and simulation by
  // simulation
  const auto& outputOrder = ewOutputData_->outputOrder;

  const auto& elementwise = elementwise_.value();
  const double writeInterval = elementwise.printTimeIntervalSec;

  const auto orderPre = elementwise.vtkorder;

  const auto order = static_cast<uint32_t>(std::max(0, orderPre));

  const auto dataCount =
      io::instance::geometry::numPoints(order, io::instance::geometry::Shape::Triangle);
  const auto pointCount = io::instance::geometry::numPoints(
      std::max(order, 1U), io::instance::geometry::Shape::Triangle);

  const auto format = orderPre < 0 ? io::instance::geometry::WriterFormat::Xdmf
                                   : io::instance::geometry::WriterFormat::Vtk;

  const auto config = io::instance::geometry::WriterConfig{
      .order = order,
      .format = format,
      .backend = writer_.backend == seissol::initializer::parameters::XdmfBackend::Posix
                     ? io::instance::geometry::WriterBackend::Binary
                     : io::instance::geometry::WriterBackend::Hdf5,
      .time = io::instance::geometry::supportedWriterGroup(
          writerGroupOf(elementwise.timeSeries), format, "elementwise fault"),
      .compress = writer_.compress};

  auto writer = io::instance::geometry::GeometryWriter(
      "fault",
      receivers.size() / dataCount / numSimulations,
      io::instance::geometry::Shape::Triangle,
      config,
      1,

      [=, &receivers, &outputOrder](double* target, std::size_t index, std::size_t) {
        if (order > 0) {
          for (std::size_t i = 0; i < pointCount; ++i) {
            for (std::size_t j = 0; j < Cell::Dim; ++j) {
              target[i * Cell::Dim + j] =
                  receivers[outputOrder[(pointCount * index + i) * numSimulations]].global[j];
            }
          }
        } else {
          const auto& triangle = receivers[outputOrder[index * numSimulations]].globalTriangle;
          for (std::size_t i = 0; i < pointCount; ++i) {
            for (std::size_t j = 0; j < Cell::Dim; ++j) {
              target[i * Cell::Dim + j] = triangle.point(i)[j];
            }
          }
        }
      });

  const auto rank = seissol::Mpi::mpi.rank();
  writer.addCellData<int>(
      "partition", {}, true, [=](int* target, std::size_t, std::size_t) { target[0] = rank; });

  writer.addCellData<int>(
      "fault-tag",
      {},
      true,
      [=, &receivers, &outputOrder](int* target, std::size_t index, std::size_t) {
        *target = faultTagOfCell(receivers, outputOrder, index, dataCount, numSimulations);
      });

  writer.addCellData<std::size_t>(
      "global-id",
      {},
      true,
      [=, &receivers, &outputOrder](std::size_t* target, std::size_t index, std::size_t) {
        *target = globalFaceIdOfCell(receivers, outputOrder, index, dataCount, numSimulations);
      });

  misc::forEach(ewOutputData_->vars, [&](const auto& var, int i) {
    if (var.isActive) {
      for (std::size_t d = 0; d < var.dim(); ++d) {
        const auto* data = var[d];
        const auto variableName = [&](std::size_t d, std::size_t s) {
          if (numSimulations > 1) {
            return VariableLabels[i][d] + "-" + std::to_string(s);
          } else {
            return VariableLabels[i][d];
          }
        };
        for (std::size_t s = 0; s < numSimulations; ++s) {
          writer.addGeometryOutput<double>(
              variableName(d, s),
              std::vector<std::size_t>(),
              false,
              [=, &outputOrder](double* target, std::size_t index, std::size_t) {
                for (std::size_t i = 0; i < dataCount; ++i) {
                  target[i] = data[outputOrder[(dataCount * index + i) * numSimulations + s]];
                }
              });
        }
      }
    }
  });

  auto& self = *this;
  writer.addHook([&](std::size_t, double time) {
    self.dofSync_.syncDofs(time);
    self.evaluateElementwiseOutput(time);
  });

  io::writer::ScheduledWriter schedWriter;
  schedWriter.interval = writeInterval;
  schedWriter.name = "fault";
  schedWriter.planWrite = writer.makeWriter();
  // the data of the faces comes from the device at the synchronization point itself
  schedWriter.capture = [&self](double /*time*/) { self.gatherElementwiseOutput(); };
  // the faces gather their data from the device themselves; the halo exchange of the unknowns
  // before (dofSync) reads those of the cells
  schedWriter.hostData = {};
  schedWriter.hostData.add<LTS::Dofs>(*wpStorage_).addVariables(*drStorage_);

  run_.outputManager().addOutput(schedWriter);
}

void OutputManager::init() {
  if (ewOutputBuilder_) {
    initElementwiseOutput();
  }
  for (auto& output : pickpointOutputs_) {
    output->init(backupTimeStamp_);
  }
}

void OutputManager::initFaceToLtsMap() {
  if (drStorage_ != nullptr) {
    faceToLtsMap_.setSize(meshReader_->getFault().size(), wpBackmap_->copyCount());

    const auto* globalFaceInformation = drStorage_->var<DynamicRupture::FaceInformation>();
    for (auto& layer : drStorage_->leaves()) {
      const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
      for (size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {

        faceToLtsMap_.addElement(layer.id(),
                                 globalFaceInformation,
                                 faceInformation,
                                 faceInformation[ltsFace].meshFace,
                                 faceInformation[ltsFace].copy,
                                 ltsFace);
      }
    }
  }
  impl_->setFaceToLtsMap(&faceToLtsMap_);
}

void OutputManager::recordInitialPickpointOutput(double time) {
  for (auto& output : pickpointOutputs_) {
    output->recordInitial(time, runtime_);
  }
}

void OutputManager::gatherElementwiseOutput() {
  if constexpr (isDeviceOn()) {
    if (this->ewOutputBuilder_) {
      ReceiverOutput::gatherFaultOutput(ewOutputData_, runtime_);
      runtime_.wait();
    }
  }
}

void OutputManager::evaluateElementwiseOutput(double time) {
  if (this->ewOutputBuilder_) {
    // at a synchronization point every cluster has just completed a time step ending here, so this
    // is also the time the stored friction state belongs to: the end of the step
    impl_->evaluateFaultOutput(seissol::initializer::parameters::OutputType::Elementwise,
                               drParameters_.slipRateOutputType,
                               ewOutputData_,
                               time,
                               1.0,
                               1.0);
  }
}
} // namespace seissol::dr::output
