// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PickpointOutput.h"

#include "Alignment.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Filesystem.h"
#include "Core/Common/Real.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/IO/Datatype/Inference.h"
#include "Core/IO/Instance/Point/Grouping.h"
#include "Core/IO/Instance/Point/Hdf5Table.h"
#include "Core/IO/Instance/Point/TableWriter.h"
#include "Core/IO/Writer/Writer.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
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
#include "Features/DynamicRupture/Output/Builders/PickPointBuilder.h"
#include "Features/DynamicRupture/Output/DataTypes.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Features/DynamicRupture/Output/Geometry.h"
#include "Features/DynamicRupture/Output/ReceiverBasedOutput.h"
#include "GeneratedCode/runtime.h"
#include "WavePropagation/Output/ReceiverParameters.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <limits>
#include <map>
#include <memory>
#include <mpi.h>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Core/Kernels/DeviceAux/HostCopyAux.h"
#endif

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

std::string buildIndexedMPIFileName(const std::string& namePrefix,
                                    int index,
                                    const std::string& nameSuffix) {
  std::stringstream fileName;
  fileName << namePrefix << '-' << nameSuffix << '-' << makeFormatted<int, WideFormat>(index);
  return fileName.str();
}

//! The samples of the on-fault receivers are taken on streams of their own, and have to be
//! complete before the host reads them.
void waitForSamples(
    std::unordered_map<std::size_t, std::shared_ptr<seissol::dr::ReceiverOutputData>>& outputData) {
  for (auto& [layerId, data] : outputData) {
    if (data->hostQueue != nullptr) {
      data->hostQueue->drain();
    }
    if (data->extraRuntime.has_value()) {
      data->extraRuntime->wait();
    }
  }
}

} // namespace

namespace seissol::dr::output {

PickpointOutput::PickpointOutput(std::string name,
                                 initializer::parameters::PickpointParameters parameters,
                                 ReceiverOutput& impl,
                                 RunContext& run,
                                 const initializer::parameters::DRParameters& drParameters)
    : name_(std::move(name)), parameters_(std::move(parameters)), impl_(impl), run_(run),
      drParameters_(drParameters) {}

PickpointOutput::~PickpointOutput() = default;

double PickpointOutput::interval(double stepSize) const {
  return parameters_.samplingInterval > 0 ? parameters_.samplingInterval
                                          : parameters_.printTimeInterval * stepSize;
}

void PickpointOutput::setInputParam(geometry::MeshReader& mesh, const FaultOutputSetup& setup) {
  setup_ = setup;
  logInfo() << "Enabling the on-fault receiver output" << name_.c_str();
  builder_ = std::make_unique<PickPointBuilder>();
  builder_->setMeshReader(&mesh);
  builder_->setParams(parameters_);
  builder_->setSimulationCount(setup.simulationCount, setup.fusedCount);
  // the shortest time between two samples, which sizes the cache of the samples between writes
  const double pickdt = parameters_.samplingInterval;
  builder_->setTimestep(pickdt > 0 ? pickdt : setup.minimumTimestep, setup.endTime);
}

void PickpointOutput::setLtsData(LTS::Storage& wpStorage,
                                 LTS::Backmap& wpBackmap,
                                 DynamicRupture::Storage& drStorage,
                                 initializer::StorageBackmap<1>* faceToLtsMap) {
  drStorage_ = &drStorage;
  faceToLtsMap_ = faceToLtsMap;
  builder_->setLtsData(wpStorage, wpBackmap, drStorage);
  builder_->setVariableList(impl_.getOutputVariables());
  builder_->setFaceToLtsMap(faceToLtsMap);
}

std::array<double, 6> PickpointOutput::initialStress(const ReceiverOutputData& outputData,
                                                     std::size_t index) {
  // every stress source in effect at the start, aligned with strike and dip
  const auto& receiver = outputData.receivers[index];
  // the receivers of a point are numbered contiguously, so the point is the last one whose
  // receivers start at or before this one
  const auto& receiverOffset = outputData.topology.receiverOffset;
  const auto point = static_cast<std::size_t>(std::ranges::upper_bound(receiverOffset, index) -
                                              receiverOffset.begin() - 1);
  const auto& face = outputData.topology.faces[outputData.topology.points[point].faceId];
  const auto position = faceToLtsMap_->get(receiver.faultFaceIndex.value(), receiver.copy);
  const auto sourceCount = dr::stressSourceCount(drParameters_);
  auto& drLayer = drStorage_->layer(position.color);

  std::array<double, 6> rotatedInitialStress{};
  dispatchConfig(drLayer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;
    const auto* stresses = drLayer.var<DynamicRupture::StressSourceInFaultCS>(cfg);
    const auto* onsets = drLayer.var<DynamicRupture::StressSourceOnset>(cfg);
    const auto* riseTimes = drLayer.var<DynamicRupture::StressSourceRiseTime>(cfg);
    auto unrotatedInitialStress =
        dr::stressAtTime<Cfg>(&stresses[position.cell * sourceCount],
                              &riseTimes[position.cell * sourceCount],
                              &onsets[position.cell * sourceCount],
                              sourceCount,
                              static_cast<std::uint32_t>(receiver.gpIndex),
                              static_cast<real>(0.0));

    constexpr auto Variant = configIdOf<Cfg>();
    const auto& transform = std::get<FaceTransform<Cfg>>(face.transform);
    runtime::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
    alignAlongDipAndStrikeKernel.stressRotationMatrix = runtime::init::stressRotationMatrix::view(
        Variant, transform.stressGlbToDipStrikeAligned.data());
    alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
        runtime::init::reducedFaceAlignedMatrix::view(Variant,
                                                      transform.stressFaceAlignedToGlb.data());

    std::array<real, 6> rotatedStress{};
    alignAlongDipAndStrikeKernel.initialStress =
        runtime::init::initialStress::view(Variant, unrotatedInitialStress.data());
    alignAlongDipAndStrikeKernel.rotatedStress =
        runtime::init::rotatedStress::view(Variant, rotatedStress.data());
    alignAlongDipAndStrikeKernel.execute(Variant);
    std::copy(rotatedStress.begin(), rotatedStress.end(), rotatedInitialStress.begin());
  });
  return rotatedInitialStress;
}

void PickpointOutput::init(const std::string& backupTimeStamp) {
  logInfo() << "Setting up the on-fault receivers of" << name_.c_str() << ".";
  builder_->build(outputData_);
  if (isDeviceOn() && useHostActor(run_.env())) {
    // the host actor takes the samples, from buffers that hold the derivatives, then the variables
    for (auto& [layerId, data] : outputData_) {
      if (data->deviceDataCollector == nullptr) {
        continue;
      }
      const auto aligned = [](std::size_t bytes) {
        return (bytes + Alignment - 1) / Alignment * Alignment;
      };
      // the start of the step, the derivatives, the variables
      data->dataOffset = aligned(sizeof(double));
      data->variablesOffset =
          aligned(data->dataOffset + data->deviceDataCollector->gatheredBytes());
      const auto bytes =
          data->variablesOffset +
          (data->deviceVariables != nullptr ? data->deviceVariables->gatheredBytes() : 0);
      data->hostQueue =
          std::make_unique<parallel::runtime::HostActor::Queue>(run_.hostActor(), bytes);
      data->extraRuntime.reset();
    }
  }

  writer_.enable(run_.modules(), parameters_.writeInterval);
  writer_.setupWriter([&] { takeText(); }, [&] { writeText(); });

  if (parameters_.format == initializer::parameters::ReceiverOutputFormat::Hdf5) {
    initTable();
  } else {
    initFiles(backupTimeStamp);
  }
}

void PickpointOutput::initFiles(const std::string& backupTimeStamp) {
  const auto numSimulations = setup_.simulationCount;
  // the files of the output named "faultreceivers" keep the names they have always had
  const std::string stem = name_ == "faultreceivers" ? std::string("faultreceiver") : name_;

  for (auto& [id, outputData] : outputData_) {
    // A file holds one receiver. With fused simulations, a sample takes a row per simulation, in
    // the order of their index, which a column of its own names -- the layout of the HDF5 table.
    std::stringstream variables;
    variables << "VARIABLES = \"Time\"";
    if (numSimulations > 1) {
      variables << " ,\"SimulationIndex\"";
    }
    misc::forEach(outputData->vars, [&variables](const auto& var, int i) {
      if (var.isActive) {
        for (std::size_t dim = 0; dim < var.dim(); ++dim) {
          variables << " ,\"" << VariableLabels[i][dim] << '\"';
        }
      }
    });

    std::map<std::size_t, std::vector<std::size_t>> simulationsOfReceiver;
    for (std::size_t i = 0; i < outputData->receivers.size(); ++i) {
      simulationsOfReceiver[outputData->receivers[i].globalReceiverIndex.value()].push_back(i);
    }

    auto& files = files_[id];
    files.clear();
    // (with OpenMP enabled, Clang 18 does not let a lambda capture the structured binding
    // outputData itself)
    const auto& receivers = outputData->receivers;
    for (auto& [index, simulations] : simulationsOfReceiver) {
      std::ranges::sort(simulations, [&receivers](std::size_t a, std::size_t b) {
        return receivers[a].simIndex < receivers[b].simIndex;
      });

      auto fileName =
          buildIndexedMPIFileName(run_.outputPrefix(), static_cast<int>(index + 1), stem);
      seissol::generateBackupFileIfNecessary(fileName, "dat", {backupTimeStamp});
      fileName += ".dat";
      files.push_back(File{.fileName = fileName, .indices = simulations});

      if (std::filesystem::exists(fileName)) {
        continue;
      }
      std::ofstream file(fileName, std::ios_base::out);
      if (!file.is_open()) {
        logError() << "Cannot open fault receiver file" << fileName;
      }

      const auto& receiver = outputData->receivers[simulations.front()];
      const auto& point = receiver.global;
      file << "TITLE = \"Temporal Signal for fault receiver number " << (index + 1) << "\"\n";
      file << variables.str() << '\n';
      file << "# x1\t" << makeFormatted(point[0]) << '\n';
      file << "# x2\t" << makeFormatted(point[1]) << '\n';
      file << "# x3\t" << makeFormatted(point[2]) << '\n';
      file << "# face-global-id\t" << receiver.globalFaultFaceId() << '\n';
      file << "# plus-cell-global-id\t" << receiver.elementGlobalIndex.value() << '\n';
      // sides are int8_t; print them as numbers, not as characters
      file << "# plus-face-side\t" << static_cast<int>(receiver.localFaceSideId.value()) << '\n';
      file << "# minus-cell-global-id\t" << receiver.elementNeighborGlobalIndex.value() << '\n';
      file << "# minus-face-side\t" << static_cast<int>(receiver.localNeighborFaceSideId.value())
           << '\n';

      // a value per simulation, in the order of the rows of a sample
      std::vector<std::array<double, 6>> stresses;
      stresses.reserve(simulations.size());
      for (const auto simulation : simulations) {
        stresses.push_back(initialStress(*outputData, simulation));
      }
      using namespace misc::quantity_indices;
      const std::array<std::pair<const char*, QuantityIndices>, 3> components{
          std::pair{"P_0", XX}, std::pair{"T_s", XY}, std::pair{"T_d", XZ}};
      for (const auto& [label, component] : components) {
        file << "# " << label;
        for (const auto& stress : stresses) {
          file << '\t' << makeFormatted(stress[component]);
        }
        file << '\n';
      }
      file << "# run-id\t" << run_.runId() << '\n';
    }
  }
}

void PickpointOutput::initTable() {
  // A receiver of the fault output is one point of one simulation, so what it records is the
  // time plus the active components -- the per-simulation column names the text files carry are
  // a row of their own here.
  std::vector<io::instance::point::TableQuantity> quantitySet;
  quantitySet.push_back(io::instance::point::TableQuantity{
      .name = "Time", .datatype = io::datatype::inferDatatype<double>()});
  // A rank without on-fault receivers has no point to describe; the table learns the quantity
  // sets of the other ranks when it groups the points.
  if (!outputData_.empty()) {
    misc::forEach(outputData_.begin()->second->vars, [&](const auto& var, int i) {
      if (var.isActive) {
        for (std::size_t dim = 0; dim < var.dim(); ++dim) {
          quantitySet.push_back(io::instance::point::TableQuantity{
              .name = VariableLabels[i][dim], .datatype = io::datatype::inferDatatype<double>()});
        }
      }
    });
  }

  // the rows of a rank are its receivers by the number they have in the receiver file, and the
  // simulations of one of them next to each other
  tableRows_.clear();
  for (const auto& [layerId, outputData] : outputData_) {
    for (std::size_t point = 0; point < outputData->receivers.size(); ++point) {
      tableRows_.emplace_back(layerId, point);
    }
  }
  std::ranges::sort(tableRows_, [this](const auto& a, const auto& b) {
    const auto& left = outputData_.at(a.first)->receivers[a.second];
    const auto& right = outputData_.at(b.first)->receivers[b.second];
    return std::make_tuple(left.globalReceiverIndex.value(), left.simIndex) <
           std::make_tuple(right.globalReceiverIndex.value(), right.simIndex);
  });

  const std::vector<std::vector<io::instance::point::TableQuantity>> pointQuantities(
      tableRows_.size(), quantitySet);
  table_ = std::make_unique<io::instance::point::Hdf5Table>(
      name_, pointQuantities, seissol::Mpi::mpi.comm(), parameters_.samplechunk);

  // what the header of a text file states about a receiver, as a column each
  std::vector<std::uint64_t> receiverIds;
  std::vector<std::uint64_t> simulations;
  std::vector<std::uint64_t> faceIds;
  std::vector<std::int64_t> plusCells;
  std::vector<std::int64_t> plusSides;
  std::vector<std::int64_t> minusCells;
  std::vector<std::int64_t> minusSides;
  std::vector<double> coordinates;
  for (const auto& [layerId, point] : tableRows_) {
    const auto& receiver = outputData_.at(layerId)->receivers[point];
    receiverIds.push_back(static_cast<std::uint64_t>(receiver.globalReceiverIndex.value()));
    simulations.push_back(static_cast<std::uint64_t>(receiver.simIndex));
    faceIds.push_back(static_cast<std::uint64_t>(receiver.globalFaultFaceId()));
    plusCells.push_back(static_cast<std::int64_t>(receiver.elementGlobalIndex.value()));
    plusSides.push_back(receiver.localFaceSideId.value());
    minusCells.push_back(static_cast<std::int64_t>(receiver.elementNeighborGlobalIndex.value()));
    minusSides.push_back(receiver.localNeighborFaceSideId.value());
    for (std::size_t dimension = 0; dimension < Cell::Dim; ++dimension) {
      coordinates.push_back(receiver.global[dimension]);
    }
  }
  table_->addPointData("ReceiverId", {}, receiverIds);
  table_->addPointData("SimulationIndex", {}, simulations);
  table_->addPointData("FaceId", {}, faceIds);
  table_->addPointData("PlusCellId", {}, plusCells);
  table_->addPointData("PlusFaceSide", {}, plusSides);
  table_->addPointData("MinusCellId", {}, minusCells);
  table_->addPointData("MinusFaceSide", {}, minusSides);
  table_->addPointData("Coordinates", {Cell::Dim}, coordinates);

  // the stress the fault starts out under, as the P_0, T_s and T_d of the header of a text file
  std::vector<double> initialStresses;
  initialStresses.reserve(tableRows_.size() * 3);
  for (const auto& [layerId, point] : tableRows_) {
    const auto stress = initialStress(*outputData_.at(layerId), point);
    using namespace misc::quantity_indices;
    initialStresses.push_back(stress[XX]);
    initialStresses.push_back(stress[XY]);
    initialStresses.push_back(stress[XZ]);
  }
  table_->addPointData("InitialStress", {3}, initialStresses);

  io::writer::ScheduledWriter scheduled;
  scheduled.name = name_;
  scheduled.interval = parameters_.writeInterval;
  // the fault receivers take their samples into buffers of their own
  scheduled.hostData = {};
  // the samples get taken into the table at the synchronization point, before the next epoch
  // adds to the cache
  scheduled.capture = [this](double /*time*/) { collectSamples(); };
  scheduled.planWrite = table_->makeWriter();
  run_.outputManager().addOutput(scheduled);
}

void PickpointOutput::recordInitial(double time, parallel::runtime::StreamRuntime& runtime) {
  for (const auto& [id, _] : outputData_) {
    // at the start of the expansion the cells hold before the first step; the step size only scales
    // the time basis, and has to be positive for that
    record(id, time, 1.0, 0.0, runtime);
    finishStep(id, true, runtime);
  }
}

bool PickpointOutput::hasReceivers(std::size_t layerId) const {
  return outputData_.contains(layerId);
}

parallel::runtime::StreamRuntime* PickpointOutput::runtime(std::size_t layerId) {
  const auto& outputData = outputData_.at(layerId);
  return outputData->extraRuntime.has_value() ? &outputData->extraRuntime.value() : nullptr;
}

parallel::runtime::StreamRuntime&
    PickpointOutput::gatherData(std::size_t layerId,
                                [[maybe_unused]] const double* deviceClock,
                                [[maybe_unused]] double* stepStart,
                                [[maybe_unused]] double* stepStartOnDevice,
                                [[maybe_unused]] double hostTime,
                                parallel::runtime::StreamRuntime& runtime) {
  const auto& outputData = outputData_.at(layerId);
#ifdef ACL_DEVICE
  if (outputData->extraRuntime.has_value()) {
    runtime.eventSync(outputData->extraRuntime->eventRecord());
  }
  // by a kernel, as the copy engines may be busy with a large copy of another stream
  kernels::device::aux::hostcopy::copyToHost(
      stepStartOnDevice, deviceClock, sizeof(double), runtime.stream());
  ReceiverOutput::gatherFaultOutput(outputData, runtime);
  if (outputData->extraRuntime.has_value()) {
    outputData->extraRuntime->eventSync(runtime.eventRecord());
  }
#else
  *stepStart = hostTime;
#endif
  return outputData->extraRuntime.has_value() ? outputData->extraRuntime.value() : runtime;
}

bool PickpointOutput::postAtRunTime(std::size_t layerId,
                                    const double* deviceClock,
                                    parallel::runtime::StreamRuntime& runtime,
                                    std::function<void(double)> record) {
  const auto& outputData = outputData_.at(layerId);
  if (!isDeviceOn() || outputData->hostQueue == nullptr) {
    return false;
  }
  const auto staging = ReceiverOutput::gatherForHost(outputData, runtime, true, deviceClock);
  outputData->hostQueue->post(
      runtime, staging, [outputData, record = std::move(record)](const void* data) {
        outputData->gathered = data;
        std::invoke(record, *static_cast<const double*>(data));
        outputData->gathered = nullptr;
      });
  return true;
}

void PickpointOutput::evaluate(std::size_t layerId, double time, double stepSize, double inStep) {
  const auto& outputData = outputData_.at(layerId);
  ensureCacheLevel(*outputData);
  impl_.evaluateFaultOutput(initializer::parameters::OutputType::AtPickpoint,
                            drParameters_.slipRateOutputType,
                            outputData,
                            time,
                            stepSize,
                            inStep);
}

void PickpointOutput::record(std::size_t layerId,
                             double time,
                             double stepSize,
                             double inStep,
                             parallel::runtime::StreamRuntime& runtime) {
  if (!hasReceivers(layerId)) {
    return;
  }
  const auto& outputData = outputData_.at(layerId);
  if (outputData->currentCacheLevel >= outputData->maxCacheLevel) {
    // On CUDA and SYCL, calcFaultOutput leaves the samples to a host function on a stream, which
    // runs once the friction law is done, and stores into whatever memory the cache has at that
    // time. Growing the cache moves it, and a sample stored into the old memory while it is being
    // moved is lost. Hence, all samples taken so far have to be in first. Likewise for the host
    // actor.
    if (outputData->hostQueue != nullptr) {
      outputData->hostQueue->drain();
    } else {
      auto& sampleRuntime =
          outputData->extraRuntime.has_value() ? outputData->extraRuntime.value() : runtime;
      sampleRuntime.wait();
    }
  }
  ensureCacheLevel(*outputData);

  impl_.calcFaultOutput(initializer::parameters::OutputType::AtPickpoint,
                        drParameters_.slipRateOutputType,
                        outputData,
                        runtime,
                        time,
                        stepSize,
                        inStep);
}

void PickpointOutput::finishStep(std::size_t layerId,
                                 bool sampled,
                                 parallel::runtime::StreamRuntime& runtime) {
  if (!hasReceivers(layerId)) {
    return;
  }
  const auto& outputData = outputData_.at(layerId);
  if (isDeviceOn() && outputData->hostQueue != nullptr) {
    // the friction state at the end of the step, in a buffer of its own for the host actor
    const auto staging = ReceiverOutput::gatherForHost(outputData, runtime, false);
    outputData->hostQueue->post(runtime, staging, [this, outputData](const void* data) {
      outputData->gathered = data;
      impl_.keepStepStart(outputData);
      outputData->gathered = nullptr;
    });
    return;
  }
  if constexpr (isDeviceOn()) {
    if (!sampled) {
      // in the order of the streams that calcFaultOutput() keeps
      if (outputData->extraRuntime.has_value()) {
        runtime.eventSync(outputData->extraRuntime->eventRecord());
      }
      ReceiverOutput::gatherFaultState(outputData, runtime);
      if (outputData->extraRuntime.has_value()) {
        outputData->extraRuntime->eventSync(runtime.eventRecord());
      }
    }
  }
  auto& callRuntime =
      outputData->extraRuntime.has_value() ? outputData->extraRuntime.value() : runtime;
  callRuntime.enqueueHost([this, outputData] { impl_.keepStepStart(outputData); });
}

void PickpointOutput::finishStepNow(std::size_t layerId) {
  impl_.keepStepStart(outputData_.at(layerId));
}

void PickpointOutput::collectSamples() {
  waitForSamples(outputData_);
  const auto& grouping = table_->grouping();

  // How far a table grows has to be the same everywhere, and the clusters a rank holds do not all
  // cache the same number of samples between two writes.
  std::vector<std::size_t> samples(grouping.groupCount(), 0);
  for (std::size_t row = 0; row < tableRows_.size(); ++row) {
    const auto group = grouping.group[row];
    samples[group] =
        std::max(samples[group], outputData_.at(tableRows_[row].first)->currentCacheLevel);
  }
  if (!samples.empty()) {
    seissol::Mpi::mpi.allreduceContainer(samples, MPI_MAX);
  }

  std::vector<double*> storage(grouping.groupCount(), nullptr);
  for (std::size_t group = 0; group < grouping.groupCount(); ++group) {
    auto* prepared = table_->prepare(group, samples[group]);
    storage[group] = reinterpret_cast<double*>(prepared);
    // A receiver that cached fewer samples than the longest one of its table leaves the rest of
    // its column unset, and a zero there is a value a reader cannot tell from a measurement.
    const auto values = samples[group] * table_->localPointCount(group) *
                        table_->sampleSize(group) / sizeof(double);
    std::fill_n(storage[group], values, std::numeric_limits<double>::quiet_NaN());
  }

  for (std::size_t row = 0; row < tableRows_.size(); ++row) {
    // named rather than bound, since the lambda below captures the point, which Clang 18 does not
    // allow for a structured binding with OpenMP enabled
    const auto layerId = tableRows_[row].first;
    const auto point = tableRows_[row].second;
    auto& outputData = *outputData_.at(layerId);
    const auto group = grouping.group[row];
    const auto column = table_->localRow(row);
    const auto points = table_->localPointCount(group);
    const auto components = table_->sampleSize(group) / sizeof(double);

    for (std::size_t level = 0; level < std::min(outputData.currentCacheLevel, samples[group]);
         ++level) {
      auto* target = storage[group] + (level * points + column) * components;
      std::size_t position = 0;
      target[position++] = outputData.cachedTime[level];
      misc::forEach(outputData.vars, [&](const auto& var, int) {
        if (var.isActive) {
          for (std::size_t dim = 0; dim < var.dim(); ++dim) {
            target[position++] = var(dim, level, point);
          }
        }
      });
    }
  }

  for (auto& [layerId, outputData] : outputData_) {
    outputData->currentCacheLevel = 0;
  }
}

void PickpointOutput::takeText() {
  if (table_ != nullptr) {
    // the tables are filled and handed over by the writer they are registered with
    return;
  }

  waitForSamples(outputData_);
  const auto numSimulations = setup_.simulationCount;
  for (auto& [layerId, outputData] : outputData_) {
    for (const auto& file : files_.at(layerId)) {
      std::stringstream data;
      for (size_t level = 0; level < outputData->currentCacheLevel; ++level) {
        // a row per simulation, in the order the file holds them in; the output variables are
        // indexed per receiver
        for (const std::size_t receiverId : file.indices) {
          data << makeFormatted(outputData->cachedTime[level]) << '\t';
          if (numSimulations > 1) {
            data << outputData->receivers[receiverId].simIndex << '\t';
          }
          auto recordResults = [receiverId, level, &data](const auto& var, int) {
            if (var.isActive) {
              for (std::size_t dim = 0; dim < var.dim(); ++dim) {
                data << makeFormatted(var(dim, level, receiverId)) << '\t';
              }
            }
          };
          misc::forEach(outputData->vars, recordResults);
          data << '\n';
        }
      }

      takenText_.emplace_back(file.fileName, data.str());
    }
    outputData->currentCacheLevel = 0;
  }
}

void PickpointOutput::writeText() {
  for (const auto& [fileName, data] : takenText_) {
    std::ofstream file(fileName, std::ios_base::app);
    if (file.is_open()) {
      file << data;
    } else {
      logError() << "cannot open" << fileName;
    }
    file.close();
  }
  takenText_.clear();
}

void PickpointOutput::ensureCacheLevel(ReceiverOutputData& outputData) {
  if (outputData.currentCacheLevel >= outputData.maxCacheLevel) {
    // our calculation was off (maybe due to many intermediate sync points), so resize

    outputData.maxCacheLevel = outputData.currentCacheLevel + 1;
    const auto newCacheLevel = outputData.maxCacheLevel;
    outputData.cachedTime.resize(newCacheLevel);
    misc::forEach(outputData.vars,
                  [newCacheLevel](auto& var, int) { var.resizeCache(newCacheLevel); });
  }
}

} // namespace seissol::dr::output
