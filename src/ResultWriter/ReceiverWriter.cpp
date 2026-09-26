// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: John Rekoske

#include "ReceiverWriter.h"

#include "Equations/Datastructures.h"
#include "Geometry/MeshReader.h"
#include "IO/Datatype/Inference.h"
#include "IO/Instance/Point/Grouping.h"
#include "IO/Instance/Point/Hdf5Table.h"
#include "IO/Instance/Point/TableWriter.h"
#include "IO/Writer/Writer.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/PointMapper.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Kernels/Receiver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Modules/Modules.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <optional>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <utils/logger.h>
#include <vector>

namespace seissol::writer {

ReceiverWriter::ReceiverWriter(seissol::SeisSol& seissolInstance)
    : seissolInstance_(seissolInstance) {}

ReceiverWriter::~ReceiverWriter() = default;

Eigen::Vector3d parseReceiverLine(const std::string& line) {
  const std::regex rgx("\\s+");
  std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
  const std::sregex_token_iterator end;
  Eigen::Vector3d coordinates{};
  std::uint32_t numberOfCoordinates = 0;
  for (; iter != end; ++iter, ++numberOfCoordinates) {
    if (numberOfCoordinates >= coordinates.size()) {
      throw std::runtime_error("Too many coordinates in line " + line + ".");
    }
    coordinates[numberOfCoordinates] = std::stod(*iter);
  }
  if (numberOfCoordinates != coordinates.size()) {
    throw std::runtime_error("To few coordinates in line " + line + ".");
  }
  return coordinates;
}

std::vector<Eigen::Vector3d> parseReceiverFile(const std::string& receiverFileName) {
  std::vector<Eigen::Vector3d> points{};

  std::ifstream file{receiverFileName};
  std::string line{};
  while (std::getline(file, line)) {
    const bool onlyWhiteSpace = std::all_of(line.begin(), line.end(), [](auto& c) {
      return std::isspace(c); // This lambda is needed (opt. argument)
    });
    if (!onlyWhiteSpace) {
      points.emplace_back(parseReceiverLine(line));
    }
  }
  return points;
}

std::string ReceiverWriter::fileName(std::size_t pointId) const {
  std::stringstream fns;
  fns << std::setfill('0') << fileNamePrefix_ << "-receiver-" << std::setw(5) << (pointId + 1);
  fns << ".dat";
  return fns.str();
}

std::vector<std::string> ReceiverWriter::quantityNames() const {
  std::vector<std::string> names(seissol::model::MaterialT::Quantities.begin(),
                                 seissol::model::MaterialT::Quantities.end());
  for (const auto& derived : derivedQuantities_) {
    auto derivedNames = derived->quantities();
    names.insert(names.end(), derivedNames.begin(), derivedNames.end());
  }
  return names;
}

void ReceiverWriter::writeHeader(std::size_t pointId,
                                 const Eigen::Vector3d& point,
                                 std::size_t globalId) {
  auto name = fileName(pointId);

  /// \todo Find a nicer solution that is not so hard-coded.
  struct stat fileStat{};
  // Write header if file does not exist
  if (stat(name.c_str(), &fileStat) != 0) {
    std::ofstream file;
    file.open(name);
    file << "TITLE = \"Temporal Signal for receiver number " << std::setfill('0') << std::setw(5)
         << (pointId + 1) << "\"" << '\n';
    // with fused simulations, a sample takes a row per simulation, which a column of its own names
    file << "VARIABLES = \"Time\"";
    if constexpr (seissol::multisim::MultisimEnabled) {
      file << ",\"SimulationIndex\"";
    }
    for (const auto& name : quantityNames()) {
      file << ",\"" << name << "\"";
    }
    file << '\n';

    // metadata header
    for (int d = 0; d < 3; ++d) {
      file << "# x" << (d + 1) << "       " << std::scientific << std::setprecision(12) << point[d]
           << '\n';
    }
    file << "# cell-global-id " << globalId << '\n';
    file.close();
  }
}

void ReceiverWriter::init(
    const std::string& fileNamePrefix,
    double endTime,
    const seissol::initializer::parameters::ReceiverOutputParameters& parameters) {
  fileNamePrefix_ = fileNamePrefix;
  receiverFileName_ = parameters.fileName;
  samplingInterval_ = parameters.samplingInterval;
  endTime_ = endTime;
  format_ = parameters.format;
  sampleChunk_ = parameters.samplechunk;

  if (parameters.computeRotation) {
    derivedQuantities_.push_back(std::make_shared<kernels::ReceiverRotation>());
  }
  if (parameters.computeStrain) {
    derivedQuantities_.push_back(std::make_shared<kernels::ReceiverStrain>());
  }

  setSyncInterval(std::min(endTime, parameters.writeInterval));
  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  Modules::registerHook(*this, ModuleHook::Shutdown);
}

void ReceiverWriter::addPoints(const seissol::geometry::MeshReader& mesh,
                               const LTS::Backmap& backmap,
                               const CompoundGlobalData& global) {
  std::vector<Eigen::Vector3d> points;
  // Only parse if we have a receiver file
  if (!receiverFileName_.empty()) {
    points = parseReceiverFile(receiverFileName_);
    logInfo() << "Record points read from" << receiverFileName_;
    logInfo() << "Number of record points =" << points.size();
  } else {
    logInfo() << "No record points read.";
  }

  const auto numberOfPoints = points.size();
  std::vector<std::size_t> meshIds(numberOfPoints);

  // We want to plot all quantities except for the memory variables
  std::vector<std::size_t> quantities(seissol::model::MaterialT::Quantities.size());
  std::iota(quantities.begin(), quantities.end(), 0);

  logInfo() << "Finding meshIds for receivers...";
  const auto contained =
      initializer::findUniqueMeshIds(points.data(), mesh, numberOfPoints, meshIds.data(), 1e-3);

  std::vector<short> globalContained(contained.size());
  for (std::size_t i = 0; i < contained.size(); ++i) {
    globalContained[i] = contained[i] ? 1 : 0;
  }

  MPI_Allreduce(MPI_IN_PLACE,
                globalContained.data(),
                globalContained.size(),
                MPI_SHORT,
                MPI_MAX,
                seissol::Mpi::mpi.comm());

  bool receiversMissing = false;
  for (std::size_t i = 0; i < numberOfPoints; ++i) {
    if (globalContained[i] == 0) {
      logWarning() << "Receiver point" << i << "could not be found. Coordinates:" << points[i](0)
                   << points[i](1) << points[i](2);
      receiversMissing = true;
    }
  }
  if (receiversMissing) {
    logError() << "Some receivers could not be found. Aborting simulation.";
  }

  logInfo() << "Mapping receivers to LTS cells...";
  receiverClusters_.clear();

  for (std::size_t point = 0; point < numberOfPoints; ++point) {
    if (contained[point]) {
      const std::size_t meshId = meshIds[point];
      const auto id = backmap.get(meshId).color;

      // Make sure that needed empty clusters are initialized.
      for (std::size_t c = receiverClusters_.size(); c <= id; ++c) {
        receiverClusters_.emplace_back(
            std::make_shared<kernels::ReceiverCluster>(global,
                                                       quantities,
                                                       samplingInterval_,
                                                       syncInterval(),
                                                       derivedQuantities_,
                                                       seissolInstance_));
      }

      if (format_ == seissol::initializer::parameters::ReceiverOutputFormat::Csv) {
        writeHeader(point, points[point], mesh.getElements()[meshId].globalId);
      }

      receiverClusters_[id]->addReceiver(meshId, point, points[point], mesh, backmap);
    }
  }

  if (format_ == seissol::initializer::parameters::ReceiverOutputFormat::Hdf5) {
    // A row of the table is one receiver of one simulation, as on the fault, so what a sample
    // records is the time and the quantities of that simulation. What a receiver records follows
    // from the material of the element it sits in, so the table is told for every row and gathers
    // those that agree into a table of their own.
    std::vector<io::instance::point::TableQuantity> quantitySet;
    quantitySet.push_back(
        io::instance::point::TableQuantity{"Time", io::datatype::inferDatatype<real>()});
    for (const auto& name : quantityNames()) {
      quantitySet.push_back(
          io::instance::point::TableQuantity{name, io::datatype::inferDatatype<real>()});
    }

    const auto receivers = orderedReceivers();
    const std::vector<std::vector<io::instance::point::TableQuantity>> pointQuantities(
        receivers.size() * seissol::multisim::NumSimulations, quantitySet);

    table_ = std::make_unique<io::instance::point::Hdf5Table>(
        "receivers", pointQuantities, seissol::Mpi::mpi.comm(), sampleChunk_);

    // which receiver of the file and which simulation a row belongs to, and where it sits
    std::vector<std::uint64_t> pointIds;
    std::vector<std::uint64_t> simulations;
    std::vector<double> coordinates;
    for (const auto& entry : receivers) {
      for (std::size_t simulation = 0; simulation < seissol::multisim::NumSimulations;
           ++simulation) {
        pointIds.push_back(static_cast<std::uint64_t>(entry.receiver->pointId));
        simulations.push_back(static_cast<std::uint64_t>(simulation));
        for (int dimension = 0; dimension < 3; ++dimension) {
          coordinates.push_back(entry.receiver->position[dimension]);
        }
      }
    }
    table_->addPointData("PointId", {}, pointIds);
    table_->addPointData("SimulationIndex", {}, simulations);
    table_->addPointData("Coordinates", {3}, coordinates);

    io::writer::ScheduledWriter scheduled;
    scheduled.name = "receivers";
    scheduled.interval = syncInterval();
    scheduled.planWrite = [this, plan = table_->makeWriter()](
                              const std::string& prefix, std::size_t counter, double time) {
      stopwatch_.start();
      collectSamples();
      auto writer = plan(prefix, counter, time);
      logInfo() << "Collected receivers in" << stopwatch_.stop() << "seconds.";
      return writer;
    };
    seissolInstance_.outputManager().addOutput(scheduled);
  }
}

std::vector<ReceiverWriter::OrderedReceiver> ReceiverWriter::orderedReceivers() {
  std::vector<OrderedReceiver> receivers;
  for (auto& cluster : receiverClusters_) {
    for (auto& receiver : *cluster) {
      receivers.push_back(OrderedReceiver{&receiver, cluster->ncols()});
    }
  }
  // the rows of a rank are its receivers in the order of the file they were read from, which is
  // the order the point map and the coordinates are written in as well
  std::sort(
      receivers.begin(), receivers.end(), [](const OrderedReceiver& a, const OrderedReceiver& b) {
        return a.receiver->pointId < b.receiver->pointId;
      });
  return receivers;
}

void ReceiverWriter::collectSamples() {
  constexpr auto Simulations = seissol::multisim::NumSimulations;
  const auto receivers = orderedReceivers();
  const auto& grouping = table_->grouping();

  // How far a table grows is the same on every rank, so the sample count is agreed on rather than
  // taken from what this rank happens to hold -- a rank without receivers holds none at all.
  std::vector<std::size_t> samples(grouping.groupCount(), 0);
  for (std::size_t i = 0; i < receivers.size(); ++i) {
    const auto held = receivers[i].receiver->output.size() / receivers[i].columns;
    for (std::size_t simulation = 0; simulation < Simulations; ++simulation) {
      const auto group = grouping.group[i * Simulations + simulation];
      samples[group] = std::max(samples[group], held);
    }
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
    storage[group] = reinterpret_cast<real*>(table_->prepare(group, samples[group]));
    // A receiver with fewer samples than the longest one of its table leaves the rest of its
    // column unset, and a zero there is a value a reader cannot tell from a measurement.
    const auto values =
        samples[group] * table_->localPointCount(group) * table_->sampleSize(group) / sizeof(real);
    std::fill_n(storage[group], values, std::numeric_limits<real>::quiet_NaN());
  }

  for (std::size_t i = 0; i < receivers.size(); ++i) {
    auto& receiver = *receivers[i].receiver;
    const auto columns = receivers[i].columns;
    const auto held = receiver.output.size() / columns;
    // a sample holds the time, then the quantities of one simulation after the other
    const auto quantities = (columns - 1) / Simulations;

    for (std::size_t simulation = 0; simulation < Simulations; ++simulation) {
      const auto point = i * Simulations + simulation;
      const auto group = grouping.group[point];
      const auto row = table_->localRow(point);
      const auto points = table_->localPointCount(group);
      for (std::size_t sample = 0; sample < std::min(held, samples[group]); ++sample) {
        const auto* source = receiver.output.data() + sample * columns;
        auto* target = storage[group] + (sample * points + row) * (1 + quantities);
        target[0] = source[0];
        std::copy_n(source + 1 + simulation * quantities, quantities, target + 1);
      }
    }
    receiver.output.clear();
  }
}

// --------------------------------------------------------------------------
void ReceiverWriter::syncPoint(double /*currentTime*/) {
  // the HDF5 table is filled by the scheduled writer registered in addPoints
  if (format_ != seissol::initializer::parameters::ReceiverOutputFormat::Csv ||
      receiverClusters_.empty()) {
    return;
  }

  stopwatch_.start();

  constexpr auto Simulations = seissol::multisim::NumSimulations;
  for (auto& cluster : receiverClusters_) {
    const auto ncols = cluster->ncols();
    // a sample holds the time, then the quantities of one simulation after the other
    const auto quantities = (ncols - 1) / Simulations;
    for (auto& receiver : *cluster) {
      assert(receiver.output.size() % ncols == 0);
      const std::size_t nSamples = receiver.output.size() / ncols;

      std::ofstream file;
      file.open(fileName(receiver.pointId), std::ios::app);
      file << std::scientific << std::setprecision(15);
      for (std::size_t i = 0; i < nSamples; ++i) {
        const auto* sample = receiver.output.data() + i * ncols;
        // with fused simulations, a row per simulation, which the second column names
        for (std::size_t simulation = 0; simulation < Simulations; ++simulation) {
          file << "  " << sample[0];
          if constexpr (seissol::multisim::MultisimEnabled) {
            file << "  " << simulation;
          }
          for (std::size_t q = 0; q < quantities; ++q) {
            file << "  " << sample[1 + simulation * quantities + q];
          }
          file << '\n';
        }
      }
      file.close();
      receiver.output.clear();
    }
  }

  const auto time = stopwatch_.stop();
  logInfo() << "Wrote receivers in" << time << "seconds.";
}

// --------------------------------------------------------------------------
void ReceiverWriter::simulationStart(std::optional<double> /*checkpointTime*/) {
  for (auto& cluster : receiverClusters_) {
    cluster->allocateData();
  }
}

// --------------------------------------------------------------------------
void ReceiverWriter::shutdown() {
  for (auto& cluster : receiverClusters_) {
    cluster->freeData();
  }
  table_.reset();
}

// --------------------------------------------------------------------------
kernels::ReceiverCluster* ReceiverWriter::receiverCluster(std::size_t id) {
  if (id < receiverClusters_.size()) {
    return receiverClusters_[id].get();
  }
  return nullptr;
}

} // namespace seissol::writer
