// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "VtkHdf.h"

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Datatype/MPIType.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Metadata/Pvd.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <optional>
#include <string>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::io::instance::mesh {
VtkHdfWriter::VtkHdfWriter(const std::string& name,
                           std::size_t localElementCount,
                           geometry::Shape shape,
                           std::size_t targetDegree,
                           bool temporal,
                           std::int32_t compress,
                           bool constFile)
    : name_(name), localElementCount_(localElementCount), globalElementCount_(localElementCount),
      pointsPerElement_(
          geometry::numPoints(std::max(targetDegree, static_cast<std::size_t>(1)), shape)),
      type_(geometry::vtkType(shape)), targetDegree_(targetDegree), temporal_(temporal),
      compress_(compress), constFile_(constFile) {
  MPI_Exscan(&localElementCount,
             &elementOffset_,
             1,
             datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
             MPI_SUM,
             seissol::Mpi::mpi.comm());
  MPI_Allreduce(&localElementCount,
                &globalElementCount_,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::Mpi::mpi.comm());
  pointOffset_ = elementOffset_ * pointsPerElement_;
  localPointCount_ = localElementCount * pointsPerElement_;
  globalPointCount_ = globalElementCount_ * pointsPerElement_;

  const auto version = temporal ? std::vector<int64_t>{2, 0} : std::vector<int64_t>{1, 0};

  addData("Type",
          {},
          temporal,
          writer::WriteInline::create("UnstructuredGrid",
                                      std::make_shared<datatype::StringDatatype>(16)),
          true);
  addData("Version",
          {},
          temporal,
          writer::WriteInline::createArray<int64_t>({version.size()}, version),
          true);

  // to capture by value
  const auto selfGlobalElementCount = globalElementCount_;
  const auto selfLocalElementCount = localElementCount_;
  const auto selfGlobalPointCount = globalPointCount_;
  const auto selfLocalPointCount = localPointCount_;
  const auto selfPointOffset = pointOffset_;
  const auto selfPointsPerElement = pointsPerElement_;
  const auto selfType = type_;

  // TODO: auto-generate using a managed buffer maybe?

  addData("NumberOfCells",
          {},
          temporal,
          writer::WriteInline::createArray<int64_t>(
              {1}, {static_cast<int64_t>(selfGlobalElementCount)}));
  addData(
      "NumberOfConnectivityIds",
      {},
      temporal,
      writer::WriteInline::createArray<int64_t>({1}, {static_cast<int64_t>(selfGlobalPointCount)}));
  addData(
      "NumberOfPoints",
      {},
      temporal,
      writer::WriteInline::createArray<int64_t>({1}, {static_cast<int64_t>(selfGlobalPointCount)}));

  const bool isLastRank = Mpi::mpi.size() == Mpi::mpi.rank() + 1;
  addData("Offsets",
          {},
          true,
          writer::GeneratedBuffer::createElementwise<int64_t>(
              selfLocalElementCount + (isLastRank ? 1 : 0),
              1,
              std::vector<std::size_t>(),
              [=](int64_t* target, std::size_t index) {
                target[0] = index * selfPointsPerElement + selfPointOffset;
              }));
  addData("Types",
          {},
          true,
          writer::GeneratedBuffer::createElementwise<uint8_t>(
              selfLocalElementCount,
              1,
              std::vector<std::size_t>(),
              [=](uint8_t* target, std::size_t /*index*/) { target[0] = selfType; }));
  addData("Connectivity",
          {},
          true,
          writer::GeneratedBuffer::createElementwise<int64_t>(
              selfLocalPointCount,
              1,
              std::vector<std::size_t>(),
              [=](int64_t* target, std::size_t index) { target[0] = index + selfPointOffset; }));

  if (temporal) {
    // https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html
    // The mesh itself is written once and every step reads it again, so all of the geometry
    // offsets stay at zero ("Offset value can be repeated for static data"); only the attribute
    // data grows.
    instructions_.emplace_back(
        [](const std::string& filename, std::size_t counter, double /*time*/) {
          return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
              writer::instructions::Hdf5Location(filename, {GroupName, StepsName}),
              "NSteps",
              writer::WriteInline::createArray<int64_t>({}, {static_cast<int64_t>(counter) + 1}));
        });

    instructions_.emplace_back(
        [](const std::string& filename, std::size_t /*counter*/, double time) {
          const auto data = writer::WriteInline::createArray<double>({}, {time});
          return std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, {GroupName, StepsName}),
              "Values",
              data,
              data->datatype(),
              true);
        });

    // NumberOfPoints and friends hold one entry for the single part the mesh was written as, not
    // one per step, so the part count has to be given explicitly rather than inferred.
    instructions_.emplace_back(
        [](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
          const auto data = writer::WriteInline::createArray<int64_t>({}, {1});
          return std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, {GroupName, StepsName}),
              "NumberOfParts",
              data,
              data->datatype(),
              true);
        });

    addStepOffset("PartOffsets", {}, 0);
    addStepOffset("PointOffsets", {}, 0);
    // CellOffsets and ConnectivityIdOffsets are (NSteps, NTopologies), and an unstructured grid
    // has one topology
    addStepOffset("CellOffsets", {}, 0);
    addStepOffset("ConnectivityIdOffsets", {}, 0);

    // the two field data entries makeWriter always adds
    addStepOffset("Time", FieldDataName + "Offsets", 1);
    addStepOffset("Index", FieldDataName + "Offsets", 1);
  }
}

void VtkHdfWriter::addStepOffset(const std::string& name,
                                 const std::optional<std::string>& group,
                                 std::size_t perStep) {
  std::vector<std::string> groups{GroupName, StepsName};
  if (group.has_value()) {
    groups.emplace_back(group.value());
  }
  const std::vector<std::size_t> shape = (name == "CellOffsets" || name == "ConnectivityIdOffsets")
                                             ? std::vector<std::size_t>{1}
                                             : std::vector<std::size_t>{};

  instructions_.emplace_back([=](const std::string& filename,
                                 std::size_t counter,
                                 double /*time*/) {
    const auto data = writer::WriteInline::createArray<uint64_t>(
        shape, {static_cast<uint64_t>(counter * perStep)});
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, groups), name, data, data->datatype(), true);
  });
}

void VtkHdfWriter::addData(const std::string& name,
                           const std::optional<std::string>& group,
                           bool isConst,
                           const std::shared_ptr<writer::DataSource>& data,
                           bool attribute) {
  auto& instrarray = isConst && (constFile_ || temporal_) ? instructionsConst_ : instructions_;

  std::vector<std::string> groups{GroupName};
  if (group.has_value()) {
    groups.emplace_back(group.value());
  }

  const auto append = !isConst && temporal_;

  const auto compress = this->compress_;

  if (attribute) {
    instrarray.emplace_back(
        [=](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
          return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
              writer::instructions::Hdf5Location(filename, groups), name, data);
        });
  } else {
    instrarray.emplace_back(
        [=](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
          return std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, groups),
              name,
              data,
              data->datatype(),
              append,
              compress);
        });

    if (isConst && constFile_ && !temporal_) {
      instructionsConstLink_.emplace_back(
          [=](const std::string& filename, const std::string& filenameConst) {
            return std::make_shared<writer::instructions::Hdf5LinkExternalWrite>(
                writer::instructions::Hdf5Location(filename, groups),
                name,
                writer::instructions::Hdf5Location(filenameConst, groups, name));
          });
    }
  }
}

void VtkHdfWriter::addHook(const std::function<void(std::size_t, double)>& hook) {
  hooks_.push_back(hook);
}

std::function<writer::Writer(const std::string&, std::size_t, double)> VtkHdfWriter::makeWriter() {
  logInfo() << "Adding VTK writer" << name_ << "of order" << targetDegree_;
  const auto self = *this;
  return [self, pvu = std::vector<metadata::PvuEntry>()](const std::string& prefix,
                                                         std::size_t counter,
                                                         double time) mutable -> writer::Writer {
    for (const auto& hook : self.hooks_) {
      hook(counter, time);
    }

    const auto lastPrefix = utils::StringUtils::split(prefix, '/');
    // a time series is one file holding every step; a snapshot is one file per step
    const auto suffix = self.temporal_ ? std::string() : "-" + std::to_string(counter);
    const auto filename = prefix + "-" + self.name_ + suffix + ".vtkhdf";
    const auto filenameFile = lastPrefix.back() + "-" + self.name_ + suffix + ".vtkhdf";
    const auto filenameConst = prefix + "-" + self.name_ + "-const.vtkhdf";
    const auto filenameConstFile = lastPrefix.back() + "-" + self.name_ + "-const.vtkhdf";
    const auto filenamePvu = prefix + "-" + self.name_ + ".pvd";
    pvu.emplace_back(metadata::PvuEntry{filenameFile, time});
    auto writer = writer::Writer();

    const auto fullWrite = counter == 0;
    if (fullWrite) {
      // with a const file the unchanging data lives next to the snapshots; in a time series it
      // goes into the one file, once
      const auto& constTarget = self.constFile_ ? filenameConst : filename;
      for (const auto& instruction : self.instructionsConst_) {
        writer.addInstruction(instruction(constTarget, counter, time));
      }
    }
    for (const auto& instruction : self.instructionsConstLink_) {
      writer.addInstruction(instruction(filename, filenameConstFile));
    }
    for (const auto& instruction : self.instructions_) {
      writer.addInstruction(instruction(filename, counter, time));
    }
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
        "Time",
        writer::WriteInline::createArray<double>({1}, {time}),
        datatype::inferDatatype<decltype(time)>(),
        self.temporal_));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
        "Index",
        writer::WriteInline::createArray<std::size_t>({1}, {counter}),
        datatype::inferDatatype<decltype(counter)>(),
        self.temporal_));
    writer.addInstructions(metadata::makePvu(pvu).instructions(filenamePvu));
    return writer;
  };
}
} // namespace seissol::io::instance::mesh
