// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "VtkHdf.h"

#include <IO/Datatype/Datatype.h>
#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <IO/Instance/Geometry/Typedefs.h>
#include <IO/Instance/Metadata/Pvd.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <IO/Writer/Writer.h>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

#include "utils/logger.h"

namespace seissol::io::instance::mesh {
VtkHdfWriter::VtkHdfWriter(const std::string& name,
                           std::size_t localElementCount,
                           geometry::Shape shape,
                           std::size_t targetDegree)
    : name(name), localElementCount(localElementCount), globalElementCount(localElementCount),
      pointsPerElement(geometry::numPoints(targetDegree, shape)), type(geometry::vtkType(shape)),
      targetDegree(targetDegree) {
  MPI_Exscan(&localElementCount,
             &elementOffset,
             1,
             datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
             MPI_SUM,
             seissol::MPI::mpi.comm());
  MPI_Allreduce(&localElementCount,
                &globalElementCount,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::MPI::mpi.comm());
  pointOffset = elementOffset * pointsPerElement;
  localPointCount = localElementCount * pointsPerElement;
  globalPointCount = globalElementCount * pointsPerElement;

  addData("Type",
          {},
          false,
          writer::WriteInline::create("UnstructuredGrid",
                                      std::make_shared<datatype::StringDatatype>(16)));
  addData("Version", {}, false, writer::WriteInline::createArray<int64_t>({2}, {1, 0}));

  // to capture by value
  auto selfGlobalElementCount = globalElementCount;
  auto selfLocalElementCount = localElementCount;
  auto selfGlobalPointCount = globalPointCount;
  auto selfLocalPointCount = localPointCount;
  auto selfPointOffset = pointOffset;
  auto selfPointsPerElement = pointsPerElement;
  auto selfType = type;

  // TODO: move the following arrays into a "common" HDF5 file
  // also, auto-generate them using a managed buffer
  addData("NumberOfCells",
          {},
          true,
          writer::WriteInline::createArray<int64_t>(
              {1}, {static_cast<int64_t>(selfGlobalElementCount)}));
  addData(
      "NumberOfConnectivityIds",
      {},
      true,
      writer::WriteInline::createArray<int64_t>({1}, {static_cast<int64_t>(selfGlobalPointCount)}));
  addData(
      "NumberOfPoints",
      {},
      true,
      writer::WriteInline::createArray<int64_t>({1}, {static_cast<int64_t>(selfGlobalPointCount)}));

  const bool isLastRank = MPI::mpi.size() == MPI::mpi.rank() + 1;
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
              [=](uint8_t* target, std::size_t index) { target[0] = selfType; }));
  addData("Connectivity",
          {},
          true,
          writer::GeneratedBuffer::createElementwise<int64_t>(
              selfLocalPointCount,
              1,
              std::vector<std::size_t>(),
              [=](int64_t* target, std::size_t index) { target[0] = index + selfPointOffset; }));
}

void VtkHdfWriter::addData(const std::string& name,
                           const std::optional<std::string>& group,
                           bool isConst,
                           const std::shared_ptr<writer::DataSource>& data) {
  auto& instrarray = isConst ? instructionsConst : instructions;

  std::vector<std::string> groups{GroupName};
  if (group.has_value()) {
    groups.emplace_back(group.value());
  }

  instrarray.emplace_back([=](const std::string& filename, double time) {
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, groups), name, data, data->datatype());
  });

  if (isConst) {
    instructionsConstLink.emplace_back(
        [=](const std::string& filename, const std::string& filenameConst) {
          return std::make_shared<writer::instructions::Hdf5LinkExternalWrite>(
              writer::instructions::Hdf5Location(filename, groups),
              name,
              writer::instructions::Hdf5Location(filenameConst, groups));
        });
  }
}

void VtkHdfWriter::addHook(const std::function<void(std::size_t, double)>& hook) {
  hooks.push_back(hook);
}

std::function<writer::Writer(const std::string&, std::size_t, double)> VtkHdfWriter::makeWriter() {
  logInfo() << "Adding VTK writer" << name << "of order" << targetDegree;
  const auto self = *this;
  return [self, pvu = std::vector<metadata::PvuEntry>()](const std::string& prefix,
                                                         std::size_t counter,
                                                         double time) mutable -> writer::Writer {
    for (const auto& hook : self.hooks) {
      hook(counter, time);
    }
    const auto filename = prefix + "-" + self.name + "-" + std::to_string(counter) + ".vtkhdf";
    const auto filenameConst = prefix + "-" + self.name + "-const.vtkhdf";
    const auto filenamePvu = prefix + "-" + self.name + ".pvu";
    pvu.emplace_back(metadata::PvuEntry{filename, time});
    auto writer = writer::Writer();

    const auto fullWrite = counter == 0;
    if (fullWrite) {
      for (const auto& instruction : self.instructionsConst) {
        writer.addInstruction(instruction(filenameConst, time));
      }
    }
    for (const auto& instruction : self.instructionsConstLink) {
      writer.addInstruction(instruction(filename, filenameConst));
    }
    for (const auto& instruction : self.instructions) {
      writer.addInstruction(instruction(filename, time));
    }
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
        "Time",
        writer::WriteInline::createArray<double>({1}, {time}),
        datatype::inferDatatype<decltype(time)>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
        "Index",
        writer::WriteInline::createArray<std::size_t>({1}, {counter}),
        datatype::inferDatatype<decltype(counter)>()));
    writer.addInstructions(metadata::makePvu(pvu).instructions(filenamePvu));
    return writer;
  };
}
} // namespace seissol::io::instance::mesh
