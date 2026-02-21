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
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::instance::mesh {
VtkHdfWriter::VtkHdfWriter(const std::string& name,
                           std::size_t localElementCount,
                           std::size_t dimension,
                           std::size_t targetDegree)
    : name_(name), localElementCount_(localElementCount), globalElementCount_(localElementCount),
      pointsPerElement_(dimension == 2
                            ? ((targetDegree + 1) * (targetDegree + 2)) / 2
                            : ((targetDegree + 1) * (targetDegree + 2) * (targetDegree + 3)) / 6),
      type_(dimension == 2 ? 69 : 71), targetDegree_(targetDegree) {
  // 69: Lagrange triangle
  // 71: Lagrange tetrahedron

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

  instructions_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Type",
        writer::WriteInline::create("UnstructuredGrid",
                                    std::make_shared<datatype::StringDatatype>(16)));
  });
  instructions_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Version",
        writer::WriteInline::createArray<int64_t>({2}, {1, 0}));
  });

  // to capture by value
  auto selfGlobalElementCount = globalElementCount_;
  auto selfLocalElementCount = localElementCount;
  auto selfGlobalPointCount = globalPointCount_;
  auto selfLocalPointCount = localPointCount_;
  auto selfPointOffset = pointOffset_;
  auto selfPointsPerElement = pointsPerElement_;
  auto selfType = type_;

  // TODO: move the following arrays into a "common" HDF5 file
  // also, auto-generate them using a managed buffer
  instructionsConst_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "NumberOfCells",
        writer::WriteInline::createArray<int64_t>({1},
                                                  {static_cast<int64_t>(selfGlobalElementCount)}),
        datatype::inferDatatype<int64_t>());
  });
  instructionsConst_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "NumberOfConnectivityIds",
        writer::WriteInline::createArray<int64_t>({1},
                                                  {static_cast<int64_t>(selfGlobalPointCount)}),
        datatype::inferDatatype<int64_t>());
  });
  instructionsConst_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "NumberOfPoints",
        writer::WriteInline::createArray<int64_t>({1},
                                                  {static_cast<int64_t>(selfGlobalPointCount)}),
        datatype::inferDatatype<int64_t>());
  });

  const bool isLastRank = Mpi::mpi.size() == Mpi::mpi.rank() + 1;
  instructionsConst_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Offsets",
        writer::GeneratedBuffer::createElementwise<int64_t>(
            selfLocalElementCount + (isLastRank ? 1 : 0),
            1,
            std::vector<std::size_t>(),
            [=](int64_t* target, std::size_t index) {
              target[0] = index * selfPointsPerElement + selfPointOffset;
            }),
        datatype::inferDatatype<int64_t>());
  });
  instructionsConst_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Types",
        writer::GeneratedBuffer::createElementwise<uint8_t>(
            selfLocalElementCount,
            1,
            std::vector<std::size_t>(),
            [=](uint8_t* target, std::size_t /*index*/) { target[0] = selfType; }),
        datatype::inferDatatype<uint8_t>());
  });
  instructionsConst_.emplace_back([=](const std::string& filename, double /*time*/) {
    return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Connectivity",
        writer::GeneratedBuffer::createElementwise<int64_t>(
            selfLocalPointCount,
            1,
            std::vector<std::size_t>(),
            [=](int64_t* target, std::size_t index) { target[0] = index + selfPointOffset; }),
        datatype::inferDatatype<int64_t>());
  });
}

void VtkHdfWriter::addHook(const std::function<void(std::size_t, double)>& hook) {
  hooks_.push_back(hook);
}

std::function<writer::Writer(const std::string&, std::size_t, double)> VtkHdfWriter::makeWriter() {
  logInfo() << "Adding VTK writer" << name_ << "of order" << targetDegree_;
  auto self = *this;
  return [self](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
    for (const auto& hook : self.hooks_) {
      hook(counter, time);
    }
    const auto filename = prefix + "-" + self.name_ + "-" + std::to_string(counter) + ".vtkhdf";
    auto writer = writer::Writer();
    for (const auto& instruction : self.instructionsConst_) {
      writer.addInstruction(instruction(filename, time));
    }
    for (const auto& instruction : self.instructions_) {
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
    return writer;
  };
}
} // namespace seissol::io::instance::mesh
