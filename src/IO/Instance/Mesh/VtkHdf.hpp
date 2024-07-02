#pragma once

#include "utils/logger.h"
#include <IO/Datatype/Datatype.hpp>
#include <IO/Datatype/Inference.hpp>
#include <IO/Datatype/MPIType.hpp>
#include <IO/Instance/SeisSolMemoryHelper.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <IO/Writer/Instructions/Instruction.hpp>
#include <IO/Writer/Writer.hpp>
#include <Initializer/MemoryManager.h>
#include <functional>
#include <memory>

namespace seissol::io::instance::mesh {
class VtkHdfWriter {
  public:
  VtkHdfWriter(const std::string& name,
               std::size_t localElementCount,
               std::size_t dimension,
               std::size_t targetDegree)
      : localElementCount(localElementCount), globalElementCount(localElementCount),
        elementOffset(0), name(name) {
    // 69: Lagrange triangle
    // 71: Lagrange tetrahedron

    type = dimension == 2 ? 69 : 71;
    pointsPerElement = dimension == 2
                           ? ((targetDegree + 1) * (targetDegree + 2)) / 2
                           : ((targetDegree + 1) * (targetDegree + 2) * (targetDegree + 3)) / 6;

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

    instructions.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "Type",
          writer::WriteInline::create("UnstructuredGrid",
                                      std::make_shared<datatype::StringDatatype>(16)));
    });
    instructions.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "Version",
          writer::WriteInline::createArray<int64_t>({2}, {1, 0}));
    });

    // to capture by value
    auto selfGlobalElementCount = globalElementCount;
    auto selfLocalElementCount = localElementCount;
    auto selfElementOffset = elementOffset;
    auto selfGlobalPointCount = globalPointCount;
    auto selfLocalPointCount = localPointCount;
    auto selfPointOffset = pointOffset;
    auto selfPointsPerElement = pointsPerElement;
    auto selfType = type;

    // TODO: move the following arrays into a "common" HDF5 file
    // also, auto-generate them using a managed buffer
    instructionsConst.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "NumberOfCells",
          writer::WriteInline::createArray<int64_t>({1},
                                                    {static_cast<int64_t>(selfGlobalElementCount)}),
          datatype::inferDatatype<int64_t>());
    });
    instructionsConst.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "NumberOfConnectivityIds",
          writer::WriteInline::createArray<int64_t>({1},
                                                    {static_cast<int64_t>(selfGlobalPointCount)}),
          datatype::inferDatatype<int64_t>());
    });
    instructionsConst.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "NumberOfPoints",
          writer::WriteInline::createArray<int64_t>({1},
                                                    {static_cast<int64_t>(selfGlobalPointCount)}),
          datatype::inferDatatype<int64_t>());
    });

    bool isLastRank = MPI::mpi.size() == MPI::mpi.rank() + 1;
    instructionsConst.emplace_back([=](const std::string& filename, double time) {
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
    instructionsConst.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "Types",
          writer::GeneratedBuffer::createElementwise<uint8_t>(
              selfLocalElementCount,
              1,
              std::vector<std::size_t>(),
              [=](uint8_t* target, std::size_t index) { target[0] = selfType; }),
          datatype::inferDatatype<uint8_t>());
    });
    instructionsConst.emplace_back([=](const std::string& filename, double time) {
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

  template <typename F>
  void addPointProjector(F&& projector) {
    auto selfLocalElementCount = localElementCount;
    auto selfPointsPerElement = pointsPerElement;

    instructionsConst.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "Points",
          writer::GeneratedBuffer::createElementwise<double>(
              selfLocalElementCount, selfPointsPerElement, std::vector<std::size_t>{3}, projector),
          datatype::inferDatatype<double>());
    });
  }

  template <typename T, typename F>
  void addPointData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    F&& pointMapper) {
    auto selfLocalElementCount = localElementCount;
    auto selfPointsPerElement = pointsPerElement;

    instructions.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName, PointDataName}),
          name,
          writer::GeneratedBuffer::createElementwise<T>(
              selfLocalElementCount, selfPointsPerElement, dimensions, pointMapper),
          datatype::inferDatatype<T>());
    });
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   F&& cellMapper) {
    auto selfLocalElementCount = localElementCount;

    instructions.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName, CellDataName}),
          name,
          writer::GeneratedBuffer::createElementwise<T>(
              selfLocalElementCount, 1, dimensions, cellMapper),
          datatype::inferDatatype<T>());
    });
  }

  template <typename T>
  void addFieldData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    const std::vector<T>& data) {
    instructions.emplace_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
          name,
          writer::WriteInline::createArray(dimensions, data),
          datatype::inferDatatype<T>());
    });
  }

  void addHook(const std::function<void(std::size_t, double)>& hook) { hooks.push_back(hook); }

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() {
    auto self = *this;
    return [self](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
      for (const auto& hook : self.hooks) {
        hook(counter, time);
      }
      const auto filename = prefix + "-" + self.name + "-" + std::to_string(counter) + ".vtkhdf";
      auto writer = writer::Writer();
      for (auto& instruction : self.instructionsConst) {
        writer.addInstruction(instruction(filename, time));
      }
      for (auto& instruction : self.instructions) {
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

  private:
  std::string name;
  std::size_t localElementCount;
  std::size_t globalElementCount;
  std::size_t elementOffset;
  std::size_t localPointCount;
  std::size_t globalPointCount;
  std::size_t pointOffset;
  std::size_t pointsPerElement;
  std::vector<std::function<void(std::size_t, double)>> hooks;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, double)>>
      instructionsConst;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, double)>>
      instructions;
  std::size_t type;
  const static inline std::string GroupName = "VTKHDF";
  const static inline std::string FieldDataName = "FieldData";
  const static inline std::string CellDataName = "CellData";
  const static inline std::string PointDataName = "PointData";
};
} // namespace seissol::io::instance::mesh
