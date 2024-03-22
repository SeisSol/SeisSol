#pragma once

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
  VtkHdfWriter(std::size_t localElementCount, std::size_t dimension, std::size_t targetDegree)
    : localElementCount(localElementCount), globalElementCount(localElementCount), elementOffset(0) {
    // 69: Lagrange triangle
    // 71: Lagrange tetrahedron

    type = dimension == 2 ? 69 : 71;
    pointsPerElement =
        dimension == 2 ? ((targetDegree + 1) * (targetDegree + 2)) / 2
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

    instructions.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "Type",
          writer::WriteInline::create("UnstructuredGrid"));
    });
    instructions.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
          writer::instructions::Hdf5Location(filename, {GroupName}),
          "Version",
          writer::WriteInline::createArray<uint64_t>({2}, {1, 0}));
    });

    // TODO: move the following arrays into a "common" HDF5 file
    // also, auto-generate them using a managed buffer
    instructionsConst.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "NumberOfCells",
        writer::WriteInline::createArray<uint64_t>({1}, {globalElementCount}),
        datatype::inferDatatype<uint64_t>());
    });
    instructionsConst.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "NumberOfConnectivityIds",
        writer::WriteInline::createArray<uint64_t>({1}, {globalPointCount}),
        datatype::inferDatatype<uint64_t>());
    });
    instructionsConst.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "NumberOfPoints",
        writer::WriteInline::createArray<uint64_t>({1}, {globalPointCount}),
        datatype::inferDatatype<uint64_t>());
    });
    instructionsConst.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Offsets",
        std::make_shared<writer::ElementBuffer<uint64_t, void(uint64_t* target, std::size_t index)>>(
            localElementCount, 1, std::vector<std::size_t>(), [&](uint64_t* target, std::size_t index) {
              target[0] = index * pointsPerElement + pointOffset;
            }));
    });
    /*writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Points",
        std::make_shared<writer::ElementBuffer<double>>(cellCount,
                                                        pointsPerCell,
                                                        std::vector<std::size_t>{dimension},
                                                        std::forward<F>(pointIndexer))));*/
    instructionsConst.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Types",
        std::make_shared<writer::ElementBuffer<uint64_t, void(uint64_t* target, std::size_t index)>>(
            localElementCount, 1, std::vector<std::size_t>(), [&](uint64_t* target, std::size_t index) {
              target[0] = type;
            }));
    });
    instructionsConst.push_back([=](const std::string& filename, double time) {
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {GroupName}),
        "Connectivity",
        std::make_shared<writer::ElementBuffer<uint64_t, void(uint64_t* target, std::size_t index)>>(
            globalPointCount, 1, std::vector<std::size_t>(), [&](uint64_t* target, std::size_t index) {
              target[0] = index + elementOffset;
            }));
    });
  }

  template <typename F>
  void addPointProjector(F&& projector) {}

  template <typename T, typename F>
  void addPointData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    F&& pointMapper) {
    instructions.push_back([=](const std::string& filename, double time) {
      return writer::instructions::Hdf5DataWrite(
          writer::instructions::Hdf5Location(filename, {GroupName, PointDataName}),
          name,
          writer::ElementBuffer<T, F>(localElementCount, pointsPerElement, dimensions, std::forward<F>(pointMapper)));
    });
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   F&& cellMapper) {
    instructions.push_back([=](const std::string& filename, double time) {
      return writer::instructions::Hdf5DataWrite(
          writer::instructions::Hdf5Location(filename, {GroupName, CellDataName}),
          name,
          writer::ElementBuffer<T, F>(localElementCount, 1, dimensions, std::forward<F>(cellMapper)));
    });
  }

  template <typename T>
  void addFieldData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    const std::vector<T>& data) {
    instructions.push_back([=](const std::string& filename, double time) {
      return writer::instructions::Hdf5DataWrite(
          writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
          name,
          writer::WriteInline::createArray(dimensions, data));
    });
  }

  std::function<writer::Writer(double)> makeWriter() {
    auto self = *this;
    return [self](double time) -> writer::Writer {
      std::string filename = "";
      auto writer = writer::Writer();
      for (auto& instruction : self.instructionsConst) {
        writer.addInstruction(instruction(filename, time));
      }
      for (auto& instruction : self.instructions) {
        writer.addInstruction(instruction(filename, time));
      }
      return writer;
    };
  }

  private:
  std::size_t localElementCount;
  std::size_t globalElementCount;
  std::size_t elementOffset;
  std::size_t localPointCount;
  std::size_t globalPointCount;
  std::size_t pointOffset;
  std::size_t pointsPerElement;
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
