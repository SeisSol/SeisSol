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
#include <type_traits>
namespace seissol::io::instance::mesh {
    class VtkHdfWriter {
    public:
        VtkHdfWriter(std::size_t dimension, std::size_t targetDegree) {
            const std::size_t type = dimension == 2 ? 69 : 71;
            const std::size_t pointsPerCell = dimension == 2 ? ((targetDegree + 1) * (targetDegree + 2)) / 2 : ((targetDegree + 1) * (targetDegree + 2) * (targetDegree + 3)) / 6;

            const std::size_t cellCount = meshReader.getElements().size();
            std::size_t totalCellCount = meshReader.getElements().size();
            std::size_t cellOffset = 0;
            MPI_Exscan( &cellCount,&cellOffset,
                        1,
                        datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                        MPI_SUM,
                        seissol::MPI::mpi.comm());
            MPI_Allreduce(&cellCount, &totalCellCount,
                            1,
                            datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                            MPI_SUM,
                            seissol::MPI::mpi.comm());
            const std::size_t pointOffset = cellOffset * pointsPerCell;
            const std::size_t pointCount = cellCount * pointsPerCell;
            const std::size_t totalPointCount = totalCellCount * pointsPerCell;
        }

        template<typename F>
        void addPointProjector(F&& projector) {

        }

        template<typename G>
        void addPointData(const std::string& name, const std::vector<std::size_t>& dimensions, G&& pointMapper) {

        }

        template<typename G>
        void addCellData(const std::string& name, const std::vector<std::size_t>& dimensions, G&& cellMapper) {

        }

        template<typename T>
        void addFieldData(const std::string& name, const std::vector<std::size_t>& dimensions, const std::vector<T>& data) {

        }

        std::function<writer::Writer(double)> makeWriter() {
            auto self = *this;
            return [self](double time) -> writer::Writer {
                auto writer = writer::Writer();
                if (time == 0) {
                    for (auto& instruction : self.instructionsConst) {
                        writer.addInstruction(instruction);
                    }
                }
                for (auto& instruction : self.instructions) {
                    writer.addInstruction(instruction);
                }
                return writer;
            };
        }
    private:
        std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(const std::string&, double)>> instructionsConst;
        std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(const std::string&, double)>> instructions;
        const std::size_t type = 0;
    };
template <typename F, typename G>
std::function<writer::Writer(double)> makeVTKHDFWriter(geometry::MeshReader& meshReader,
                                                       initializer::MemoryManager& memoryManager,
                                                       std::size_t pointsPerCell,
                                                       std::size_t dimension,
                                                       F&& pointIndexer,
                                                       G&& projector) {

// 69: Lagrange triangle
// 71: Lagrange tetrahedron
  const std::size_t type = dimension == 2 ? 69 : 71;

  const std::size_t cellCount = meshReader.getElements().size();
  std::size_t totalCellCount = meshReader.getElements().size();
  std::size_t cellOffset = 0;
  MPI_Exscan(&cellCount,
             &cellOffset,
             1,
             datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
             MPI_SUM,
             seissol::MPI::mpi.comm());
  MPI_Allreduce(&cellCount,
                &totalCellCount,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::MPI::mpi.comm());
  const std::size_t pointOffset = cellOffset * pointsPerCell;
  const std::size_t pointCount = cellCount * pointsPerCell;
  const std::size_t totalPointCount = totalCellCount * pointsPerCell;

  return [&](double time) -> writer::Writer {
    writer::Writer writer;
    auto mainGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF"});
    auto cellDataGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF", "CellData"});
    auto pointDataGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF", "PointData"});
    auto fieldDataGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF", "FieldData"});

    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        mainGroup, "Type", writer::WriteInline::create("UnstructuredGrid")));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        mainGroup, "Version", writer::WriteInline::createArray<uint64_t>({2}, {1, 0})));
    
    // TODO: move the following arrays into a "common" HDF5 file
    // also, auto-generate them using a managed buffer
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "NumberOfCells",
        writer::WriteInline::createArray<uint64_t>({1}, {totalCellCount}),
        datatype::inferDatatype<uint64_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "NumberOfConnectivityIds",
        writer::WriteInline::createArray<uint64_t>({1}, {totalPointCount}),
        datatype::inferDatatype<uint64_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "NumberOfPoints",
        writer::WriteInline::createArray<uint64_t>({1}, {totalPointCount}),
        datatype::inferDatatype<uint64_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Offsets",
        std::make_shared<writer::ElementBuffer<uint64_t>>(cellCount,
                                        1,
                                        std::vector<std::size_t>(),
                                        [&](uint64_t* target, std::size_t index) {
                                          target[0] = index * pointsPerCell + pointOffset;
                                        }),
        datatype::inferDatatype<uint64_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Points",
        std::make_shared<writer::ElementBuffer<double>>(cellCount,
                                      pointsPerCell,
                                      std::vector<std::size_t>{dimension},
                                      std::forward<F>(pointIndexer)),
        datatype::inferDatatype<double>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Types",
        std::make_shared<writer::ElementBuffer<uint64_t>>(
            cellCount,
            1,
            std::vector<std::size_t>(),
            [&](uint64_t* target, std::size_t index) { target[0] = type; }),
        datatype::inferDatatype<uint64_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Connectivity",
        std::make_shared<writer::ElementBuffer<uint64_t>>(
            pointCount,
            1,
            std::vector<std::size_t>(),
            [&](uint64_t* target, std::size_t index) { target[0] = index + cellOffset; }),
        datatype::inferDatatype<uint64_t>()));

    // the following definitely need to stay in the current file
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        fieldDataGroup, "Time", writer::WriteInline::create<double>(time)));
    for (int quantity; quantity < NUMBER_OF_QUANTITIES; ++quantity) {
      std::string quantityName = std::string("q") + std::to_string(quantity + 1);
      // needs to use a managed buffer as well :/
      writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          pointDataGroup,
          quantityName,
          std::make_shared<writer::ElementBuffer<double>>(cellCount,
                                        pointsPerCell,
                                        std::vector<std::size_t>(),
                                        std::forward<G>(projector)),
          datatype::inferDatatype<double>()));
    }
    return writer;
  };
}
} // namespace seissol::io::instance
