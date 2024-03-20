#pragma once

#include <IO/Datatype/Inference.hpp>
#include <IO/Instance/SeisSolMemoryHelper.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <IO/Writer/Writer.hpp>
#include <Initializer/MemoryManager.h>
#include <functional>
#include <type_traits>
namespace seissol::io::instance {
inline std::function<writer::Writer(double)>
    makeCheckpointWriter(initializer::MemoryManager& memoryManager) {
  return [&](double time) -> writer::Writer {
    writer::Writer writer;
    auto mainGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF"});
    auto cellDataGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF", "CellData"});
    auto pointDataGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF", "PointData"});
    auto fieldDataGroup = writer::instructions::Hdf5Location("vtkhdf.h5", {"VTKHDF", "FieldData"});

    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        mainGroup, "Type", writer::WriteInline::create("UnstructuredGrid")));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        mainGroup, "Version", writer::WriteInline::create<uint64_t[2]>({1, 0})));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "NumberOfCells",
        writer::WriteInline::create<uint64_t[3]>({1, 0, 0}),
        datatype::inferDatatype<uint64_t[3]>(),
        std::vector<std::size_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "NumberOfConnectivityIds",
        writer::WriteInline::create<uint64_t[3]>({1, 0, 0}),
        datatype::inferDatatype<uint64_t[3]>(),
        std::vector<std::size_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "NumberOfPoints",
        writer::WriteInline::create<uint64_t[3]>({1, 0, 0}),
        datatype::inferDatatype<uint64_t[3]>(),
        std::vector<std::size_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Offsets",
        writer::WriteInline::create<uint64_t[3]>({1, 0, 0}),
        datatype::inferDatatype<uint64_t[3]>(),
        std::vector<std::size_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Points",
        writer::WriteInline::create<uint64_t[3]>({1, 0, 0}),
        datatype::inferDatatype<uint64_t[3]>(),
        std::vector<std::size_t>()));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        mainGroup,
        "Types",
        writer::WriteInline::create<uint64_t[3]>({1, 0, 0}),
        datatype::inferDatatype<uint64_t[3]>(),
        std::vector<std::size_t>()));

    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        fieldDataGroup,
        "Time",
        writer::WriteInline::create<double>(time),
        datatype::inferDatatype<double>(),
        std::vector<std::size_t>()));
    for (int quantity; quantity < NUMBER_OF_QUANTITIES; ++quantity) {
      std::string quantityName = std::string("q") + std::to_string(quantity + 1);
      writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          pointDataGroup,
          quantityName,
          writer::WriteInline::create<uint64_t[3]>({1, 0, 0}),
          datatype::inferDatatype<double>(),
          std::vector<std::size_t>{2}));
    }
    return writer;
  };
}
} // namespace seissol::io::instance
