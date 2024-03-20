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
    if (time == 0) {
      /*writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location("checkpoint.h5", {"checkpoint"}),
          "globalId",
          writer::WriteBuffer::create(
              memoryManager.getLtsTree()->var(memoryManager.getLts()->dofs),
              memoryManager.getLtsTree()->getNumberOfCells(memoryManager.getLts()->dofs.mask)),
          datatype::inferDatatype<std::remove_pointer_t<decltype(memoryManager.getLtsTree()->var(
              memoryManager.getLts()->dofs))>>(),
          std::vector<std::size_t>()));*/
    }
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location("checkpoint.h5", {"checkpoint"}),
        "dofs",
        variableBuffer(*memoryManager.getLtsTree(), memoryManager.getLts()->dofs),
        datatype::inferDatatype<std::remove_pointer_t<decltype(memoryManager.getLtsTree()->var(
            memoryManager.getLts()->dofs))>>(),
        std::vector<std::size_t>()));
    if constexpr (kernels::size<tensor::Qane>() > 0) {
      writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location("checkpoint.h5", {"checkpoint"}),
          "dofsAne",
          variableBuffer(*memoryManager.getLtsTree(), memoryManager.getLts()->dofsAne),
          datatype::inferDatatype<std::remove_pointer_t<decltype(memoryManager.getLtsTree()->var(
              memoryManager.getLts()->dofsAne))>>(),
          std::vector<std::size_t>()));
    }
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location("checkpoint.h5", {"checkpoint"}),
        "pstrain",
        variableBuffer(*memoryManager.getLtsTree(), memoryManager.getLts()->pstrain),
        datatype::inferDatatype<std::remove_pointer_t<decltype(memoryManager.getLtsTree()->var(
            memoryManager.getLts()->pstrain))>>(),
        std::vector<std::size_t>()));
    return writer;
  };
}
} // namespace seissol::io::instance
