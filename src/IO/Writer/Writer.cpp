// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Writer.hpp"

#include "Instructions/Instruction.hpp"
#include "async/ExecInfo.h"
#include <IO/Writer/File/BinaryWriter.hpp>
#include <IO/Writer/File/Hdf5Writer.hpp>
#include <IO/Writer/Instructions/Binary.hpp>
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

WriteInstance::WriteInstance(MPI_Comm comm) : hdf5(comm), binary(comm) {}

void WriteInstance::write(const async::ExecInfo& info,
                          std::shared_ptr<instructions::WriteInstruction> instruction) {
  if (dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()) != nullptr) {
    hdf5.writeData(info, *dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()));
  }
  if (dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()) != nullptr) {
    hdf5.writeAttribute(info, *dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()));
  }
  if (dynamic_cast<instructions::BinaryWrite*>(instruction.get()) != nullptr) {
    binary.write(info, *dynamic_cast<instructions::BinaryWrite*>(instruction.get()));
  }
}

void WriteInstance::close() {
  hdf5.finalize();
  binary.finalize();
}

Writer::Writer() = default;

Writer::Writer(const std::string& data) {
  const YAML::Node plan = YAML::Load(data);
  for (const auto& instruction : plan) {
    instructions.push_back(instructions::WriteInstruction::deserialize(instruction));
  }
}

void Writer::addInstruction(std::shared_ptr<instructions::WriteInstruction> instruction) {
  instructions.push_back(instruction);
}

std::string Writer::serialize() {
  std::stringstream sstr;
  {
    YAML::Emitter output(sstr);
    output << YAML::BeginSeq;
    for (const auto& instruction : instructions) {
      output << instruction->serialize();
    }
    output << YAML::EndSeq;
  }
  return sstr.str();
}

WriteInstance Writer::beginWrite(const async::ExecInfo& info) {
  WriteInstance instance(MPI_COMM_WORLD);
  for (auto instruction : instructions) {
    instance.write(info, instruction);
  }
  return instance;
}

void Writer::endWrite() {}

const std::vector<std::shared_ptr<instructions::WriteInstruction>>&
    Writer::getInstructions() const {
  return instructions;
}

} // namespace seissol::io::writer
