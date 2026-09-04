// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Writer.h"

#include "IO/Writer/File/BinaryWriter.h"
#include "IO/Writer/File/Hdf5Writer.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "Instructions/Instruction.h"

#include <async/ExecInfo.h>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

WriteInstance::WriteInstance(MPI_Comm comm) : hdf5_(comm), binary_(comm) {}

void WriteInstance::write(const async::ExecInfo& info,
                          const std::shared_ptr<instructions::WriteInstruction>& instruction) {
  if (dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()) != nullptr) {
    hdf5_.writeData(info, *dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()));
  }
  if (dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()) != nullptr) {
    hdf5_.writeAttribute(info, *dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()));
  }
  if (dynamic_cast<instructions::BinaryWrite*>(instruction.get()) != nullptr) {
    binary_.write(info, *dynamic_cast<instructions::BinaryWrite*>(instruction.get()));
  }
}

void WriteInstance::close() {
  hdf5_.finalize();
  binary_.finalize();
}

Writer::Writer() = default;

Writer::Writer(const std::string& data) {
  const YAML::Node plan = YAML::Load(data);
  for (const YAML::Node& instruction : plan) {
    instructions_.push_back(instructions::WriteInstruction::deserialize(instruction));
  }
}

void Writer::addInstruction(const std::shared_ptr<instructions::WriteInstruction>& instruction) {
  instructions_.push_back(instruction);
}

std::string Writer::serialize() {
  std::stringstream sstr;
  {
    YAML::Emitter output(sstr);
    output << YAML::BeginSeq;
    for (const auto& instruction : instructions_) {
      output << instruction->serialize();
    }
    output << YAML::EndSeq;
  }
  return sstr.str();
}

WriteInstance Writer::beginWrite(const async::ExecInfo& info, MPI_Comm comm) {
  WriteInstance instance(comm);
  for (const auto& instruction : instructions_) {
    instance.write(info, instruction);
  }
  return instance;
}

void Writer::endWrite() {}

const std::vector<std::shared_ptr<instructions::WriteInstruction>>&
    Writer::getInstructions() const {
  return instructions_;
}

} // namespace seissol::io::writer
