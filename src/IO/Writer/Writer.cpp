// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Writer.h"

#include "IO/Writer/File/BinaryWriter.h"
#include "IO/Writer/File/Hdf5Writer.h"
#include "IO/Writer/File/RunFiles.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "Instructions/Instruction.h"

#include <async/ExecInfo.h>
#include <memory>
#include <mpi.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

WriteInstance::WriteInstance(MPI_Comm comm, file::RunFiles* runFiles)
    : hdf5_(comm, runFiles), binary_(comm, runFiles) {}

void WriteInstance::write(
    const async::ExecInfo& info,
    const std::vector<std::shared_ptr<instructions::WriteInstruction>>& plan) {
  // A binary file that takes no distributed data in the plan only takes what every rank holds
  // alike, which rank 0 writes alone -- the XML of an Xdmf or PVD output. The others need not open
  // it at all; opening is collective otherwise, and costs every rank a file system operation.
  std::set<std::string> distributedFiles;
  for (const auto& instruction : plan) {
    if (const auto* binary = dynamic_cast<instructions::BinaryWrite*>(instruction.get());
        binary != nullptr && binary->dataSource->distributed()) {
      distributedFiles.insert(binary->filename);
    }
  }
  binary_.beginWrite(distributedFiles);
  for (const auto& instruction : plan) {
    write(info, instruction);
  }
  hdf5_.endWrite();
  binary_.endWrite();
}

void WriteInstance::write(const async::ExecInfo& info,
                          const std::shared_ptr<instructions::WriteInstruction>& instruction) {
  if (dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()) != nullptr) {
    hdf5_.writeData(info, *dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()));
  }
  if (dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()) != nullptr) {
    hdf5_.writeAttribute(info, *dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()));
  }
  if (dynamic_cast<instructions::Hdf5LinkExternalWrite*>(instruction.get()) != nullptr) {
    hdf5_.writeLinkExternal(info,
                            *dynamic_cast<instructions::Hdf5LinkExternalWrite*>(instruction.get()));
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

void Writer::write(const async::ExecInfo& info, WriteInstance& instance) const {
  instance.write(info, instructions_);
}

const std::vector<std::shared_ptr<instructions::WriteInstruction>>&
    Writer::getInstructions() const {
  return instructions_;
}

} // namespace seissol::io::writer
