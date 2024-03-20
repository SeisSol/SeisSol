#pragma once

#include "Instructions/Instruction.hpp"
#include "async/ExecInfo.h"
#include <IO/Writer/File/Hdf5Writer.hpp>
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <memory>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

class WriteInstance {
  public:
  WriteInstance(MPI_Comm comm) : hdf5(comm) {}

  void write(const async::ExecInfo& info,
             std::shared_ptr<instructions::WriteInstruction> instruction) {
    if (dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()) != nullptr) {
      hdf5.writeData(info, *dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()));
    }
    if (dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()) != nullptr) {
      hdf5.writeAttribute(info,
                          *dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get()));
    }
  }

  void close() { hdf5.finalize(); }

  private:
  file::Hdf5Writer hdf5;
};

class Writer {
  public:
  Writer() = default;

  explicit Writer(const std::string& data) {
    YAML::Node plan = YAML::Load(data);
    for (const auto& instruction : plan) {
      instructions.push_back(instructions::WriteInstruction::deserialize(instruction));
    }
  }

  void addInstruction(std::shared_ptr<instructions::WriteInstruction> instruction) {
    instructions.push_back(instruction);
  }

  std::string serialize() {
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

  WriteInstance beginWrite(const async::ExecInfo& info) {
    WriteInstance instance(MPI_COMM_WORLD);
    for (auto instruction : instructions) {
      instance.write(info, instruction);
    }
    return instance;
  }

  void endWrite() {}

  const std::vector<std::shared_ptr<instructions::WriteInstruction>>& getInstructions() const {
    return instructions;
  }

  private:
  std::vector<std::shared_ptr<instructions::WriteInstruction>> instructions;
};

struct ScheduledWriter {
  std::string name;
  double interval;
  std::function<Writer(double)> planWrite;
};

} // namespace seissol::io::writer
