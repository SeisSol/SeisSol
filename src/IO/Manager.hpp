#pragma once

#include "Writer/Module/WriterModule.hpp"
#include <IO/Writer/Writer.hpp>
#include <Modules/Module.h>
#include <Modules/Modules.h>
#include <memory>
#include <vector>

namespace seissol::io {

class OutputManager : public seissol::Module {
  public:
  OutputManager() = default;

  void setup() { Modules::registerHook(*this, ModuleHook::PostMPIInit); }

  void postMPIInit() override {}

  void addOutput(const writer::ScheduledWriter& writer) {
    modules.emplace_back(std::make_unique<seissol::io::writer::module::WriterModule>(writer));
    modules.back()->startup();
  }

  void loadCheckpoint(const std::string& path) {}

  private:
  std::vector<std::unique_ptr<seissol::io::writer::module::WriterModule>> modules;
};

} // namespace seissol::io
