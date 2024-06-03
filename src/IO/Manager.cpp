#include "Manager.hpp"

#include "Writer/Module/WriterModule.hpp"
#include <IO/Writer/Writer.hpp>
#include <SeisSol.h>
#include <memory>
#include <vector>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::io {

OutputManager::OutputManager(SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

void OutputManager::addOutput(const writer::ScheduledWriter& writer) {
  modules.emplace_back(std::make_unique<seissol::io::writer::module::WriterModule>(
      seissolInstance.getSeisSolParameters().output.prefix, writer, seissolInstance.getPinning()));
  modules.back()->startup();
}

} // namespace seissol::io
