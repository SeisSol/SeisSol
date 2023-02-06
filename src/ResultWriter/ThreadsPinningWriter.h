#pragma once

#include "Parallel/Pin.h"
#include <vector>
#include <string>

namespace seissol::writer {
class ThreadsPinningWriter {
  public:
  ThreadsPinningWriter(const std::string& outputDirectory) : outputDirectory(outputDirectory) {}
  void write(seissol::parallel::Pinning& pinning);

  private:
  std::string outputDirectory;
};
} // namespace seissol::writer