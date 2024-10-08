#pragma once

#include "Parallel/Pin.h"
#include <string>
#include <vector>

namespace seissol::writer {
class ThreadsPinningWriter {
  public:
  ThreadsPinningWriter(const std::string& outputDirectory) : outputDirectory(outputDirectory) {}
  void write(const seissol::parallel::Pinning& pinning);

  private:
  std::string outputDirectory;
};
} // namespace seissol::writer
