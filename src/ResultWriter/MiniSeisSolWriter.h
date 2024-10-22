#pragma once

#include <string>
#include <vector>

namespace seissol::writer {
class MiniSeisSolWriter {
  public:
  MiniSeisSolWriter(const char* outputDirectory) : outputDirectory(outputDirectory) {}
  void write(double elapsedTime, double weight);

  private:
  std::string outputDirectory;
};
} // namespace seissol::writer
