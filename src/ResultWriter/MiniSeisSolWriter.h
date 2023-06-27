#pragma once

#include <vector>
#include <string>

namespace seissol::writer {
class MiniSeisSolWriter {
  public:
  MiniSeisSolWriter(const char* outputDirectory) : outputDirectory(outputDirectory) {}
  void write(double elapsedTime, double weight);

  private:
  std::string outputDirectory;
};
} // namespace seissol::writer