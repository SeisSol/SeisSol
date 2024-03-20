#pragma once

#include <string>
namespace seissol::io::writer::file {
class File {
  public:
  std::string filename() { return name; }

  private:
  std::string name;
};
} // namespace seissol::io::writer::file
