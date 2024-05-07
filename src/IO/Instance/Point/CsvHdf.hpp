#pragma once

#include <string>
namespace seissol::io::instance::mesh {

class CsvHdf {
  public:
  CsvHdf();

  void addQuantity(const std::string& name);

  private:
};

} // namespace seissol::io::instance::mesh
