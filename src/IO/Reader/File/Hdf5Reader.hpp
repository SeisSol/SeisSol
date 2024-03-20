#pragma once

#include <IO/Datatype/Datatype.hpp>
#include <hdf5.h>
#include <memory>
#include <stack>
#include <string>

namespace seissol::io::reader::file {
class Hdf5Reader {
  public:
  Hdf5Reader(MPI_Comm comm);
  void openFile(MPI_Comm comm, const std::string& name);
  void openGroup(const std::string& name);
  void readAttribute(const std::string& name);
  void readData(const std::string& name,
                std::size_t count,
                std::shared_ptr<datatype::Datatype> targetType);
  void closeGroup();
  void closeFile();

  private:
  std::stack<hid_t> handles;
};
} // namespace seissol::io::reader::file
