#pragma once

#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <vector>

namespace seissol::io::reader::file {
class BinaryReader {
  public:
  BinaryReader(MPI_Comm comm);
  void openFile(const std::string& name);

  void readDataCollectiveRaw(void* data,
                             const std::string& name,
                             std::size_t count,
                             std::shared_ptr<datatype::Datatype> targetType,
                             bool seek);
  void readDataGlobalRaw(void* data,
                         const std::string& name,
                         std::size_t count,
                         std::shared_ptr<datatype::Datatype> targetType,
                         bool seek);
  void closeGroup();
  void closeFile();

  private:
  MPI_Comm comm;
};
} // namespace seissol::io::reader::file
