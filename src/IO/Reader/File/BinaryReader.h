// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_READER_FILE_BINARYREADER_HPP_
#define SEISSOL_SRC_IO_READER_FILE_BINARYREADER_HPP_

#include <IO/Datatype/Datatype.hpp>
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
  MPI_File file;
};
} // namespace seissol::io::reader::file

#endif // SEISSOL_SRC_IO_READER_FILE_BINARYREADER_HPP_
