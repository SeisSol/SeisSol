#pragma once
#include <IO/Writer/Instructions/Binary.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::io::writer::file {

class BinaryFile {
  public:
  BinaryFile(MPI_Comm comm);
  void openFile(const std::string& name);
  void writeGlobal(const void* data, std::size_t size);
  void writeDistributed(const void* data, std::size_t size);
  void closeFile();

  private:
  MPI_Comm comm;
  MPI_File file;
};

class BinaryWriter {
  public:
  BinaryWriter(MPI_Comm comm);

  void write(const async::ExecInfo& info, const instructions::BinaryWrite& write);

  void finalize();

  private:
  std::unordered_map<std::string, std::unique_ptr<BinaryFile>> openFiles;
  MPI_Comm comm;
};
} // namespace seissol::io::writer::file
