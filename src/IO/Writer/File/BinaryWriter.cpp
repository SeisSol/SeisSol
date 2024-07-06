#include "BinaryWriter.hpp"
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

BinaryFile::BinaryFile(MPI_Comm comm) {}
void BinaryFile::openFile(const std::string& name) {}
void BinaryFile::writeGlobal(void* data, std::size_t size) {}
void BinaryFile::writeDistributed(void* data, std::size_t size) {}
void BinaryFile::closeFile() {}

BinaryWriter::BinaryWriter(MPI_Comm comm) {}

void BinaryWriter::write(const async::ExecInfo& info, const instructions::BinaryWrite& write) {}

void BinaryWriter::finalize() {}

} // namespace seissol::io::writer::file
