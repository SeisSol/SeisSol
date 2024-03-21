#pragma once

#include <IO/Datatype/Inference.hpp>
#include <IO/Datatype/MPIType.hpp>
#include <cstddef>
#include <mpi.h>
#include <vector>

namespace seissol::io::reader {
class Distributor {
  public:
  Distributor(MPI_Comm comm);

  void setup(std::size_t globalCount,
             const std::vector<std::size_t>& sourceIds,
             const std::vector<std::size_t>& targetIds);

  template <typename T>
  void distribute(T* target,
                  const T* source,
                  MPI_Datatype datatype = seissol::io::datatype::convertToMPI(
                      seissol::io::datatype::inferDatatype<T>())) {
    return distributeInternal(target, source, datatype);
  }

  private:
  // distributes data. Note that in-place operations are supported.
  void distributeInternal(void* target, const void* source, MPI_Datatype datatype);

  std::vector<std::size_t> sendOffsets;
  std::vector<std::size_t> recvOffsets;
  std::vector<std::size_t> sendReorder;
  std::vector<std::size_t> recvReorder;
  MPI_Comm comm;
};
} // namespace seissol::io::reader
