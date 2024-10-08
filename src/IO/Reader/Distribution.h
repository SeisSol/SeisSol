// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_READER_DISTRIBUTION_H_
#define SEISSOL_SRC_IO_READER_DISTRIBUTION_H_

#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <cstddef>
#include <functional>
#include <mpi.h>
#include <vector>

namespace seissol::io::reader {

class Distributor {
  public:
  class DistributionInstance {
public:
    DistributionInstance(const std::function<void()>& completion);
    void complete();

private:
    bool completed;
    std::function<void()> completion;
  };

  Distributor(MPI_Comm comm);

  void setup(const std::vector<std::size_t>& sourceIds, const std::vector<std::size_t>& targetIds);

  template <typename T>
  DistributionInstance distribute(T* target,
                                  const T* source,
                                  MPI_Datatype datatype = seissol::io::datatype::convertToMPI(
                                      seissol::io::datatype::inferDatatype<T>())) {
    return distributeInternal(target, source, datatype);
  }

  private:
  // distributes data. Note that in-place operations are supported.
  DistributionInstance distributeInternal(void* target, const void* source, MPI_Datatype datatype);

  std::vector<std::size_t> sendOffsets;
  std::vector<std::size_t> recvOffsets;
  std::vector<std::size_t> sendReorder;
  std::vector<std::size_t> recvReorder;
  MPI_Comm comm;
};
} // namespace seissol::io::reader

#endif // SEISSOL_SRC_IO_READER_DISTRIBUTION_H_
