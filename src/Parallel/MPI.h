/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * MPI Wrapper
 */

#ifndef MPI_H
#define MPI_H

#ifndef USE_MPI
#include "MPIDummy.h"
#else // USE_MPI

#include <mpi.h>
#include "utils/logger.h"
#include "MPIBasic.h"
#include <numeric>
#include <algorithm>
#include <string>
#include <optional>

#endif // USE_MPI

namespace seissol {

#ifndef USE_MPI
typedef MPIDummy MPI;
#else // USE_MPI

/**
 * MPI handling.
 *
 * Make sure only one instance of this class exists!
 */
class MPI : public MPIBasic {
  private:
  MPI_Comm m_comm;
  MPI_Comm m_sharedMemComm;
  MPI() : m_comm(MPI_COMM_NULL) {}

  std::vector<std::string> hostNames{};

  public:
  ~MPI() {}

#ifdef ACL_DEVICE
  /**
   * @brief Inits Device(s).
   *
   * Some MPI implementations create a so-called context between GPUs and OS Processes inside of
   * MPI_Init(...). It results in allocating some memory buffers in memory attached to the nearest
   * NUMA domain of a core where a process is running. In case of somebody wants to bind a processes
   * in a different way, e.g. move a process closer to a GPU, it must be done before calling
   * MPI_Init(...) using env. variables or hwloc library.
   *
   * Currently, the function does a simple binding, i.e. it binds to the first visible device.
   * The user is responsible for the correct binding on a multi-gpu setup.
   * One can use a wrapper script and manipulate with CUDA_VISIBLE_DEVICES/HIP_VISIBLE_DEVICES and
   * OMPI_COMM_WORLD_LOCAL_RANK env. variables
   * */
  void bindAcceleratorDevice();
#endif // ACL_DEVICE

  /**
   * Initialize MPI
   */
  void init(int& argc, char**& argv);

  void setComm(MPI_Comm comm);

  template <typename T>
  MPI_Datatype castToMpiType() {
    if constexpr (std::is_same_v<T, double>) {
      return MPI_DOUBLE;
    } else if constexpr (std::is_same_v<T, float>) {
      return MPI_FLOAT;
    } else if constexpr (std::is_same_v<T, unsigned>) {
      return MPI_UNSIGNED;
    } else if constexpr (std::is_same_v<T, int>) {
      return MPI_INT;
    } else if constexpr (std::is_same_v<T, char>) {
      return MPI_CHAR;
    } else if constexpr (std::is_same_v<T, bool>) {
      return MPI_C_BOOL;
    }
  }

  /**
   * @return a vector of aggregated values from all processes
   *
   * Note: the vector is available only in the master process.
   * Method supports only basic types
   */
  template <typename T>
  auto collect(T value, std::optional<MPI_Comm> comm = {}) {
    auto collect = std::vector<T>(m_size);
    auto type = castToMpiType<T>();
    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(m_comm);
    }
    MPI_Gather(&value, 1, type, collect.data(), 1, type, 0, comm.value());
    return collect;
  }

  /**
   * @return a vector of aggregated containers from all processes
   *
   * Note: the vector is available only in the master process
   * Method supports only basic types
   */
  template <typename ContainerType>
  std::vector<ContainerType> collectContainer(const ContainerType& container,
                                              std::optional<MPI_Comm> comm = {}) {
    using InternalType = typename ContainerType::value_type;

    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(m_comm);
    }
    int commSize{};
    MPI_Comm_size(comm.value(), &commSize);

    auto length = static_cast<int>(container.size());
    auto lengths = collect(length);

    std::vector<int> displacements(commSize);
    std::exclusive_scan(lengths.begin(), lengths.end(), displacements.begin(), 0);

    const auto recvBufferSize = std::accumulate(lengths.begin(), lengths.end(), 0);
    std::vector<InternalType> recvBuffer(recvBufferSize);

    auto mpiType = castToMpiType<InternalType>();
    MPI_Gatherv(container.data(),
                container.size(),
                mpiType,
                const_cast<InternalType*>(recvBuffer.data()),
                lengths.data(),
                displacements.data(),
                mpiType,
                0,
                comm.value());

    displacements.push_back(recvBufferSize);
    std::vector<ContainerType> collected(commSize);
    for (int rank = 0; rank < commSize; ++rank) {
      const auto begin = displacements[rank];
      const auto end = displacements[rank + 1];
      collected[rank] = ContainerType(&recvBuffer[begin], &recvBuffer[end]);
    }
    return collected;
  }

  /**
   * @return The main communicator for the application
   */
  MPI_Comm comm() const { return m_comm; }

  /**
   * @return The node communicator (shared memory) for the application
   */
  MPI_Comm sharedMemComm() const { return m_sharedMemComm; }

  /**
   * @return hostnames for all ranks in the communicator of the application
   */
  const auto& getHostNames() { return hostNames; }

  void barrier(MPI_Comm comm) const { MPI_Barrier(comm); }

  /**
   * Finalize MPI
   */
  void finalize() {
    fault.finalize();
    MPI_Finalize();
  }

  /** The only instance of the class */
  static MPI mpi;
};

#endif // USE_MPI

} // namespace seissol

#endif // MPI_H
