// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_PARALLEL_MPI_H_
#define SEISSOL_SRC_PARALLEL_MPI_H_

#include <functional>
#ifndef USE_MPI
#include "MPIDummy.h"
#else // USE_MPI

#include "MPIBasic.h"
#include "utils/logger.h"
#include <algorithm>
#include <mpi.h>
#include <numeric>
#include <optional>
#include <string>

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
  public:
  ~MPI() override = default;

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

  void printAcceleratorDeviceInfo();

  /**
   * Initialize MPI
   */
  void init(int& argc, char**& argv);

  void setComm(MPI_Comm comm);

  template <typename T>
  [[nodiscard]] MPI_Datatype castToMpiType() const {
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
  [[nodiscard]] auto collect(T value, std::optional<MPI_Comm> comm = {}) const {
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
  [[nodiscard]] std::vector<ContainerType>
      collectContainer(const ContainerType& container, std::optional<MPI_Comm> comm = {}) const {
    using InternalType = typename ContainerType::value_type;

    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(m_comm);
    }
    int commSize{};
    MPI_Comm_size(comm.value(), &commSize);

    auto length = static_cast<int>(container.size());
    auto lengths = collect(length);

    std::vector<int> displacements(commSize);
    displacements[0] = 0;
    for (std::size_t i = 1; i < displacements.size(); ++i) {
      displacements[i] = displacements[i - 1] + lengths[i - 1];
    }

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

  // executes an operation on the given MPI communicator in serial order, i.e. first it is run by
  // rank 0, then by rank 1, then by rank 2 etc.
  template <typename F>
  void serialOrderExecute(F&& operation, std::optional<MPI_Comm> comm = {}) {
    if (!comm.has_value()) {
      comm = std::optional<MPI_Comm>(m_comm);
    }

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(comm.value(), &rank);
    MPI_Comm_size(comm.value(), &size);

    const int tag = 15; // TODO(David): replace by a tag allocation system one day
    char flag = 0;
    if (rank > 0) {
      MPI_Recv(&flag, 1, MPI_CHAR, rank - 1, tag, comm.value(), MPI_STATUS_IGNORE);
    }

    std::invoke(std::forward<F>(operation));

    // size >= 1 is ensured
    if (rank < size - 1) {
      MPI_Send(&flag, 1, MPI_CHAR, rank + 1, tag, comm.value());
    }
  }

  /**
   * sends a value to all processors
   */
  template <typename T>
  void broadcast(T* value, int root, std::optional<MPI_Comm> comm = {}) const {
    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(m_comm);
    }
    auto mpiType = castToMpiType<T>();
    MPI_Bcast(value, 1, mpiType, root, comm.value());
  }

  /**
   * sends a container to all processors
   */
  template <typename ContainerType>
  void broadcastContainer(ContainerType& container,
                          int root,
                          std::optional<MPI_Comm> comm = {}) const {
    using InternalType = typename ContainerType::value_type;
    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(m_comm);
    }
    auto size = static_cast<unsigned>(container.size());
    broadcast(&size, root);
    if (m_rank != root) {
      container.resize(size);
    }
    auto mpiType = castToMpiType<InternalType>();
    MPI_Bcast(const_cast<InternalType*>(container.data()), size, mpiType, root, comm.value());
  }

  /**
   * @return The main communicator for the application
   */
  [[nodiscard]] MPI_Comm comm() const { return m_comm; }

  /**
   * @return The node communicator (shared memory) for the application
   */
  [[nodiscard]] MPI_Comm sharedMemComm() const { return m_sharedMemComm; }

  /**
   * @return hostnames for all ranks in the communicator of the application
   */
  const auto& getHostNames() { return hostNames; }

  static void barrier(MPI_Comm comm) { MPI_Barrier(comm); }

  /**
   * Finalize MPI
   */
  static void finalize() { MPI_Finalize(); }

  void setDataTransferModeFromEnv();

  enum class DataTransferMode { Direct, CopyInCopyOutHost };
  DataTransferMode getPreferredDataTransferMode() { return preferredDataTransferMode; }

  /** The only instance of the class */
  static MPI mpi;

  private:
  MPI_Comm m_comm;
  MPI_Comm m_sharedMemComm{};
  MPI() : m_comm(MPI_COMM_NULL) {}
  DataTransferMode preferredDataTransferMode{DataTransferMode::Direct};
  std::vector<std::string> hostNames;
};

#endif // USE_MPI

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_MPI_H_
