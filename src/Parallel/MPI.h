// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_PARALLEL_MPI_H_
#define SEISSOL_SRC_PARALLEL_MPI_H_

#include "Common/Real.h"
#include "Kernels/Precision.h"

#include <algorithm>
#include <functional>
#include <mpi.h>
#include <numeric>
#include <optional>
#include <string>
#include <type_traits>
#include <typeindex>
#include <unordered_map>
#include <utils/logger.h>

#ifdef MPI_VERSION
#if MPI_VERSION >= 4
#define SEISSOL_MPI_LARGEINT
#endif
#endif

#ifdef SEISSOL_MPI_LARGEINT
#define SEISSOL_MPI_C(fun) fun##_c
#else
#define SEISSOL_MPI_C(fun) fun
#endif

namespace seissol {

/**
 * Wraps an MPI communicator. (and currently still a node-local communicator)
 * Provides convenience methods for communication and MPI type inference.
 *
 * Currently still mostly used as a singleton; TODO: change that.
 */
class Mpi {
  private:
  /** This rank */
  int rank_{0};

  /** Rank in the shared memory sub-communicator */
  int sharedMemMpiRank_{0};

  /** Number of processors */
  int size_{1};

  /** Number of ranks in the shared memory sub-communicator */
  int sharedMemMpiSize_{1};

  public:
  /**
   * @return The rank of this process
   */
  [[nodiscard]] int rank() const { return rank_; }

  /**
   * @return The rank within the shared memory sub-communicator
   */
  [[nodiscard]] int sharedMemMpiRank() const { return sharedMemMpiRank_; }

  /**
   * @return The total number of processes
   */
  [[nodiscard]] int size() const { return size_; }

  /**
   * @return The number of ranks within the shared memory sub-communicator
   */
  [[nodiscard]] int sharedMemMpiSize() const { return sharedMemMpiSize_; }

  [[nodiscard]] bool isSingleProcess() const { return size() == 1; }

  [[nodiscard]] bool isSingleNode() const { return size() == sharedMemMpiSize(); }

#ifdef SEISSOL_MPI_LARGEINT
  using Count = MPI_Count;
#else
  using Count = int;
#endif

  /**
   * Initialize MPI
   */
  void init(int& argc, char**& argv);

  void setComm(MPI_Comm comm);

  template <typename T>
  [[nodiscard]] static MPI_Datatype castToMpiType() {
    if constexpr (std::is_same_v<T, double>) {
      return MPI_DOUBLE;
    } else if constexpr (std::is_same_v<T, float>) {
      return MPI_FLOAT;
    } else if constexpr (std::is_same_v<T, unsigned>) {
      return MPI_UNSIGNED;
    } else if constexpr (std::is_same_v<T, int>) {
      return MPI_INT;
    } else if constexpr (std::is_same_v<T, unsigned long>) {
      return MPI_UNSIGNED_LONG;
    } else if constexpr (std::is_same_v<T, long>) {
      return MPI_LONG;
    } else if constexpr (std::is_same_v<T, unsigned long long>) {
      return MPI_UNSIGNED_LONG_LONG;
    } else if constexpr (std::is_same_v<T, long long>) {
      return MPI_LONG_LONG;
    } else if constexpr (std::is_same_v<T, short>) {
      return MPI_SHORT;
    } else if constexpr (std::is_same_v<T, char>) {
      return MPI_CHAR;
    } else if constexpr (std::is_same_v<T, unsigned short>) {
      return MPI_UNSIGNED_SHORT;
    } else if constexpr (std::is_same_v<T, unsigned char>) {
      return MPI_UNSIGNED_CHAR;
    } else if constexpr (std::is_same_v<T, bool>) {
      return MPI_C_BOOL;
    } else if constexpr (std::is_same_v<T, std::pair<float, int>>) {
      return MPI_FLOAT_INT;
    } else if constexpr (std::is_same_v<T, std::pair<long, int>>) {
      return MPI_LONG_INT;
    } else if constexpr (std::is_same_v<T, std::pair<double, int>>) {
      return MPI_DOUBLE_INT;
    } else if constexpr (std::is_same_v<T, std::pair<short, int>>) {
      return MPI_SHORT_INT;
    } else if constexpr (std::is_same_v<T, std::pair<int, int>>) {
      return MPI_2INT;
    } else {
      static_assert(sizeof(T) == 0, "Unimplemented MPI type.");
      // return something to make NVHPC happy
      return MPI_BYTE;
    }
  }

  [[nodiscard]] static MPI_Datatype precisionToMpiType(RealType type) {
    switch (type) {
    case seissol::RealType::F32:
      return MPI_FLOAT;
    case seissol::RealType::F64:
      return MPI_DOUBLE;
    default:
      return MPI_BYTE;
    }
  }

  template <typename T>
  void registerMpiType(MPI_Datatype type, bool commit = false) {
    const auto typeIndex = std::type_index(typeid(T));
    cachedTypes[typeIndex] = type;
    if (commit) {
      MPI_Type_commit(&type);
    }
  }

  template <typename T>
  MPI_Datatype getMpiType() {
    const auto typeIndex = std::type_index(typeid(T));
    const auto cacheFind = cachedTypes.find(typeIndex);
    if (cacheFind != cachedTypes.end()) {
      return cacheFind->second;
    }
    logError() << "Requested type not (dynamically) registered.";
    return MPI_BYTE;
  }

  /**
   * @return a vector of aggregated values from all processes
   *
   * Note: the vector is available only in the master process.
   * Method supports only basic types
   */
  template <typename T>
  [[nodiscard]] auto collect(T value, std::optional<MPI_Comm> comm = {}) const {
    auto collect = std::vector<T>(size_);
    auto type = castToMpiType<T>();
    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(comm_);
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
      comm = std::optional<MPI_Comm>(comm_);
    }
    int commSize{};
    MPI_Comm_size(comm.value(), &commSize);

    auto length = static_cast<int>(container.size());
    auto lengths = collect(length, comm);

    // including the total length
    std::vector<int> displacements(commSize + 1);
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

    std::vector<ContainerType> collected(commSize);
    for (int rank = 0; rank < commSize; ++rank) {
      const auto begin = displacements[rank];
      const auto end = displacements[rank + 1];
      collected[rank] = ContainerType(recvBuffer.begin() + begin, recvBuffer.begin() + end);
    }
    return collected;
  }

  // executes an operation on the given MPI communicator in serial order, i.e. first it is run by
  // rank 0, then by rank 1, then by rank 2 etc.
  template <typename F>
  void serialOrderExecute(F&& operation, std::optional<MPI_Comm> comm = {}) {
    if (!comm.has_value()) {
      comm = std::optional<MPI_Comm>(comm_);
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
      comm = std::optional<MPI_Comm>(comm_);
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
      comm = std::optional<MPI_Comm>(comm_);
    }
    auto size = static_cast<unsigned>(container.size());
    broadcast(&size, root);
    if (rank_ != root) {
      container.resize(size);
    }
    auto mpiType = castToMpiType<InternalType>();
    MPI_Bcast(const_cast<InternalType*>(container.data()), size, mpiType, root, comm.value());
  }

  template <typename T>
  T allreduce(const T& value, MPI_Op op, std::optional<MPI_Comm> comm = {}) const {
    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(comm_);
    }

    T result{};
    const auto mpiType = castToMpiType<T>();
    MPI_Allreduce(&value, &result, 1, mpiType, op, comm.value());
    return result;
  }

  template <typename ContainerType>
  void allreduceContainer(ContainerType& container,
                          MPI_Op op,
                          std::optional<MPI_Comm> comm = {}) const {
    using InternalType = typename ContainerType::value_type;
    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(comm_);
    }
    const auto size = static_cast<unsigned>(container.size());
    const auto mpiType = castToMpiType<InternalType>();
    MPI_Allreduce(
        MPI_IN_PLACE, const_cast<InternalType*>(container.data()), size, mpiType, op, comm.value());
  }

  template <typename T>
  T scan(const T& value, MPI_Op op, bool inclusive, std::optional<MPI_Comm> comm = {}) const {
    if (not comm.has_value()) {
      comm = std::optional<MPI_Comm>(comm_);
    }

    T result{};
    auto mpiType = castToMpiType<T>();
    if (inclusive) {
      MPI_Scan(&value, &result, 1, mpiType, op, comm.value());
    } else {
      MPI_Exscan(&value, &result, 1, mpiType, op, comm.value());
    }

    return result;
  }

  /**
   * @return The main communicator for the application
   */
  [[nodiscard]] MPI_Comm comm() const { return comm_; }

  /**
   * @return The node communicator (shared memory) for the application
   */
  [[nodiscard]] MPI_Comm sharedMemComm() const { return sharedMemComm_; }

  void barrier(std::optional<MPI_Comm> comm = {}) {
    if (!comm.has_value()) {
      comm = std::optional<MPI_Comm>(comm_);
    }
    MPI_Barrier(comm.value());
  }

  /**
   * Finalize MPI
   */
  void finalize() {
    for (auto& [_, type] : cachedTypes) {
      MPI_Type_free(&type);
    }
    MPI_Finalize();
  }

  /** The only instance of the class */
  static Mpi mpi;

  private:
  MPI_Comm comm_{MPI_COMM_NULL};
  MPI_Comm sharedMemComm_{};
  Mpi() = default;
  std::unordered_map<std::type_index, MPI_Datatype> cachedTypes;
};

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_MPI_H_
