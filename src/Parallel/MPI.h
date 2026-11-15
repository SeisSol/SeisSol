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

#include <algorithm>
#include <cstddef>
#include <functional>
#include <mpi.h>
#include <numeric>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

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
 * The methods that communicate take the communicator to use as their last argument; without it,
 * they use comm().
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

  /**
   * @return The rank of this process within the communicator @p comm
   */
  [[nodiscard]] static int rankOf(MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }

  /**
   * @return The number of processes in the communicator @p comm
   */
  [[nodiscard]] static int sizeOf(MPI_Comm comm) {
    int size = 0;
    MPI_Comm_size(comm, &size);
    return size;
  }

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

  /**
   * The MPI_<T>_INT pair types describe a plain struct of the two members. std::pair has no
   * guaranteed layout; passing it as such a type is only valid where it is laid out the same way.
   */
  template <typename First, typename Second>
  struct MpiPairLayout {
    First first;
    Second second;
  };

  template <typename First, typename Second>
  [[nodiscard]] static constexpr bool pairMatchesMpiLayout() {
    using PairT = std::pair<First, Second>;
    using ReferenceT = MpiPairLayout<First, Second>;
    return std::is_standard_layout_v<PairT> && sizeof(PairT) == sizeof(ReferenceT) &&
           offsetof(PairT, first) == offsetof(ReferenceT, first) &&
           offsetof(PairT, second) == offsetof(ReferenceT, second);
  }

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
      static_assert(pairMatchesMpiLayout<float, int>());
      return MPI_FLOAT_INT;
    } else if constexpr (std::is_same_v<T, std::pair<long, int>>) {
      static_assert(pairMatchesMpiLayout<long, int>());
      return MPI_LONG_INT;
    } else if constexpr (std::is_same_v<T, std::pair<double, int>>) {
      static_assert(pairMatchesMpiLayout<double, int>());
      return MPI_DOUBLE_INT;
    } else if constexpr (std::is_same_v<T, std::pair<short, int>>) {
      static_assert(pairMatchesMpiLayout<short, int>());
      return MPI_SHORT_INT;
    } else if constexpr (std::is_same_v<T, std::pair<int, int>>) {
      static_assert(pairMatchesMpiLayout<int, int>());
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

  /**
   * @return a vector of aggregated values from all processes
   *
   * Note: the vector is available only in the master process.
   * Method supports only basic types
   */
  template <typename T>
  [[nodiscard]] auto collect(T value, std::optional<MPI_Comm> comm = {}) const {
    const auto& useComm = comm.value_or(this->comm());
    auto collect = std::vector<T>(sizeOf(useComm));
    auto type = castToMpiType<T>();
    MPI_Gather(&value, 1, type, collect.data(), 1, type, 0, useComm);
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
    using InternalType = ContainerType::value_type;

    const auto& useComm = comm.value_or(this->comm());
    const int commSize = sizeOf(useComm);

    const auto length = static_cast<int>(container.size());
    auto lengths = collect(length, useComm);

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
                length,
                mpiType,
                recvBuffer.data(),
                lengths.data(),
                displacements.data(),
                mpiType,
                0,
                useComm);

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
  void serialOrderExecute(F&& operation, std::optional<MPI_Comm> comm = {}) const {
    const auto& useComm = comm.value_or(this->comm());
    const int rank = rankOf(useComm);
    const int size = sizeOf(useComm);

    const int tag = 15; // TODO(David): replace by a tag allocation system one day
    char flag = 0;
    if (rank > 0) {
      MPI_Recv(&flag, 1, MPI_CHAR, rank - 1, tag, useComm, MPI_STATUS_IGNORE);
    }

    std::invoke(std::forward<F>(operation));

    // size >= 1 is ensured
    if (rank < size - 1) {
      MPI_Send(&flag, 1, MPI_CHAR, rank + 1, tag, useComm);
    }
  }

  /**
   * sends a value to all processors
   */
  template <typename T>
  void broadcast(T* value, int root, std::optional<MPI_Comm> comm = {}) const {
    auto mpiType = castToMpiType<T>();
    MPI_Bcast(value, 1, mpiType, root, comm.value_or(this->comm()));
  }

  /**
   * sends a container to all processors
   */
  template <typename ContainerType>
  void broadcastContainer(ContainerType& container,
                          int root,
                          std::optional<MPI_Comm> comm = {}) const {
    using InternalType = ContainerType::value_type;
    const auto& useComm = comm.value_or(this->comm());
    auto size = static_cast<Count>(container.size());
    broadcast(&size, root, useComm);
    if (rankOf(useComm) != root) {
      container.resize(size);
    }
    auto mpiType = castToMpiType<InternalType>();
    SEISSOL_MPI_C(MPI_Bcast)(container.data(), size, mpiType, root, useComm);
  }

  /**
   * @return @p value reduced with @p op over all processes
   */
  template <typename T>
  [[nodiscard]] T allreduce(const T& value, MPI_Op op, std::optional<MPI_Comm> comm = {}) const {
    T result{};
    const auto mpiType = castToMpiType<T>();
    MPI_Allreduce(&value, &result, 1, mpiType, op, comm.value_or(this->comm()));
    return result;
  }

  /**
   * Reduces the elements of @p container with @p op over all processes, in place.
   */
  template <typename ContainerType>
  void allreduceContainer(ContainerType& container,
                          MPI_Op op,
                          std::optional<MPI_Comm> comm = {}) const {
    using InternalType = ContainerType::value_type;
    const auto size = static_cast<Count>(container.size());
    const auto mpiType = castToMpiType<InternalType>();
    SEISSOL_MPI_C(MPI_Allreduce)
    (MPI_IN_PLACE, container.data(), size, mpiType, op, comm.value_or(this->comm()));
  }

  /**
   * @return @p value reduced with @p op over all processes, on @p root; a value-initialized T on
   * the others
   */
  template <typename T>
  [[nodiscard]] T
      reduce(const T& value, MPI_Op op, int root, std::optional<MPI_Comm> comm = {}) const {
    T result{};
    const auto mpiType = castToMpiType<T>();
    MPI_Reduce(&value, &result, 1, mpiType, op, root, comm.value_or(this->comm()));
    return result;
  }

  /**
   * Reduces the elements of @p container with @p op over all processes, in place on @p root; the
   * other processes keep theirs.
   */
  template <typename ContainerType>
  void reduceContainer(ContainerType& container,
                       MPI_Op op,
                       int root,
                       std::optional<MPI_Comm> comm = {}) const {
    using InternalType = ContainerType::value_type;
    const auto& useComm = comm.value_or(this->comm());
    const auto size = static_cast<Count>(container.size());
    const auto mpiType = castToMpiType<InternalType>();
    if (rankOf(useComm) == root) {
      SEISSOL_MPI_C(MPI_Reduce)
      (MPI_IN_PLACE, container.data(), size, mpiType, op, root, useComm);
    } else {
      SEISSOL_MPI_C(MPI_Reduce)(container.data(), nullptr, size, mpiType, op, root, useComm);
    }
  }

  /**
   * @return the values of all processes, on all of them
   */
  template <typename T>
  [[nodiscard]] std::vector<T> allgather(const T& value, std::optional<MPI_Comm> comm = {}) const {
    const auto& useComm = comm.value_or(this->comm());
    std::vector<T> values(sizeOf(useComm));
    const auto mpiType = castToMpiType<T>();
    MPI_Allgather(&value, 1, mpiType, values.data(), 1, mpiType, useComm);
    return values;
  }

  /**
   * @return @p value reduced with @p op over the processes up to this one (inclusive), or over the
   * processes before this one (exclusive; a value-initialized T on the first process)
   */
  template <typename T>
  [[nodiscard]] T
      scan(const T& value, MPI_Op op, bool inclusive, std::optional<MPI_Comm> comm = {}) const {
    const auto& useComm = comm.value_or(this->comm());
    T result{};
    auto mpiType = castToMpiType<T>();
    if (inclusive) {
      MPI_Scan(&value, &result, 1, mpiType, op, useComm);
    } else {
      MPI_Exscan(&value, &result, 1, mpiType, op, useComm);
      // MPI_Exscan leaves the receive buffer of the first rank undefined
      if (rankOf(useComm) == 0) {
        result = T{};
      }
    }

    return result;
  }

  /**
   * Scans the elements of @p container with @p op, in place, as scan() does each of them.
   */
  template <typename ContainerType>
  void scanContainer(ContainerType& container,
                     MPI_Op op,
                     bool inclusive,
                     std::optional<MPI_Comm> comm = {}) const {
    using InternalType = ContainerType::value_type;
    const auto& useComm = comm.value_or(this->comm());
    const auto size = static_cast<Count>(container.size());
    const auto mpiType = castToMpiType<InternalType>();
    if (inclusive) {
      SEISSOL_MPI_C(MPI_Scan)(MPI_IN_PLACE, container.data(), size, mpiType, op, useComm);
    } else {
      SEISSOL_MPI_C(MPI_Exscan)(MPI_IN_PLACE, container.data(), size, mpiType, op, useComm);
      // MPI_Exscan leaves the receive buffer of the first rank undefined
      if (rankOf(useComm) == 0) {
        std::fill(container.begin(), container.end(), InternalType{});
      }
    }
  }

  /**
   * @return The main communicator for the application; on a thread that has one of its own (see
   * setThreadComm), that one.
   */
  [[nodiscard]] MPI_Comm comm() const {
    return threadCommRef() != MPI_COMM_NULL ? threadCommRef() : comm_;
  }

  /**
   * @return The communicator of the outputs, a duplicate of the main one: the thread of the outputs
   * (see OutputWorker) talks over it, so that its collectives cannot mix with those of the main
   * thread.
   */
  [[nodiscard]] MPI_Comm outputComm() const { return outputComm_; }

  /**
   * Makes comm() return `comm` on the calling thread; MPI_COMM_NULL restores the main communicator.
   */
  static void setThreadComm(MPI_Comm comm) { threadCommRef() = comm; }

  /**
   * @return The node communicator (shared memory) for the application
   */
  [[nodiscard]] MPI_Comm sharedMemComm() const { return sharedMemComm_; }

  void barrier(std::optional<MPI_Comm> comm = {}) const {
    MPI_Barrier(comm.value_or(this->comm()));
  }

  /**
   * Finalize MPI
   */
  void finalize();

  /** The only instance of the class */
  static Mpi mpi;

  private:
  static MPI_Comm& threadCommRef() {
    thread_local MPI_Comm threadComm = MPI_COMM_NULL;
    return threadComm;
  }

  MPI_Comm comm_{MPI_COMM_NULL};
  MPI_Comm outputComm_{MPI_COMM_NULL};
  MPI_Comm sharedMemComm_{MPI_COMM_NULL};
  Mpi() = default;
};

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_MPI_H_
