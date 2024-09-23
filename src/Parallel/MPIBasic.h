// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_PARALLEL_MPIBASIC_H_
#define SEISSOL_SRC_PARALLEL_MPIBASIC_H_

namespace seissol {

/**
 * Basic MPI abstraction
 */
class MPIBasic {
  protected:
  /** This rank */
  int m_rank;

  /** Rank in the shared memory sub-communicator */
  int m_sharedMemMpiRank;

  /** Number of processors */
  int m_size;

  /** Number of ranks in the shared memory sub-communicator */
  int m_sharedMemMpiSize;

  /** Requires threadsafe MPI */
  bool m_threadsafe;

  MPIBasic() : m_rank(0), m_size(1) {}

  public:
  virtual ~MPIBasic() = default;

  /**
   * @return The rank of this process
   */
  int rank() const { return m_rank; }

  /**
   * @return The rank within the shared memory sub-communicator
   */
  int sharedMemMpiRank() const { return m_sharedMemMpiRank; }

  /**
   * @return The total number of processes
   */
  int size() const { return m_size; }

  /**
   * @return The number of ranks within the shared memory sub-communicator
   */
  int sharedMemMpiSize() const { return m_sharedMemMpiSize; }

  bool isSingleProcess() const { return size() == 1; }

  bool isSingleNode() const { return size() == sharedMemMpiSize(); }
};
} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_MPIBASIC_H_
