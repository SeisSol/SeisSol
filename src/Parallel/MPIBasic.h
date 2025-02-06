// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_PARALLEL_MPIBASIC_H_
#define SEISSOL_SRC_PARALLEL_MPIBASIC_H_

namespace seissol {

/**
 * Basic MPI abstraction
 */
class MPIBasic {
  protected:
  /** This rank */
  int m_rank{0};

  /** Rank in the shared memory sub-communicator */
  int m_sharedMemMpiRank{0};

  /** Number of processors */
  int m_size{1};

  /** Number of ranks in the shared memory sub-communicator */
  int m_sharedMemMpiSize{1};

  MPIBasic() = default;

  public:
  virtual ~MPIBasic() = default;

  /**
   * @return The rank of this process
   */
  [[nodiscard]] int rank() const { return m_rank; }

  /**
   * @return The rank within the shared memory sub-communicator
   */
  [[nodiscard]] int sharedMemMpiRank() const { return m_sharedMemMpiRank; }

  /**
   * @return The total number of processes
   */
  [[nodiscard]] int size() const { return m_size; }

  /**
   * @return The number of ranks within the shared memory sub-communicator
   */
  [[nodiscard]] int sharedMemMpiSize() const { return m_sharedMemMpiSize; }

  [[nodiscard]] bool isSingleProcess() const { return size() == 1; }

  [[nodiscard]] bool isSingleNode() const { return size() == sharedMemMpiSize(); }
};
} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_MPIBASIC_H_
