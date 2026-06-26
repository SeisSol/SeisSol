// SPDX-FileCopyrightText: 2015 SeisSol Group
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
class MpiBasic {
  protected:
  /** This rank */
  int rank_{0};

  /** Rank in the shared memory sub-communicator */
  int sharedMemMpiRank_{0};

  /** Number of processors */
  int size_{1};

  /** Number of ranks in the shared memory sub-communicator */
  int sharedMemMpiSize_{1};

  MpiBasic() = default;

  public:
  virtual ~MpiBasic() = default;

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
};
} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_MPIBASIC_H_
