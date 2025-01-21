// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_PARALLEL_MPIDUMMY_H_
#define SEISSOL_SRC_PARALLEL_MPIDUMMY_H_

#include "Parallel/MPIBasic.h"

namespace seissol {

// MPI_Comm substitute, so we don't have to disambiguate when it's not there
#ifndef USE_MPI
using MPI_Comm = int;
constexpr MPI_Comm MPI_COMM_NULL = 0;
constexpr MPI_Comm MPI_COMM_SELF = 1;
constexpr MPI_Comm MPI_COMM_WORLD = 2;
#endif

/**
 * Dummy class when running without MPI
 */
class MPIDummy : public MPIBasic {
  private:
  MPIDummy() = default;

  public:
  ~MPIDummy() override = default;

  /**
   * Does nothing
   */
  void init(int& argc, char**& argv) {}

  /**
   * @return Dummy communicator
   */
  [[nodiscard]] int comm() const { return 0; }

  /**
   * Does nothing
   */
  void barrier(int comm) const {}

  /**
   * Does nothing
   */
  void finalize() {}

  /** The only instance of the class */
  static MPIDummy mpi;
};

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_MPIDUMMY_H_
