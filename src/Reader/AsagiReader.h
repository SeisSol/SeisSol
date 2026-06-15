// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_READER_ASAGIREADER_H_
#define SEISSOL_SRC_READER_ASAGIREADER_H_

#ifdef USE_ASAGI

#include "Parallel/MPI.h"
#include "easi/util/AsagiReader.h"

namespace asagi {
class Grid;
} // namespace asagi

namespace seissol::asagi {
enum class NumaCacheMode { Off, On, Cache };

class AsagiReader : public easi::AsagiReader {
  private:
  /** Number of threads used by ASAGI */
  unsigned int asagiThreads_{0};

  /** MPI communicator used by ASAGI */
  MPI_Comm comm_;

  public:
  explicit AsagiReader(MPI_Comm comm = seissol::Mpi::mpi.comm());

  ::asagi::Grid* open(const char* file, const char* varname) override;
  [[nodiscard]] unsigned numberOfThreads() const override;

  private:
  static NumaCacheMode getNumaMode();
};

} // namespace seissol::asagi

#endif

#endif // SEISSOL_SRC_READER_ASAGIREADER_H_
