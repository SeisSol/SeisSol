// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_READER_DATAFIELD_ASAGIREADER_H_
#define SEISSOL_SRC_READER_DATAFIELD_ASAGIREADER_H_

#include "Parallel/MPI.h"

namespace asagi {
class Grid;
} // namespace asagi

#ifdef USE_EASI
#include "easi/util/AsagiReader.h"
#else
// class interface replacement if not linked against easi
namespace easi {
class AsagiReader {
  public:
  virtual ~AsagiReader() = default;
  virtual ::asagi::Grid* open(const char* file, const char* varname) = 0;
  [[nodiscard]] virtual unsigned numberOfThreads() const { return 1; }
};
} // namespace easi
#endif

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

#endif // SEISSOL_SRC_READER_DATAFIELD_ASAGIREADER_H_
