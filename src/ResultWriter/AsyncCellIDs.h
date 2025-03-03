// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_ASYNCCELLIDS_H_
#define SEISSOL_SRC_RESULTWRITER_ASYNCCELLIDS_H_

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "SeisSol.h"

namespace seissol {

/**
 * This class can fix cells (vertex ids) in asynchronous mode.
 *
 * Sicne cells assume local vertex ids, we have to add an additional
 * when using the asynchronous MPI mode.
 *
 * @tparam CellVertices Number of vertices per cell
 */
template <int CellVertices>
class AsyncCellIDs {
  private:
  /** Null, if MPI is not enabled */
  std::vector<unsigned> localCells;

  const unsigned int* constCells;

  public:
  AsyncCellIDs(unsigned int nCells,
               unsigned int nVertices,
               const unsigned int* cells,
               seissol::SeisSol& seissolInstance) {
#ifdef USE_MPI
    // Add the offset to the cells
    MPI_Comm groupComm = seissolInstance.asyncIO().groupComm();
    unsigned int offset = nVertices;
    MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
    offset -= nVertices;

    // Add the offset to all cells
    localCells.resize(nCells * CellVertices);
    for (unsigned int i = 0; i < nCells * CellVertices; i++) {
      localCells[i] = cells[i] + offset;
    }
    constCells = localCells.data();
#else  // USE_MPI
    constCells = cells;
#endif // USE_MPI
  }

  [[nodiscard]] const unsigned int* cells() const { return constCells; }
};

} // namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_ASYNCCELLIDS_H_
