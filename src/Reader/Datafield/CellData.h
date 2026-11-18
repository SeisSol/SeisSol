// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_READER_DATAFIELD_CELLDATA_H_
#define SEISSOL_SRC_READER_DATAFIELD_CELLDATA_H_

#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <map>
#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::reader::datafield {

/// Data on the mesh, of the datasets of the kind cells (GridKind::Cells): a row of values per cell
/// of the mesh, in the order of the cells in the mesh file, in an HDF5 dataset whose first axis
/// they are -- the mesh file itself or another one. A consumer reads the rows of its cells by their
/// ids in the file (Element::globalId, the gid of a cell before the partition), collectively, and
/// binds the value of the cell of every point to the input channel `d.c` of a component.
///
/// The store keeps, per dataset (by its file and variable), the rows read last, of the cells of
/// this rank: a later read of the same dataset replaces them, so that the cells before the
/// partition give way to those after it.
class CellDataStore {
  public:
  /// Reads the rows of the cells `ids` of the dataset `desc` -- every rank its own cells, the ghost
  /// cells too --, collectively over `comm`: every rank reads a part of the rows of the file, which
  /// go to the ranks of their cells. An error where the file or its variable is missing, where a
  /// row has another number of values than the dataset has components, and where the dataset has
  /// another number of rows than the mesh has cells (the greatest id of a cell, plus one).
  void load(const GridDesc& desc, const std::vector<std::size_t>& ids, MPI_Comm comm);

  /// Whether the rows of `desc` were read.
  [[nodiscard]] bool loaded(const GridDesc& desc) const;

  /// The values of the cell `id` of the dataset `desc`, one per component; an error where its rows
  /// were not read, or not those of the cell.
  [[nodiscard]] const double* row(const GridDesc& desc, std::size_t id) const;

  /// Forgets every dataset.
  void clear() { datasets_.clear(); }

  private:
  struct Rows {
    std::size_t width{0};
    std::unordered_map<std::size_t, std::size_t> position;
    std::vector<double> values;
  };
  [[nodiscard]] static std::string keyOf(const GridDesc& desc);

  std::map<std::string, Rows> datasets_;
};

/// The process-wide store, which the initialization fills before it evaluates the materials.
[[nodiscard]] CellDataStore& sharedCellDataStore();

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_CELLDATA_H_
