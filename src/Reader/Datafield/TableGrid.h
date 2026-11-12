// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_TABLEGRID_H_
#define SEISSOL_SRC_READER_DATAFIELD_TABLEGRID_H_

// A grid of a text table (GridKind::Table). The text: a first line with the names of the columns,
// then a row per line, its fields separated by spaces or commas; `#` starts a comment, and a line
// without fields is none. The numbers are spelled as sderiv spells them -- decimal or hexadecimal
// floating point, nan and inf, each with a sign -- and rounded correctly. The points of `sderiv
// eval` are such a text, too.
//
// The columns a declaration names as its coordinates are the axes of the grid, axis 0 first: the
// knots of an axis are the distinct values of its column, in increasing order, finite. Every row is
// a node, and every node of the grid the knots span is given once, in any order. The columns it
// names as its components are the components of the grid, in its order; the other columns are not
// read. The table becomes the data of an inline grid (GridData), and is sampled as one: every axis
// is knots.

#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <string>
#include <vector>

namespace seissol::reader::datafield {

/// The columns of a text table by their names.
struct TextTable {
  std::vector<std::string> names;
  /// [column][row]
  std::vector<std::vector<double>> columns;
  /// the line of the text of every row
  std::vector<std::size_t> lines;

  [[nodiscard]] std::size_t rows() const { return lines.size(); }
  /// The column `name`, or null.
  [[nodiscard]] const std::vector<double>* column(const std::string& name) const;
};

/// The table of the text `text`, of the file `path` that a message names. Throws
/// std::runtime_error, with the line, for a text that is no table.
[[nodiscard]] TextTable parseTextTable(const std::string& text, const std::string& path);

/// The table of the file `path`. Throws std::runtime_error where it cannot be read or is no table.
[[nodiscard]] TextTable readTextTable(const std::string& path);

/// The data of the grid of `table` (of the file `path`) whose axes are the columns `coordinates`
/// and whose components are the columns `components`. Throws std::runtime_error where a column is
/// missing, a coordinate is not finite, or the rows are not the nodes of a grid, each once.
[[nodiscard]] GridData tableGridData(const TextTable& table,
                                     const std::vector<std::string>& coordinates,
                                     const std::vector<std::string>& components,
                                     const std::string& path);

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_TABLEGRID_H_
