// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_POINT_GROUPING_H_
#define SEISSOL_SRC_IO_INSTANCE_POINT_GROUPING_H_

#include "IO/Instance/Point/TableWriter.h"

#include <cstddef>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::io::instance::point {

/**
 * @brief How a set of points is laid out when not all of them record the same quantities.
 *
 * Which quantities a point has follows from the material of the element it sits in, so it is
 * fixed at setup and never changes during a run. Points are therefore gathered into groups that
 * share a quantity set, and every group becomes a dense table of its own: no filler values that a
 * reader cannot tell from a measurement, and no per-sample price for something that is static.
 *
 * Inside a group the points are renumbered so that every rank owns one contiguous range. That
 * keeps the layout expressible as a single distributed dimension -- the alternative, letting the
 * ranks own scattered indices, would need a per-point index in the write plan. The renumbering is
 * a permutation of the input order, so the mapping has to be written alongside the data.
 */
struct Grouping {
  //! The quantity set of every group. Known on every rank, including for groups it has no points
  //! of, since a rank still takes part in declaring the dataset.
  std::vector<std::vector<TableQuantity>> quantities;
  //! For every local point, the group it belongs to.
  std::vector<std::size_t> group;
  //! For every local point, its index inside its group: globally unique, contiguous per rank.
  std::vector<std::size_t> index;
  //! For every group, how many points it holds over all ranks.
  std::vector<std::size_t> globalCount;

  [[nodiscard]] std::size_t groupCount() const { return quantities.size(); }
};

/**
 * @brief Groups the local points by their quantity set and renumbers them inside each group.
 *
 * Collective on @p comm : the groups and their order have to be the same everywhere, so the
 * quantity sets are exchanged rather than derived from what a rank happens to hold.
 */
Grouping groupPoints(const std::vector<std::vector<TableQuantity>>& pointQuantities, MPI_Comm comm);

//! @brief The text a quantity set is identified by; also what it is reconstructed from.
std::string quantitySetKey(const std::vector<TableQuantity>& quantities);
std::vector<TableQuantity> quantitySetFromKey(const std::string& key);

} // namespace seissol::io::instance::point

#endif // SEISSOL_SRC_IO_INSTANCE_POINT_GROUPING_H_
