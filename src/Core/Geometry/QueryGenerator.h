// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_CORE_GEOMETRY_QUERYGENERATOR_H_
#define SEISSOL_SRC_CORE_GEOMETRY_QUERYGENERATOR_H_

#include "Core/Common/Constants.h"
#include "Core/Data/DataTable.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace seissol::initializer {

/// The points at which a model is queried: the columns x, y, z and group of a table.
class QueryGenerator {
  public:
  virtual ~QueryGenerator() = default;
  [[nodiscard]] virtual data::DataTable generate() const = 0;
  [[nodiscard]] virtual std::size_t outputPerCell() const { return 1; }
};

/// Binds the point set of a query, computed a range of points at a time: `coordinates(first,
/// count, xyz)` writes the coordinates of the points [first, first + count) to xyz[i * 3 + d], and
/// `groups(first, count, out)` their groups. Nothing is materialised; the table refers to whatever
/// the two callbacks capture.
template <typename CoordinateFn, typename GroupFn>
void bindPointSet(data::DataTable& table, CoordinateFn coordinates, GroupFn groups) {
  const auto shared = std::make_shared<CoordinateFn>(std::move(coordinates));
  for (std::size_t d = 0; d < Cell::Dim; ++d) {
    table.bindComputedBatch<double>(std::string(1, "xyz"[d]),
                                    [shared, d](std::size_t first, std::size_t count, double* out) {
                                      thread_local std::vector<double> xyz;
                                      xyz.resize(count * Cell::Dim);
                                      (*shared)(first, count, xyz.data());
                                      for (std::size_t i = 0; i < count; ++i) {
                                        out[i] = xyz[i * Cell::Dim + d];
                                      }
                                    });
  }
  table.bindComputedBatch<std::int32_t>("group", std::move(groups));
}

} // namespace seissol::initializer

#endif // SEISSOL_SRC_CORE_GEOMETRY_QUERYGENERATOR_H_
