// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_DEDUPLICATE_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_DEDUPLICATE_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace seissol::io::instance::geometry {

/**
 * @brief The unique points of a point list, and where each input point went.
 */
struct PointMap {
  //! The unique points, in the order in which they first occur.
  std::vector<double> points;
  //! For every input point, the index of the unique point it became.
  std::vector<std::size_t> indices;

  [[nodiscard]] std::size_t pointCount() const { return points.size() / 3; }
};

/**
 * @brief Merges points that coincide, keeping the order of first occurrence.
 *
 * The output of a mesh writer repeats every vertex once per cell that touches it, which for
 * tetrahedra is around twenty times. Merging them locally -- within one rank -- shrinks the point
 * array accordingly without any communication: a VTKHDF file is read partition by partition, so
 * points that two ranks share may stay duplicated.
 *
 * Coordinates are snapped to a grid to find candidates and then compared with @p tolerance, which
 * is relative to the extent of the input. Snapping alone would separate two points that fall on
 * either side of a cell boundary, so the neighbouring cells are searched as well; comparing
 * afterwards means the grid only has to be fine enough to keep distinct points apart.
 */
inline PointMap deduplicatePoints(const std::vector<double>& coordinates,
                                  double tolerance = 1e-12) {
  PointMap result;
  const auto count = coordinates.size() / 3;
  result.indices.resize(count);
  if (count == 0) {
    return result;
  }

  std::array<double, 3> minimum{coordinates[0], coordinates[1], coordinates[2]};
  std::array<double, 3> maximum = minimum;
  for (std::size_t point = 1; point < count; ++point) {
    for (std::size_t d = 0; d < 3; ++d) {
      minimum[d] = std::min(minimum[d], coordinates[point * 3 + d]);
      maximum[d] = std::max(maximum[d], coordinates[point * 3 + d]);
    }
  }
  double extent = 0;
  for (std::size_t d = 0; d < 3; ++d) {
    extent = std::max(extent, maximum[d] - minimum[d]);
  }
  const auto epsilon = extent > 0 ? extent * tolerance : tolerance;
  // coarse enough that coinciding points land in the same cell or a neighbouring one
  const auto spacing = epsilon > 0 ? epsilon * 16 : 1.0;

  struct Key {
    std::int64_t x, y, z;
    bool operator==(const Key& other) const { return x == other.x && y == other.y && z == other.z; }
  };
  struct KeyHash {
    std::size_t operator()(const Key& key) const {
      // three odd primes, so that neighbouring cells do not collide systematically
      return static_cast<std::size_t>(key.x * 73856093L ^ key.y * 19349663L ^ key.z * 83492791L);
    }
  };

  std::unordered_map<Key, std::vector<std::size_t>, KeyHash> buckets;
  buckets.reserve(count);

  const auto cellOf = [&](const double* point) {
    return Key{static_cast<std::int64_t>(std::floor(point[0] / spacing)),
               static_cast<std::int64_t>(std::floor(point[1] / spacing)),
               static_cast<std::int64_t>(std::floor(point[2] / spacing))};
  };

  for (std::size_t point = 0; point < count; ++point) {
    const auto* current = &coordinates[point * 3];
    const auto cell = cellOf(current);

    std::size_t found = result.pointCount();
    for (std::int64_t dx = -1; dx <= 1 && found == result.pointCount(); ++dx) {
      for (std::int64_t dy = -1; dy <= 1 && found == result.pointCount(); ++dy) {
        for (std::int64_t dz = -1; dz <= 1 && found == result.pointCount(); ++dz) {
          const auto neighbour = buckets.find(Key{cell.x + dx, cell.y + dy, cell.z + dz});
          if (neighbour == buckets.end()) {
            continue;
          }
          for (const auto candidate : neighbour->second) {
            const auto* other = &result.points[candidate * 3];
            if (std::abs(current[0] - other[0]) <= epsilon &&
                std::abs(current[1] - other[1]) <= epsilon &&
                std::abs(current[2] - other[2]) <= epsilon) {
              found = candidate;
              break;
            }
          }
        }
      }
    }

    if (found == result.pointCount()) {
      buckets[cell].push_back(found);
      result.points.insert(result.points.end(), current, current + 3);
    }
    result.indices[point] = found;
  }

  return result;
}

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_DEDUPLICATE_H_
