// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Core/Data/DataTable.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::unit_test {

// a struct with virtual functions, as the materials are: not standard-layout
struct DataTablePolymorphic {
  virtual ~DataTablePolymorphic() = default;
  [[nodiscard]] virtual double sum() const { return a + b; }
  double a{0.0};
  double b{0.0};
};

TEST_CASE("DataTable: a member view writes the member of a struct with virtual functions" *
          doctest::test_suite("reader")) {
  std::vector<DataTablePolymorphic> structs(3);
  data::DataTable table(structs.size());
  table.bindMemberView("b", data::Direction::Out, structs.data(), &DataTablePolymorphic::b);

  const auto& entry = table.dataEntries().at(0);
  for (std::size_t i = 0; i < structs.size(); ++i) {
    entry.setValueAs(i, 1.5 * static_cast<double>(i));
  }
  for (std::size_t i = 0; i < structs.size(); ++i) {
    CHECK(structs[i].a == 0.0);
    CHECK(structs[i].b == 1.5 * static_cast<double>(i));
  }

  // the view says the same as the setter
  REQUIRE(entry.view.has_value());
  const auto& view = entry.view.value();
  const auto* bytes = static_cast<const char*>(view.base) + view.byteOffset;
  CHECK(reinterpret_cast<const double*>(bytes + view.displacement(2)) == &structs[2].b);
}

TEST_CASE("DataTable: a cell point view with an index writes the points of the cells it names" *
          doctest::test_suite("reader")) {
  // 4 cells of 2 points with 3 values each (e.g. fused simulations); the table holds the cells 3
  // and 1, in that order, and writes the value 1 of their points
  constexpr std::size_t PointsPerCell = 2;
  constexpr std::size_t Values = 3;
  std::vector<float> memory(4 * PointsPerCell * Values, -1.0F);
  const std::array<std::uint32_t, 2> cells{3, 1};

  data::DataTable table(cells.size() * PointsPerCell);
  table.bindCellPointView("v",
                          data::Direction::Out,
                          memory.data(),
                          PointsPerCell,
                          PointsPerCell * Values,
                          Values,
                          1,
                          cells.data());
  const auto& entry = table.dataEntries().at(0);
  for (std::size_t point = 0; point < table.numPoints(); ++point) {
    entry.setValueAs(point, static_cast<double>(point));
  }

  for (std::size_t cell = 0; cell < 4; ++cell) {
    for (std::size_t point = 0; point < PointsPerCell; ++point) {
      for (std::size_t value = 0; value < Values; ++value) {
        CAPTURE(cell);
        CAPTURE(point);
        CAPTURE(value);
        const auto written = memory[(cell * PointsPerCell + point) * Values + value];
        if (value == 1 && cell == 3) {
          CHECK(written == static_cast<float>(point));
        } else if (value == 1 && cell == 1) {
          CHECK(written == static_cast<float>(PointsPerCell + point));
        } else {
          CHECK(written == -1.0F);
        }
      }
    }
  }
}

TEST_CASE("DataTable: a column computed a range at a time" * doctest::test_suite("reader")) {
  std::size_t calls = 0;
  data::DataTable table(10);
  table.bindComputedBatch<double>("x", [&calls](std::size_t first, std::size_t count, double* out) {
    ++calls;
    for (std::size_t i = 0; i < count; ++i) {
      out[i] = 2.0 * static_cast<double>(first + i);
    }
  });
  const auto& entry = table.dataEntries().at(0);

  std::array<double, 4> values{};
  entry.getValues<double>(3, values.size(), values.data());
  CHECK(calls == 1);
  for (std::size_t i = 0; i < values.size(); ++i) {
    CHECK(values[i] == 2.0 * static_cast<double>(3 + i));
  }
  CHECK(entry.getValue<double>(7) == 14.0);
}

} // namespace seissol::unit_test
