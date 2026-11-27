// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Core/Data/DataTable.h"
#include "Core/IO/Easi/EasiContext.h"
#include "Core/IO/Easi/EasiReader.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("EasiReader: a model may supply more than the table asks for" *
          doctest::test_suite("reader")) {
  // the material and the plasticity in one file, of which a table asks for the material only; in
  // components that every build of easi has (a FunctionMap needs Lua or ImpalaJIT): the
  // coefficients of a PolynomialMap go from the highest power down
  const auto path = std::filesystem::temp_directory_path() / "seissol-easireader-more.yaml";
  {
    std::ofstream file(path);
    file << "!AffineMap\n"
            "matrix:\n"
            "  xi: [1.0, 0.0, 0.0]\n"
            "translation:\n"
            "  xi: 0.0\n"
            "components: !PolynomialMap\n"
            "  map:\n"
            "    rho: [1.0, 2700.0]\n"
            "    mu: [0.0, 30000000000.0]\n"
            "    plastCo: [0.0, 1000000.0]\n"
            "    bulkFriction: [0.0, 0.6]\n";
  }

  data::EasiReader reader(path.string(), {"x", "y", "z"}, data::EasiContext{});
  CHECK(reader.outputVars().size() == 4);

  constexpr std::size_t Points = 3;
  const std::array<std::array<double, 3>, Points> points{
      {{0.0, 1.0, 2.0}, {1.5, 0.0, 0.0}, {-2.0, 3.0, 1.0}}};
  std::vector<double> rho(Points, -1.0);
  std::vector<double> mu(Points, -1.0);
  data::DataTable table(Points);
  using data::Direction;
  table.bindViewConst("x", Direction::In, points.data()->data(), 3, 0);
  table.bindViewConst("y", Direction::In, points.data()->data(), 3, 1);
  table.bindViewConst("z", Direction::In, points.data()->data(), 3, 2);
  table.bindView("rho", Direction::Out, rho.data());
  table.bindView("mu", Direction::Out, mu.data());
  reader.prepare(table);
  reader.call(table);

  for (std::size_t p = 0; p < Points; ++p) {
    CAPTURE(p);
    CHECK(rho[p] == doctest::Approx(2700.0 + points[p][0]));
    CHECK(mu[p] == doctest::Approx(3.0e10));
  }
  std::filesystem::remove(path);
}

} // namespace seissol::unit_test
