// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Output/DataTypes.h"

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::dr::output;

// ---------------------------------------------------------------------------
// VarT: basic properties
// ---------------------------------------------------------------------------

TEST_CASE("VarT dim" * doctest::test_suite("dynamicrupture")) {
  Var1D v1;
  CHECK(v1.dim() == 1);
  Var2D v2;
  CHECK(v2.dim() == 2);
  Var3D v3;
  CHECK(v3.dim() == 3);
}

TEST_CASE("VarT default inactive" * doctest::test_suite("dynamicrupture")) {
  Var1D v;
  CHECK(v.isActive == false);
  CHECK(v.size == 0);
  CHECK(v.maxCacheLevel == 1);
}

// ---------------------------------------------------------------------------
// VarT: allocateData when active
// ---------------------------------------------------------------------------

TEST_CASE("VarT allocateData active" * doctest::test_suite("dynamicrupture")) {
  Var2D v;
  v.isActive = true;
  v.allocateData(100);

  CHECK(v.size == 100);
  CHECK(v.data[0].size() == 100); // 100 * maxCacheLevel(1)
  CHECK(v.data[1].size() == 100);
}

TEST_CASE("VarT allocateData inactive" * doctest::test_suite("dynamicrupture")) {
  Var2D v;
  v.isActive = false;
  v.allocateData(100);

  CHECK(v.size == 100);
  CHECK(v.data[0].empty());
  CHECK(v.data[1].empty());
}

// ---------------------------------------------------------------------------
// VarT: access operators
// ---------------------------------------------------------------------------

TEST_CASE("VarT 1D access" * doctest::test_suite("dynamicrupture")) {
  Var1D v;
  v.isActive = true;
  v.allocateData(10);

  // Write via operator(dim, level, index)
  v(0, 0, 5) = 42.0;
  CHECK(v(0, 0, 5) == doctest::Approx(42.0));

  // Write via operator(level, index) - 1D shorthand
  v(0, 3) = 99.0;
  CHECK(v(0, 3) == doctest::Approx(99.0));
}

TEST_CASE("VarT 3D access" * doctest::test_suite("dynamicrupture")) {
  Var3D v;
  v.isActive = true;
  v.allocateData(20);

  v(0, 0, 7) = 1.0;
  v(1, 0, 7) = 2.0;
  v(2, 0, 7) = 3.0;

  CHECK(v(0, 0, 7) == doctest::Approx(1.0));
  CHECK(v(1, 0, 7) == doctest::Approx(2.0));
  CHECK(v(2, 0, 7) == doctest::Approx(3.0));
}

TEST_CASE("VarT bracket operator" * doctest::test_suite("dynamicrupture")) {
  Var2D v;
  v.isActive = true;
  v.allocateData(5);

  real* ptr0 = v[0];
  real* ptr1 = v[1];
  CHECK(ptr0 != nullptr);
  CHECK(ptr1 != nullptr);
  CHECK(ptr0 != ptr1);
}

// ---------------------------------------------------------------------------
// VarT: resizeCache
// ---------------------------------------------------------------------------

TEST_CASE("VarT resizeCache" * doctest::test_suite("dynamicrupture")) {
  Var1D v;
  v.isActive = true;
  v.allocateData(10); // 10 * 1 = 10 elements
  CHECK(v.data[0].size() == 10);

  v.resizeCache(5);
  CHECK(v.maxCacheLevel == 5);
  CHECK(v.data[0].size() == 50); // 10 * 5

  v.resizeCache(1);
  CHECK(v.data[0].size() == 10); // 10 * 1
}

TEST_CASE("VarT resizeCache inactive no-op" * doctest::test_suite("dynamicrupture")) {
  Var1D v;
  v.isActive = false;
  v.allocateData(10);
  v.resizeCache(5);
  CHECK(v.data[0].empty());
}

// ---------------------------------------------------------------------------
// VarT: multi-level caching access
// ---------------------------------------------------------------------------

TEST_CASE("VarT cache level access" * doctest::test_suite("dynamicrupture")) {
  Var1D v;
  v.isActive = true;
  v.allocateData(4);
  v.resizeCache(3); // 3 cache levels, 4 points each → 12 elements

  // Write different values at different cache levels
  v(0, 0) = 10.0; // level 0, index 0
  v(1, 0) = 20.0; // level 1, index 0
  v(2, 0) = 30.0; // level 2, index 0

  CHECK(v(0, 0) == doctest::Approx(10.0));
  CHECK(v(1, 0) == doctest::Approx(20.0));
  CHECK(v(2, 0) == doctest::Approx(30.0));

  // Different index at same level
  v(0, 3) = 99.0;
  CHECK(v(0, 3) == doctest::Approx(99.0));
  CHECK(v(0, 0) == doctest::Approx(10.0)); // unchanged
}

// ---------------------------------------------------------------------------
// VariableLabels consistency
// ---------------------------------------------------------------------------

TEST_CASE("VariableLabels matches VariableID::Size" * doctest::test_suite("dynamicrupture")) {
  CHECK(VariableLabels.size() == static_cast<std::size_t>(VariableID::Size));
}

TEST_CASE("VariableLabels all non-empty" * doctest::test_suite("dynamicrupture")) {
  for (std::size_t i = 0; i < VariableLabels.size(); ++i) {
    CHECK_FALSE(VariableLabels[i].empty());
    for (const auto& label : VariableLabels[i]) {
      CHECK_FALSE(label.empty());
    }
  }
}

TEST_CASE("VariableLabels specific checks" * doctest::test_suite("dynamicrupture")) {
  CHECK(VariableLabels[VariableID::SlipRate].size() == 2);           // strike, dip
  CHECK(VariableLabels[VariableID::TransientTractions].size() == 3); // Ts, Td, Pn
  CHECK(VariableLabels[VariableID::NormalVelocity].size() == 1);
  CHECK(VariableLabels[VariableID::Slip].size() == 2);
  CHECK(VariableLabels[VariableID::RuptureTime].size() == 1);
}

// ---------------------------------------------------------------------------
// DrVarsT tuple size matches VariableID::Size
// ---------------------------------------------------------------------------

TEST_CASE("DrVarsT tuple matches VariableID count" * doctest::test_suite("dynamicrupture")) {
  CHECK(std::tuple_size_v<DrVarsT> == static_cast<std::size_t>(VariableID::Size));
}

// ---------------------------------------------------------------------------
// DirectionID and ParamID enums
// ---------------------------------------------------------------------------

TEST_CASE("DirectionID values" * doctest::test_suite("dynamicrupture")) {
  CHECK(static_cast<int>(DirectionID::Strike) == 0);
  CHECK(static_cast<int>(DirectionID::Dip) == 1);
  CHECK(static_cast<int>(DirectionID::Normal) == 2);
}

TEST_CASE("ParamID values" * doctest::test_suite("dynamicrupture")) {
  CHECK(static_cast<int>(ParamID::FrictionCoefficient) == 0);
  CHECK(static_cast<int>(ParamID::State) == 1);
}

} // namespace seissol::unit_test
