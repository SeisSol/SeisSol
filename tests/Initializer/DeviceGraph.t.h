// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/BasicTypedefs.h"
#include "Initializer/DeviceGraph.h"

namespace seissol::unit_test {

TEST_CASE("Device graph keys hash all their fields" * doctest::test_suite("initializer")) {
  const initializer::GraphKeyHash hash{};
  const initializer::GraphKey key(ComputeGraphType::AccumulatedVelocities, 0.5, false);
  const initializer::GraphKey sameKey(ComputeGraphType::AccumulatedVelocities, 0.5, false);
  const initializer::GraphKey keyWithDisplacements(
      ComputeGraphType::AccumulatedVelocities, 0.5, true);

  CHECK(key == sameKey);
  CHECK(hash(key) == hash(sameKey));

  // hashes may collide in general; but std::hash<bool> separates false and true,
  // hence these two keys only collide if the hash ignores withDisplacements
  CHECK_FALSE(key == keyWithDisplacements);
  CHECK(hash(key) != hash(keyWithDisplacements));
}

} // namespace seissol::unit_test
