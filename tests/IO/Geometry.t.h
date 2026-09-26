// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Instance/Geometry/Geometry.h"

namespace seissol::unit_test {

TEST_CASE("IO/Geometry: a format that cannot hold a time series says so" *
          doctest::test_suite("io")) {
  using namespace seissol::io::instance::geometry;

  SUBCASE("VTKHDF keeps what it was asked for") {
    CHECK(supportedWriterGroup(WriterGroup::Monolith, WriterFormat::Vtk, "wavefield") ==
          WriterGroup::Monolith);
    CHECK(supportedWriterGroup(WriterGroup::IncrementalSnapshot, WriterFormat::Vtk, "wavefield") ==
          WriterGroup::IncrementalSnapshot);
    CHECK(supportedWriterGroup(WriterGroup::FullSnapshot, WriterFormat::Vtk, "wavefield") ==
          WriterGroup::FullSnapshot);
  }

  SUBCASE("Xdmf falls back to a snapshot per step") {
    CHECK(supportedWriterGroup(WriterGroup::Monolith, WriterFormat::Xdmf, "wavefield") ==
          WriterGroup::FullSnapshot);
    CHECK(supportedWriterGroup(WriterGroup::IncrementalSnapshot, WriterFormat::Xdmf, "wavefield") ==
          WriterGroup::FullSnapshot);
    CHECK(supportedWriterGroup(WriterGroup::FullSnapshot, WriterFormat::Xdmf, "wavefield") ==
          WriterGroup::FullSnapshot);
  }
}

} // namespace seissol::unit_test
