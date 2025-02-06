// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ResultWriter/ReceiverWriter.h"
#include <string>

namespace seissol::unit_test {

TEST_CASE("Parses line correctly") {
  const auto lines = std::vector<std::string>{
      "1\t0.1  10",
      "10\t2\t0.2",
      "1e4   1e3  1e-10",
      "-1e4   1e3  1e-10",
      "1e4\t\t-1e3 \t 1e-10",
  };
  const auto expectedPoints = std::vector<Eigen::Vector3d>{
      {1, 0.1, 10}, {10, 2, 0.2}, {1e4, 1e3, 1e-10}, {-1e4, 1e3, 1e-10}, {1e4, -1e3, 1e-10}};

  for (unsigned i = 0; i < lines.size(); ++i) {
    REQUIRE(seissol::writer::parseReceiverLine(lines[i]) == expectedPoints[i]);
  }
}

TEST_CASE("Throws expected exceptions for incorrect number of columns") {
  CHECK_THROWS_AS(seissol::writer::parseReceiverLine("0.0 0.1"), std::runtime_error);
  CHECK_THROWS_AS(seissol::writer::parseReceiverLine("0.0 0.1 0.2 0.3"), std::runtime_error);
}

TEST_CASE("Throws expected exceptions for conversion errors") {
  CHECK_THROWS_AS(seissol::writer::parseReceiverLine("0.1 center 0.3"), std::invalid_argument);
}

TEST_CASE("Parses receiver file correctly") {
  const std::string receiverFileName = "Testing/receiver_correct.dat";
  const auto points = seissol::writer::parseReceiverFile(receiverFileName);

  const auto expectedPoints = std::vector<Eigen::Vector3d>{{1, 0.1, 10}, {10, 2, 0.2}};

  REQUIRE(points.size() == expectedPoints.size());
  for (auto i = 0U; i < points.size(); ++i) {
    REQUIRE(points[i] == expectedPoints[i]);
  }
}
} // namespace seissol::unit_test
