// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_UNIT_H_
#define SEISSOL_SRC_MONITORING_UNIT_H_

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace seissol {
struct SIUnit {
  public:
  SIUnit(const std::string& unit, bool binary);

  [[nodiscard]] std::string formatTime(double value, bool exact = true, int digits = 4) const;

  [[nodiscard]] std::string formatPrefix(double value, int digits = 4) const;

  [[nodiscard]] static std::string formatScientific(double value, int digits = 4);

  private:
  std::string unit;
  bool binary;
};

const inline SIUnit UnitTime = SIUnit("s", false);
const inline SIUnit UnitFlop = SIUnit("FLOP", false);
const inline SIUnit UnitFlopPerS = SIUnit("FLOP/s", false);
const inline SIUnit UnitByte = SIUnit("B", true);
} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_UNIT_H_
