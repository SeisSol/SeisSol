#pragma once

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace seissol {
struct SIUnit {
  public:
  SIUnit(const std::string& unit, bool binary);

  std::string formatTime(double value, bool exact = true, int digits = 4) const;

  std::string formatPrefix(double value, int digits = 4) const;

  std::string formatScientific(double value, int digits = 4) const;

  private:
  std::string unit;
  bool binary;
};

const inline SIUnit UnitTime = SIUnit("s", false);
const inline SIUnit UnitFlop = SIUnit("FLOP", false);
const inline SIUnit UnitFlopPerS = SIUnit("FLOP/s", false);
const inline SIUnit UnitByte = SIUnit("B", true);
} // namespace seissol
