// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Unit.h"

#include <cmath>
#include <ios>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace {
const std::vector<std::string> PositivePrefixes = {
    "k", "M", "G", "T", "P", "E", "Z", "Y", "R", "Q"};
const std::vector<std::string> PositiveBytePrefixes = {
    "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi"};
const std::vector<std::string> NegativePrefixes = {
    "m", "µ", "n", "p", "f", "a", "z", "y", "r", "q"};
} // namespace

namespace seissol {
SIUnit::SIUnit(const std::string& unit, bool binary) : unit(unit), binary(binary) {}

std::string SIUnit::formatTime(double value, bool exact, int digits) const {
  const double byDay = std::floor(value / (60 * 60 * 24));
  const double hours = value - byDay * (60 * 60 * 24);
  const double byHour = std::floor(hours / (60 * 60));
  const double minutes = hours - byHour * (60 * 60);
  const double byMinute = std::floor(minutes / 60);
  const double seconds = minutes - byMinute * (60);

  std::ostringstream stream;
  stream.precision(0);
  stream << std::fixed;
  bool done = false;
  if (byDay > 0 && (!done || exact)) {
    stream << byDay << " d";
    done = true;
  }
  if (byHour > 0 && (!done || exact)) {
    if (done) {
      stream << " ";
    }
    stream << byHour << " h";
    done = true;
  }
  if (byMinute > 0 && (!done || exact)) {
    if (done) {
      stream << " ";
    }
    stream << byMinute << " min";
    done = true;
  }
  if (seconds != 0 || !done || exact) {
    if (done) {
      stream << " ";
    }
    stream << formatPrefix(seconds, {}, digits);
    done = true;
  }
  return stream.str();
}

std::string SIUnit::formatPrefix(double value, std::optional<double> error, int digits) const {
  double mantissa = std::abs(value);
  double errorOrZero = error.value_or(0);
  int position = 0;
  const double skip = binary ? 1024 : 1000;
  const double sign = value < 0 ? -1 : 1;
  // only one of the following two while loops should be triggered at any time
  // the 100 is rather arbitrary to prevent a loop forever
  if (mantissa != 0) {
    while (mantissa < 1 && position > -100) {
      mantissa *= skip;
      errorOrZero *= skip;
      --position;
    }
  }
  while (mantissa >= skip && position < 100) {
    mantissa /= skip;
    errorOrZero /= skip;
    ++position;
  }

  if ((binary && position > static_cast<int>(PositiveBytePrefixes.size())) ||
      (!binary && position > static_cast<int>(PositivePrefixes.size())) ||
      -position > static_cast<int>(NegativePrefixes.size())) {
    // out of range, default to scientific notation
    return formatScientific(value, error, digits);
  } else {
    const std::string prefix = [&]() {
      if (position < 0) {
        return NegativePrefixes[-position - 1];
      } else if (position > 0) {
        if (binary) {
          return PositiveBytePrefixes[position - 1];
        } else {
          return PositivePrefixes[position - 1];
        }
      } else {
        return std::string("");
      }
    }();

    std::ostringstream stream;
    stream.precision(digits);
    stream << std::fixed;
    if (error.has_value()) {
      stream << "(";
    }
    stream << sign * mantissa;
    if (error.has_value()) {
      stream << " ± " << errorOrZero << ")";
    }
    stream << " " << prefix << unit;
    return stream.str();
  }
}

std::string SIUnit::formatScientific(double value, std::optional<double> error, int digits) const {
  std::ostringstream stream;
  stream.precision(digits);
  stream << std::scientific;
  if (error.has_value()) {
    stream << "(";
  }
  stream << value;
  if (error.has_value()) {
    stream << " ± " << error.value() << ")";
  }
  stream << " " << unit;
  return stream.str();
}

} // namespace seissol
