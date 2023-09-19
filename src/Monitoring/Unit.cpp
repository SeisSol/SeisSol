#include "Unit.hpp"

#include <string>
#include <vector>
#include <cmath>
#include <sstream>

namespace {
const static std::vector<std::string> PositivePrefixes = {
    "k", "M", "G", "T", "P", "E", "Z", "Y", "R", "Q"};
const static std::vector<std::string> PositiveBytePrefixes = {
    "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi"};
const static std::vector<std::string> NegativePrefixes = {
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
    stream << formatPrefix(seconds, digits);
    done = true;
  }
  return stream.str();
}

std::string SIUnit::formatPrefix(double value, int digits) const {
  double adjValue = std::abs(value);
  int position = 0;
  const double skip = binary ? 1024 : 1000;
  const double sign = value < 0 ? -1 : 1;
  // only one of the following two while loops should be triggered at any time
  if (adjValue != 0) {
    while (adjValue < 1) {
      adjValue *= skip;
      --position;
    }
  }
  while (adjValue >= skip) {
    adjValue /= skip;
    ++position;
  }

  if ((binary && position > static_cast<int>(PositiveBytePrefixes.size())) ||
      (!binary && position > static_cast<int>(PositivePrefixes.size())) ||
      -position > static_cast<int>(NegativePrefixes.size())) {
    // out of range, default to scientific notation
    return formatScientific(value, digits);
  }

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
  stream << sign * adjValue << " " << prefix << unit;
  return stream.str();
}

std::string SIUnit::formatScientific(double value, int digits) const {
  const int exponential = [&]() {
    if (value == 0) {
      return 0.0;
    } else {
      const double log = std::log10(std::abs(value));
      return std::floor(log);
    }
  }();
  const double rest = value / std::pow(10, exponential);
  std::ostringstream stream;
  stream.precision(digits);
  stream << std::fixed;
  stream << rest << " * 10^";
  if (exponential < 0) {
    stream << "(" << exponential << ")";
  } else {
    stream << exponential;
  }
  stream << " " << unit;
  return stream.str();
}

} // namespace seissol
