// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Json.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <locale>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <utils/logger.h>
#include <variant>

namespace seissol::json {

namespace {
void newline(std::ostream& stream, int indent, int depth) {
  if (indent >= 0) {
    stream << '\n' << std::string(static_cast<std::size_t>(indent) * depth, ' ');
  }
}
} // namespace

std::string quote(std::string_view text) {
  std::string result;
  result.reserve(text.size() + 2);
  result += '"';
  for (const char character : text) {
    switch (character) {
    case '"':
      result += "\\\"";
      break;
    case '\\':
      result += "\\\\";
      break;
    case '\b':
      result += "\\b";
      break;
    case '\f':
      result += "\\f";
      break;
    case '\n':
      result += "\\n";
      break;
    case '\r':
      result += "\\r";
      break;
    case '\t':
      result += "\\t";
      break;
    default:
      if (static_cast<unsigned char>(character) < 0x20) {
        // the remaining control characters have no short form
        constexpr std::string_view HexDigits = "0123456789abcdef";
        const auto code = static_cast<unsigned char>(character);
        result += "\\u00";
        result += HexDigits[code / 16];
        result += HexDigits[code % 16];
      } else {
        // everything else, UTF-8 included, is written as it is
        result += character;
      }
    }
  }
  result += '"';
  return result;
}

std::string formatNumber(double value) {
  if (!std::isfinite(value)) {
    return "null";
  }
  const auto format = [&](int digits) {
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << std::setprecision(digits) << value;
    return stream.str();
  };
  // 15 significant digits are exact for most values a person would write down; the 17 of
  // max_digits10 always are, but would show 0.1 as 0.10000000000000001
  for (const int digits : {15, 16}) {
    auto text = format(digits);
    if (std::strtod(text.c_str(), nullptr) == value) {
      return text;
    }
  }
  return format(17);
}

bool Value::isNumber() const {
  return std::holds_alternative<std::int64_t>(data_) ||
         std::holds_alternative<std::uint64_t>(data_) || std::holds_alternative<double>(data_);
}

Value& Value::operator[](std::string_view key) {
  if (isNull()) {
    data_ = Object{};
  }
  if (!isObject()) {
    logError() << "A JSON value that is not an object has no member" << std::string(key);
  }
  auto& object = std::get<Object>(data_);
  for (auto& [name, value] : object) {
    if (name == key) {
      return value;
    }
  }
  return object.emplace_back(std::string(key), Value()).second;
}

const Value* Value::find(std::string_view key) const {
  if (!isObject()) {
    return nullptr;
  }
  for (const auto& [name, value] : std::get<Object>(data_)) {
    if (name == key) {
      return &value;
    }
  }
  return nullptr;
}

Value& Value::push(Value value) {
  if (isNull()) {
    data_ = Array{};
  }
  if (!isArray()) {
    logError() << "Only a JSON array can be appended to.";
  }
  return std::get<Array>(data_).emplace_back(std::move(value));
}

Value& Value::at(std::size_t index) {
  if (!isArray() || index >= std::get<Array>(data_).size()) {
    logError() << "The JSON value has no element" << index;
  }
  return std::get<Array>(data_)[index];
}

std::size_t Value::size() const {
  if (isArray()) {
    return std::get<Array>(data_).size();
  }
  if (isObject()) {
    return std::get<Object>(data_).size();
  }
  return 0;
}

const std::string& Value::asString() const {
  if (!isString()) {
    logError() << "The JSON value is not a string.";
  }
  return std::get<std::string>(data_);
}

double Value::asDouble() const {
  if (std::holds_alternative<std::int64_t>(data_)) {
    return static_cast<double>(std::get<std::int64_t>(data_));
  }
  if (std::holds_alternative<std::uint64_t>(data_)) {
    return static_cast<double>(std::get<std::uint64_t>(data_));
  }
  if (!std::holds_alternative<double>(data_)) {
    logError() << "The JSON value is not a number.";
  }
  return std::get<double>(data_);
}

bool Value::asBool() const {
  if (!isBool()) {
    logError() << "The JSON value is not a boolean.";
  }
  return std::get<bool>(data_);
}

const Value::Array& Value::asArray() const {
  if (!isArray()) {
    logError() << "The JSON value is not an array.";
  }
  return std::get<Array>(data_);
}

const Value::Object& Value::asObject() const {
  if (!isObject()) {
    logError() << "The JSON value is not an object.";
  }
  return std::get<Object>(data_);
}

void Value::write(std::ostream& stream, int indent) const { write(stream, indent, 0); }

std::string Value::dump(int indent) const {
  std::ostringstream stream;
  stream.imbue(std::locale::classic());
  write(stream, indent);
  return stream.str();
}

void Value::write(std::ostream& stream, int indent, int depth) const {
  if (isNull()) {
    stream << "null";
  } else if (isBool()) {
    stream << (std::get<bool>(data_) ? "true" : "false");
  } else if (std::holds_alternative<std::int64_t>(data_)) {
    stream << std::get<std::int64_t>(data_);
  } else if (std::holds_alternative<std::uint64_t>(data_)) {
    stream << std::get<std::uint64_t>(data_);
  } else if (std::holds_alternative<double>(data_)) {
    stream << formatNumber(std::get<double>(data_));
  } else if (isString()) {
    stream << quote(std::get<std::string>(data_));
  } else if (isArray()) {
    const auto& array = std::get<Array>(data_);
    if (array.empty()) {
      stream << "[]";
      return;
    }
    stream << '[';
    for (std::size_t i = 0; i < array.size(); ++i) {
      if (i > 0) {
        stream << ',';
      }
      newline(stream, indent, depth + 1);
      array[i].write(stream, indent, depth + 1);
    }
    newline(stream, indent, depth);
    stream << ']';
  } else {
    const auto& object = std::get<Object>(data_);
    if (object.empty()) {
      stream << "{}";
      return;
    }
    stream << '{';
    for (std::size_t i = 0; i < object.size(); ++i) {
      if (i > 0) {
        stream << ',';
      }
      newline(stream, indent, depth + 1);
      stream << quote(object[i].first) << (indent >= 0 ? ": " : ":");
      object[i].second.write(stream, indent, depth + 1);
    }
    newline(stream, indent, depth);
    stream << '}';
  }
}

} // namespace seissol::json
