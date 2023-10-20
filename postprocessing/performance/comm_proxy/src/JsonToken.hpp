// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef JSONTOKEN_20231019_HPP
#define JSONTOKEN_20231019_HPP

#include <limits>
#include <string>
#include <variant>

namespace seissol {

enum class JsonTokenKind : unsigned {
  // structural
  BeginArray,
  EndArray,
  BeginObject,
  EndObject,
  NameSeparator,
  ValueSeparator,
  // value
  False,
  True,
  Null,
  Object,
  Array,
  IntegerNumber,
  FloatingNumber,
  String,
  // end of input
  EOI,
  Unknown = std::numeric_limits<unsigned>::max()
};

std::string to_string(JsonTokenKind kind);

enum class JsonLexerError {
  None,
  UnknownToken,
  IntegerOverflow,
  FloatingOutOfRange,
  InvalidString
};

struct JsonToken {
  using yysval = std::variant<long, double, std::string, JsonLexerError>;

  yysval val;
};

template <typename T>
struct to_kind;

template <>
struct to_kind<long> {
  constexpr static auto value = JsonTokenKind::IntegerNumber;
};
template <>
struct to_kind<double> {
  constexpr static auto value = JsonTokenKind::FloatingNumber;
};
template <>
struct to_kind<std::string> {
  constexpr static auto value = JsonTokenKind::String;
};

template <typename T>
inline constexpr auto to_kind_v = to_kind<T>::value;

}; // namespace seissol

#endif // JSONTOKEN_20231019_HPP
