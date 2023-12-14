// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "JsonToken.hpp"

namespace seissol {
std::string to_string(JsonTokenKind kind) {
  switch (kind) {
  case JsonTokenKind::BeginArray:
    return "BeginArray";
  case JsonTokenKind::EndArray:
    return "EndArray";
  case JsonTokenKind::BeginObject:
    return "BeginObject";
  case JsonTokenKind::EndObject:
    return "EndObject";
  case JsonTokenKind::NameSeparator:
    return "NameSeparator";
  case JsonTokenKind::ValueSeparator:
    return "ValueSeparator";
  case JsonTokenKind::False:
    return "False";
  case JsonTokenKind::True:
    return "True";
  case JsonTokenKind::Null:
    return "Null";
  case JsonTokenKind::Object:
    return "Object";
  case JsonTokenKind::Array:
    return "Array";
  case JsonTokenKind::IntegerNumber:
    return "IntegerNumber";
  case JsonTokenKind::FloatingNumber:
    return "FloatingNumber";
  case JsonTokenKind::String:
    return "String";
  case JsonTokenKind::EOI:
    return "EOI";
  case JsonTokenKind::Unknown:
    return "Unknown";
  }
  return {};
}
} // namespace seissol
