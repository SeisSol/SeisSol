// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LTSSTRUCTUREPARSER_20231019_HPP
#define LTSSTRUCTUREPARSER_20231019_HPP

#include "JsonToken.hpp"
#include "JsonLexer.hpp"
#include "LtsStructure.hpp"

#include <iosfwd>
#include <vector>
#include <optional>
#include <unordered_map>
#include <variant>

namespace seissol {

class LtsStructureParser {
  public:
  LtsStructureParser(std::string const& input, std::ostream* os = nullptr);
  auto parse() -> std::optional<std::vector<TimeCluster>>;
  auto parseCluster() -> std::optional<TimeCluster>;
  auto parseRegion() -> std::optional<Region>;

  private:
  auto getNextToken(char const* file = nullptr, int line = -1) -> JsonTokenKind;
  bool expect(JsonTokenKind kind, char const* file = nullptr, int line = -1);

  template <typename ObjT>
  auto parseArrayOfObject(std::optional<ObjT> (LtsStructureParser::*parseFun)())
      -> std::optional<std::vector<ObjT>> {
    if (!expect(JsonTokenKind::BeginArray)) {
      return std::nullopt;
    }

    auto objects = std::vector<ObjT>{};
    getNextToken(__FILE__, __LINE__);
    if (kind_ != JsonTokenKind::EndArray) {
      for (;;) {
        auto object = (this->*parseFun)();
        if (!object) {
          return std::nullopt;
        }
        objects.push_back(std::move(*object));

        if (kind_ != JsonTokenKind::ValueSeparator) {
          break;
        }
        getNextToken(__FILE__, __LINE__);
      }
    }

    if (!expect(JsonTokenKind::EndArray)) {
      return std::nullopt;
    }
    getNextToken(__FILE__, __LINE__);

    return std::make_optional(std::move(objects));
  }

  template <typename T, typename C>
  bool set(std::string const& key, std::unordered_map<std::string, T C::*> const& map, C& c) {
    if (auto it = map.find(key); it != map.end()) {
      if (!expect(to_kind_v<T>)) {
        return false;
      }
      c.*(it->second) = std::get<T>(tok_.val);
      return true;
    }
    return false;
  }

  static const std::unordered_map<std::string, long TimeCluster::*> TimeClusterIntFields;
  static const std::unordered_map<std::string, double TimeCluster::*> TimeClusterDoubleFields;
  static const std::unordered_map<std::string, long Region::*> RegionIntFields;
  static const std::unordered_map<std::string, double Region::*> RegionDoubleFields;

  JsonLexer lex_;
  std::ostream* os_;

  JsonTokenKind kind_;
  JsonToken tok_;
};

} // namespace seissol

#endif // LTSSTRUCTUREPARSER_20231019_HPP

