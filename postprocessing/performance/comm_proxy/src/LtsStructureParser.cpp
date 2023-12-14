// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "LtsStructureParser.hpp"

#include <ostream>
#include <utility>

namespace seissol {

const std::unordered_map<std::string, long TimeCluster::*>
    LtsStructureParser::TimeClusterIntFields = {
        {"time_cluster_id", &TimeCluster::timeClusterId},
        {"time_step_rate", &TimeCluster::timeStepRate},
        {"number_of_interior_cells", &TimeCluster::numberOfInteriorCells},
        {"interior_layer_size", &TimeCluster::interiorLayerSize}};
const std::unordered_map<std::string, double TimeCluster::*>
    LtsStructureParser::TimeClusterDoubleFields = {
        {"cfl_time_step_width", &TimeCluster::cflTimeStepWidth}};
const std::unordered_map<std::string, long Region::*> LtsStructureParser::RegionIntFields = {
    {"neighbor_rank", &Region::neighborRank},
    {"neighbor_time_cluster_id", &Region::neighborTimeClusterId},
    {"neighbor_time_step_rate", &Region::neighborTimeStepRate},
    {"number_of_ghost_region_cells", &Region::numberOfGhostRegionCells},
    {"number_of_ghost_region_derivatives", &Region::numberOfGhostRegionDerivatives},
    {"number_of_copy_region_cells", &Region::numberOfCopyRegionCells},
    {"number_of_communicated_copy_region_derivatives",
     &Region::numberOfCommunicatedCopyRegionDerivatives},
    {"send_identifier", &Region::sendIdentifier},
    {"receive_identifier", &Region::receiveIdentifier},
    {"copy_region_size", &Region::copyRegionSize},
    {"ghost_region_size", &Region::ghostRegionSize}};
const std::unordered_map<std::string, double Region::*> LtsStructureParser::RegionDoubleFields = {
    {"neighbor_cfl_time_step_width", &Region::neighborCflTimeStepWidth}};

LtsStructureParser::LtsStructureParser(std::string const& input, std::ostream* os)
    : lex_{input.c_str()}, os_(os) {}

auto LtsStructureParser::parse() -> std::optional<std::vector<TimeCluster>> {
  getNextToken(__FILE__, __LINE__);
  return parseArrayOfObject(&LtsStructureParser::parseCluster);
}

auto LtsStructureParser::parseCluster() -> std::optional<TimeCluster> {
  if (!expect(JsonTokenKind::BeginObject, __FILE__, __LINE__)) {
    return std::nullopt;
  }
  getNextToken(__FILE__, __LINE__);

  auto tc = TimeCluster{};

  for (;;) {
    if (!expect(JsonTokenKind::String, __FILE__, __LINE__)) {
      return std::nullopt;
    }

    auto key = std::get<std::string>(tok_.val);
    getNextToken(__FILE__, __LINE__);

    if (!expect(JsonTokenKind::NameSeparator, __FILE__, __LINE__)) {
      return std::nullopt;
    }
    getNextToken(__FILE__, __LINE__);

    if (key == "regions") {
      auto regions = parseArrayOfObject(&LtsStructureParser::parseRegion);
      if (regions) {
        tc.regions = std::move(*regions);
      } else {
        return std::nullopt;
      }
    } else {
      if (!set(key, TimeClusterIntFields, tc)) {
        if (!set(key, TimeClusterDoubleFields, tc)) {
          if (os_) {
            *os_ << "parseCluster: Key " << key << " is unknown or value type is wrong."
                 << std::endl;
          }
          return std::nullopt;
        }
      }
      getNextToken(__FILE__, __LINE__);
    }

    if (kind_ != JsonTokenKind::ValueSeparator) {
      break;
    }
    getNextToken(__FILE__, __LINE__);
  }

  if (!expect(JsonTokenKind::EndObject, __FILE__, __LINE__)) {
    return std::nullopt;
  }
  getNextToken(__FILE__, __LINE__);

  return std::make_optional(tc);
}

auto LtsStructureParser::parseRegion() -> std::optional<Region> {
  if (!expect(JsonTokenKind::BeginObject, __FILE__, __LINE__)) {
    return std::nullopt;
  }
  getNextToken(__FILE__, __LINE__);

  auto reg = Region{};

  for (;;) {
    if (!expect(JsonTokenKind::String, __FILE__, __LINE__)) {
      return std::nullopt;
    }

    auto key = std::get<std::string>(tok_.val);
    getNextToken(__FILE__, __LINE__);

    if (!expect(JsonTokenKind::NameSeparator, __FILE__, __LINE__)) {
      return std::nullopt;
    }
    getNextToken(__FILE__, __LINE__);

    if (!set(key, RegionIntFields, reg)) {
      if (!set(key, RegionDoubleFields, reg)) {
        if (os_) {
          *os_ << "parseRegion: Key " << key << " is unknown or value type is wrong." << std::endl;
        }
        return std::nullopt;
      }
    }
    getNextToken(__FILE__, __LINE__);

    if (kind_ != JsonTokenKind::ValueSeparator) {
      break;
    }
    getNextToken(__FILE__, __LINE__);
  }

  if (!expect(JsonTokenKind::EndObject, __FILE__, __LINE__)) {
    return std::nullopt;
  }
  getNextToken(__FILE__, __LINE__);

  return std::make_optional(reg);
}

auto LtsStructureParser::getNextToken(char const* file, int line) -> JsonTokenKind {
  kind_ = lex_(tok_);
  if (kind_ == JsonTokenKind::Unknown) {
    if (os_) {
      *os_ << "Lexer error: ";
      switch (std::get<JsonLexerError>(tok_.val)) {
      case JsonLexerError::None:
        *os_ << "Unknown token";
        break;
      case JsonLexerError::IntegerOverflow:
        *os_ << "Integer overflow";
        break;
      case JsonLexerError::FloatingOutOfRange:
        *os_ << "Floating out of range";
        break;
      case JsonLexerError::InvalidString:
        *os_ << "Invalid string";
        break;
      default:
        break;
      }
      if (file) {
        *os_ << " in " << file;
      }
      if (line >= 0) {
        *os_ << "@" << line;
      }
      *os_ << std::endl;
    }
  }
  return kind_;
}

bool LtsStructureParser::expect(JsonTokenKind kind, char const* file, int line) {
  if (kind_ != kind) {
    if (os_) {
      *os_ << "Expected " << to_string(kind) << " but got " << to_string(kind_);
      if (file) {
        *os_ << " in " << file;
      }
      if (line >= 0) {
        *os_ << "@" << line;
      }
      *os_ << std::endl;
    }
    return false;
  }
  return true;
}

} // namespace seissol
