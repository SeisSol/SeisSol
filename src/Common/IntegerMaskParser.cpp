// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Common/IntegerMaskParser.h"
#include <regex>

namespace seissol {
std::optional<IntegerMaskParser::MaskType> IntegerMaskParser::parse(std::string mask) {
  MaskType resultMask{};

  auto const regex = std::regex("^(\\d+|\\d+-\\d+|\\{(\\d+,?)+\\})(,|$)");
  std::smatch match;

  while (std::regex_search(mask, match, regex)) {
    std::string subString = match[1].str();
    mask = match.suffix();

    auto intRange = IntegerMaskParser::parseIntRange(subString);
    if (intRange) {
      resultMask.push_back(*intRange);
      continue;
    }

    auto intList = IntegerMaskParser::parseIntList(subString);
    if (intList) {
      resultMask.push_back(*intList);
      continue;
    }

    resultMask.push_back({std::stoi(subString)});
  }

  if (resultMask.empty()) {
    return std::nullopt;
  }
  return resultMask;
}

IntegerMaskParser::OptionalIntVectorType IntegerMaskParser::parseIntRange(const std::string& str) {
  OptionalIntVectorType result{std::nullopt};
  std::size_t separator = str.find('-');
  if (separator != std::string::npos) {
    auto start = std::stoi(std::string(str, 0, separator));
    auto end = std::stoi(std::string(str, ++separator, str.size()));

    if (start > end) {
      std::swap(start, end);
    }
    std::vector<int> integers;
    for (int i{start}; i <= end; ++i) {
      integers.push_back(i);
    }
    result = integers;
  }
  return result;
}

IntegerMaskParser::OptionalIntVectorType IntegerMaskParser::parseIntList(const std::string& str) {
  OptionalIntVectorType result{std::nullopt};
  if (str[0] == '{') {
    std::string subStr{};
    std::vector<int> integers;

    const auto strCopy = std::string(str, 1, str.size());
    for (size_t i = 0; i < strCopy.size(); ++i) {
      if ((strCopy[i] == ',') || (strCopy[i] == '}')) {
        integers.push_back(std::stoi(subStr));
        subStr = std::string();
      } else {
        subStr.push_back(strCopy[i]);
      }
    }
    result = integers;
  }
  return result;
}
} // namespace seissol

