// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "RunId.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <string_view>

namespace seissol {

namespace {
constexpr std::array<std::size_t, 4> DashPositions{8, 13, 18, 23};

bool isDashPosition(std::size_t position) {
  return std::ranges::find(DashPositions, position) != DashPositions.end();
}
} // namespace

std::string newRunId() {
  std::random_device device;
  std::array<std::uint8_t, 16> bytes{};
  for (std::size_t i = 0; i < bytes.size(); i += 4) {
    const auto word = static_cast<std::uint32_t>(device());
    for (std::size_t j = 0; j < 4; ++j) {
      bytes[i + j] = static_cast<std::uint8_t>(word >> (8 * j));
    }
  }
  // RFC 9562: version 4 (random), variant 10xx
  bytes[6] = static_cast<std::uint8_t>((bytes[6] & 0x0fU) | 0x40U);
  bytes[8] = static_cast<std::uint8_t>((bytes[8] & 0x3fU) | 0x80U);

  constexpr std::string_view Digits = "0123456789abcdef";
  std::string text;
  text.reserve(RunIdLength);
  for (std::size_t i = 0; i < bytes.size(); ++i) {
    if (isDashPosition(text.size())) {
      text += '-';
    }
    text += Digits[bytes[i] >> 4U];
    text += Digits[bytes[i] & 0x0fU];
  }
  return text;
}

std::optional<std::string> parseRunId(std::string_view text) {
  if (text.size() != RunIdLength) {
    return {};
  }
  std::string normalized(text);
  for (std::size_t i = 0; i < normalized.size(); ++i) {
    const auto character = static_cast<unsigned char>(normalized[i]);
    if (isDashPosition(i)) {
      if (character != '-') {
        return {};
      }
    } else if (std::isxdigit(character) == 0) {
      return {};
    } else {
      normalized[i] = static_cast<char>(std::tolower(character));
    }
  }
  return normalized;
}

} // namespace seissol
