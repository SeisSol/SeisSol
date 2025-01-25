// SPDX-FileCopyrightText: 2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_FNV1A_H_
#define SEISSOL_SRC_COMMON_FNV1A_H_

#include <cstdint>
#include <string_view>

namespace seissol {

constexpr auto fnv1a(const char* str, std::size_t n) -> std::uint64_t {
  constexpr std::uint64_t Basis = 0xcbf29ce484222325ULL;
  constexpr std::uint64_t Prime = 0x00000100000001b3ULL;
  return n > 0 ? (fnv1a(str, n - 1) ^ str[n - 1]) * Prime : Basis;
}
constexpr auto fnv1a(const std::string_view& str) -> std::uint64_t {
  return fnv1a(str.data(), str.size());
}

namespace literals {

constexpr auto operator""_fnv1a(const char* str, std::size_t n) -> std::uint64_t {
  return seissol::fnv1a(str, n);
}

} // namespace literals

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_FNV1A_H_
