// SPDX-FileCopyrightText: 2024 SeisSol Group
// Copyright (C) 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause


#ifndef SEISSOL_SRC_COMMON_FNV1A_H_
#define SEISSOL_SRC_COMMON_FNV1A_H_

#include <cstdint>
#include <string_view>

namespace seissol {

constexpr std::uint64_t fnv1a(char const* str, std::size_t n) {
  constexpr std::uint64_t basis = 0xcbf29ce484222325;
  constexpr std::uint64_t prime = 0x00000100000001b3;
  return n > 0 ? (fnv1a(str, n - 1) ^ str[n - 1]) * prime : basis;
}
constexpr std::uint64_t fnv1a(std::string_view str) { return fnv1a(str.data(), str.size()); }

namespace literals {

constexpr std::uint64_t operator""_fnv1a(char const* str, std::size_t n) {
  return seissol::fnv1a(str, n);
}

} // namespace literals

} // namespace seissol


#endif // SEISSOL_SRC_COMMON_FNV1A_H_

