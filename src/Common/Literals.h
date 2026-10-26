// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_LITERALS_H_
#define SEISSOL_SRC_COMMON_LITERALS_H_

#include <cstddef>
#include <cstdint>

namespace seissol {

// custom literals; use `_` first

constexpr std::int8_t operator""_i8(unsigned long long value) {
  return static_cast<std::int8_t>(value);
}

constexpr std::int16_t operator""_i16(unsigned long long value) {
  return static_cast<std::int16_t>(value);
}

constexpr std::int32_t operator""_i32(unsigned long long value) {
  return static_cast<std::int32_t>(value);
}

constexpr std::int64_t operator""_i64(unsigned long long value) {
  return static_cast<std::int64_t>(value);
}

constexpr std::uint8_t operator""_u8(unsigned long long value) {
  return static_cast<std::uint8_t>(value);
}

constexpr std::uint16_t operator""_u16(unsigned long long value) {
  return static_cast<std::uint16_t>(value);
}

constexpr std::uint32_t operator""_u32(unsigned long long value) {
  return static_cast<std::uint32_t>(value);
}

constexpr std::uint64_t operator""_u64(unsigned long long value) {
  return static_cast<std::uint64_t>(value);
}

constexpr std::int8_t operator""_I8(unsigned long long value) {
  return static_cast<std::int8_t>(value);
}

constexpr std::int16_t operator""_I16(unsigned long long value) {
  return static_cast<std::int16_t>(value);
}

constexpr std::int32_t operator""_I32(unsigned long long value) {
  return static_cast<std::int32_t>(value);
}

constexpr std::int64_t operator""_I64(unsigned long long value) {
  return static_cast<std::int64_t>(value);
}

constexpr std::uint8_t operator""_U8(unsigned long long value) {
  return static_cast<std::uint8_t>(value);
}

constexpr std::uint16_t operator""_U16(unsigned long long value) {
  return static_cast<std::uint16_t>(value);
}

constexpr std::uint32_t operator""_U32(unsigned long long value) {
  return static_cast<std::uint32_t>(value);
}

constexpr std::uint64_t operator""_U64(unsigned long long value) {
  return static_cast<std::uint64_t>(value);
}

// backport from C++23 (almost)

constexpr ptrdiff_t operator""_z(unsigned long long value) { return static_cast<ptrdiff_t>(value); }

constexpr ptrdiff_t operator""_Z(unsigned long long value) { return static_cast<ptrdiff_t>(value); }

constexpr std::size_t operator""_uz(unsigned long long value) {
  return static_cast<std::size_t>(value);
}

constexpr std::size_t operator""_UZ(unsigned long long value) {
  return static_cast<std::size_t>(value);
}

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_LITERALS_H_
