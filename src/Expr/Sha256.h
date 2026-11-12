// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_SHA256_H_
#define SEISSOL_SRC_EXPR_SHA256_H_

// SHA-256 (FIPS 180-4), for the provenance of the core format: the hashes of the sources a core
// file was generated from, and of its body (`meta source`, `depends`, `body`); and for the files of
// the grids, which the evaluation record names by their hashes.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

namespace seissol::expr {

/// The SHA-256 of data given in pieces.
class Sha256 {
  public:
  void update(const void* data, std::size_t size);
  /// The hash of all data given, as 64 lowercase hexadecimal digits. Spends the object.
  [[nodiscard]] std::string finish();

  private:
  std::array<std::uint32_t, 8> state_{0x6a09e667U,
                                      0xbb67ae85U,
                                      0x3c6ef372U,
                                      0xa54ff53aU,
                                      0x510e527fU,
                                      0x9b05688cU,
                                      0x1f83d9abU,
                                      0x5be0cd19U};
  std::array<unsigned char, 64> block_{};
  std::size_t filled_{0};
  std::uint64_t length_{0};
};

/// The SHA-256 of `data`, as 64 lowercase hexadecimal digits.
[[nodiscard]] std::string sha256(const std::string& data);

/// The SHA-256 of the file `path`: from the cache of what SeisSol derives across runs
/// (cacheDirectory()) where it holds one for the file as it is -- its path, its size and the time
/// it changed --, else read, and kept there. Empty where the file cannot be read.
[[nodiscard]] std::string sha256OfFile(const std::string& path);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_SHA256_H_
