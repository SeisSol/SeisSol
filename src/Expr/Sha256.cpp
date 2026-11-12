// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Sha256.h"

#include "Expr/RtcCpu.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <system_error>
#include <vector>

namespace seissol::expr {

namespace {

constexpr std::array<std::uint32_t, 64> RoundConstants = {
    0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U, 0x923f82a4U,
    0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU,
    0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU,
    0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU, 0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
    0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU,
    0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
    0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U, 0x19a4c116U,
    0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
    0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U, 0x90befffaU, 0xa4506cebU, 0xbef9a3f7U,
    0xc67178f2U};

constexpr std::uint32_t rotate(std::uint32_t x, unsigned n) { return (x >> n) | (x << (32U - n)); }

void compress(std::array<std::uint32_t, 8>& state, const unsigned char* block) {
  std::array<std::uint32_t, 64> w{};
  for (std::size_t i = 0; i < 16; ++i) {
    w[i] = (std::uint32_t{block[4 * i]} << 24U) | (std::uint32_t{block[4 * i + 1]} << 16U) |
           (std::uint32_t{block[4 * i + 2]} << 8U) | std::uint32_t{block[4 * i + 3]};
  }
  for (std::size_t i = 16; i < 64; ++i) {
    const std::uint32_t s0 = rotate(w[i - 15], 7) ^ rotate(w[i - 15], 18) ^ (w[i - 15] >> 3U);
    const std::uint32_t s1 = rotate(w[i - 2], 17) ^ rotate(w[i - 2], 19) ^ (w[i - 2] >> 10U);
    w[i] = w[i - 16] + s0 + w[i - 7] + s1;
  }
  std::array<std::uint32_t, 8> v = state;
  for (std::size_t i = 0; i < 64; ++i) {
    const std::uint32_t s1 = rotate(v[4], 6) ^ rotate(v[4], 11) ^ rotate(v[4], 25);
    const std::uint32_t choice = (v[4] & v[5]) ^ (~v[4] & v[6]);
    const std::uint32_t t1 = v[7] + s1 + choice + RoundConstants[i] + w[i];
    const std::uint32_t s0 = rotate(v[0], 2) ^ rotate(v[0], 13) ^ rotate(v[0], 22);
    const std::uint32_t majority = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
    const std::uint32_t t2 = s0 + majority;
    v[7] = v[6];
    v[6] = v[5];
    v[5] = v[4];
    v[4] = v[3] + t1;
    v[3] = v[2];
    v[2] = v[1];
    v[1] = v[0];
    v[0] = t1 + t2;
  }
  for (std::size_t i = 0; i < 8; ++i) {
    state[i] += v[i];
  }
}

} // namespace

void Sha256::update(const void* data, std::size_t size) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  length_ += size;
  while (size > 0) {
    if (filled_ == 0 && size >= 64) {
      compress(state_, bytes);
      bytes += 64;
      size -= 64;
      continue;
    }
    const std::size_t take = std::min<std::size_t>(64 - filled_, size);
    std::copy_n(bytes, take, block_.data() + filled_);
    filled_ += take;
    bytes += take;
    size -= take;
    if (filled_ == 64) {
      compress(state_, block_.data());
      filled_ = 0;
    }
  }
}

std::string Sha256::finish() {
  // the rest, a one bit, zeros, and the length in bits
  const std::uint64_t bits = length_ * 8U;
  std::array<unsigned char, 128> tail{};
  std::copy_n(block_.data(), filled_, tail.data());
  tail[filled_] = 0x80U;
  const std::size_t blocks = filled_ + 1 + 8 > 64 ? 2 : 1;
  for (std::size_t i = 0; i < 8; ++i) {
    tail[blocks * 64 - 1 - i] = static_cast<unsigned char>(bits >> (8U * i));
  }
  for (std::size_t b = 0; b < blocks; ++b) {
    compress(state_, tail.data() + 64 * b);
  }
  static constexpr char Digits[] = "0123456789abcdef";
  std::string hex;
  hex.reserve(64);
  for (const std::uint32_t word : state_) {
    for (int shift = 28; shift >= 0; shift -= 4) {
      hex.push_back(Digits[(word >> static_cast<unsigned>(shift)) & 0xFU]);
    }
  }
  return hex;
}

std::string sha256(const std::string& data) {
  Sha256 hash;
  hash.update(data.data(), data.size());
  return hash.finish();
}

std::string sha256OfFile(const std::string& path) {
  std::error_code error;
  const auto absolute = std::filesystem::absolute(path, error).lexically_normal().string();
  const auto size = std::filesystem::file_size(path, error);
  if (error) {
    return {};
  }
  const auto changed = std::filesystem::last_write_time(path, error).time_since_epoch().count();
  if (error) {
    return {};
  }
  // the entry of the cache: the size, the time of change and the hash, then the path
  const std::string identity = std::to_string(size) + " " + std::to_string(changed);
  const std::string& directory = cacheDirectory();
  const std::filesystem::path entry =
      directory.empty() ? std::filesystem::path()
                        : std::filesystem::path(directory) / (sha256(absolute) + ".file-hash");
  if (!directory.empty()) {
    std::ifstream cached(entry);
    std::string line;
    std::string cachedPath;
    if (std::getline(cached, line) && std::getline(cached, cachedPath) && cachedPath == absolute &&
        line.size() > 65 && line.compare(0, identity.size() + 1, identity + " ") == 0) {
      return line.substr(identity.size() + 1);
    }
  }

  std::ifstream file(path, std::ios::binary);
  if (!file) {
    return {};
  }
  Sha256 hash;
  std::vector<char> buffer(std::size_t{4} << 20);
  while (file) {
    file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    hash.update(buffer.data(), static_cast<std::size_t>(file.gcount()));
  }
  if (file.bad()) {
    return {};
  }
  const std::string hex = hash.finish();

  if (!directory.empty()) {
    // to a file of its own, renamed into place, so that a reader finds a whole entry
    const std::filesystem::path temporary =
        entry.string() + ".tmp" + std::to_string(std::random_device{}());
    {
      std::ofstream out(temporary);
      out << identity << " " << hex << "\n" << absolute << "\n";
    }
    std::filesystem::rename(temporary, entry, error);
    if (error) {
      std::filesystem::remove(temporary, error);
    }
  }
  return hex;
}

} // namespace seissol::expr
