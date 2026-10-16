// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Writer/File/BinaryWriter.h"
#include "WriterHarness.t.h"

#include <cstddef>
#include <fstream>
#include <iterator>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace binaryfiletest {

inline std::vector<char> pattern(std::size_t size, char first) {
  std::vector<char> bytes(size);
  for (std::size_t i = 0; i < size; ++i) {
    bytes[i] = static_cast<char>(first + static_cast<char>(i % 23));
  }
  return bytes;
}

inline std::vector<char> contents(const std::string& file) {
  std::ifstream stream(file, std::ios::binary);
  return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

} // namespace binaryfiletest

/**
 * MPI counts bytes in an int, so a write of more than 2 GiB has to go in pieces. With pieces of a
 * few bytes, the same logic shows on small data: the file has to hold every byte, in order.
 */
TEST_CASE("IO/BinaryFile: a write larger than one piece goes in pieces" *
          doctest::test_suite("io")) {
  using namespace binaryfiletest;
  const unit_test::io::TempDir dir;
  const auto file = dir.path + "/pieces.bin";

  const auto head = pattern(100, 'a');
  const auto body = pattern(53, 'A');
  const auto patch = pattern(30, '0');
  {
    seissol::io::writer::file::BinaryFile binary(MPI_COMM_SELF, false, 7);
    binary.openFile(file, false);
    binary.writeGlobal(head.data(), head.size());
    binary.writeDistributed(body.data(), body.size());
    binary.writeGlobalAt(patch.data(), patch.size(), 10);
    binary.closeFile();
  }

  auto expected = head;
  expected.insert(expected.end(), body.begin(), body.end());
  std::copy(patch.begin(), patch.end(), expected.begin() + 10);
  CHECK(contents(file) == expected);
}

} // namespace seissol::unit_test
