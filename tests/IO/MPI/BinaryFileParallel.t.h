// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Writer/File/BinaryWriter.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <mpi.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace seissol::unit_test {

namespace binaryfileparalleltest {

//! What rank @p rank writes: a size of its own, so that the ranks need different numbers of pieces.
inline std::vector<char> dataOf(int rank) {
  std::vector<char> bytes(static_cast<std::size_t>(11 * (rank + 1)));
  for (std::size_t i = 0; i < bytes.size(); ++i) {
    bytes[i] = static_cast<char>('A' + rank + static_cast<int>(i % 7));
  }
  return bytes;
}

} // namespace binaryfileparalleltest

/**
 * A distributed write is collective, so a rank whose data takes fewer pieces than another's still
 * has to take part in every round. The file has to hold the data of all ranks, in rank order.
 */
TEST_CASE("IO/BinaryFile: ranks with different numbers of pieces write in rounds" *
          doctest::test_suite("io")) {
  using namespace binaryfileparalleltest;

  std::vector<char> name(256, '\0');
  if (Mpi::mpi.rank() == 0) {
    std::string pattern = "/tmp/seissolbinaryXXXXXX";
    const char* created = mkdtemp(pattern.data());
    REQUIRE(created != nullptr);
    std::copy(pattern.begin(), pattern.end(), name.begin());
  }
  MPI_Bcast(name.data(), static_cast<int>(name.size()), MPI_CHAR, 0, Mpi::mpi.comm());
  const std::string dir(name.data());
  const auto file = dir + "/rounds.bin";

  {
    const auto data = dataOf(Mpi::mpi.rank());
    seissol::io::writer::file::BinaryFile binary(Mpi::mpi.comm(), false, 5);
    binary.openFile(file, false);
    binary.writeDistributed(data.data(), data.size());
    binary.closeFile();
  }

  MPI_Barrier(Mpi::mpi.comm());
  if (Mpi::mpi.rank() == 0) {
    std::vector<char> expected;
    for (int rank = 0; rank < Mpi::mpi.size(); ++rank) {
      const auto data = dataOf(rank);
      expected.insert(expected.end(), data.begin(), data.end());
    }
    std::ifstream stream(file, std::ios::binary);
    const std::vector<char> written{std::istreambuf_iterator<char>(stream),
                                    std::istreambuf_iterator<char>()};
    CHECK(written == expected);
    std::error_code error;
    std::filesystem::remove_all(dir, error);
  }
}

} // namespace seissol::unit_test
