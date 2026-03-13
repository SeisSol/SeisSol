// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/MPIType.h"

#include <memory>
#include <mpi.h>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::io::datatype;

TEST_CASE("convertToMPI F32" * doctest::test_suite("io")) {
  auto dt = std::make_shared<F32Datatype>();
  MPI_Datatype mpit = convertToMPI(dt);
  CHECK(mpit == MPI_FLOAT);
}

TEST_CASE("convertToMPI F64" * doctest::test_suite("io")) {
  auto dt = std::make_shared<F64Datatype>();
  MPI_Datatype mpit = convertToMPI(dt);
  CHECK(mpit == MPI_DOUBLE);
}

TEST_CASE("convertToMPI F80" * doctest::test_suite("io")) {
  auto dt = std::make_shared<F80Datatype>();
  MPI_Datatype mpit = convertToMPI(dt);
  CHECK(mpit == MPI_LONG_DOUBLE);
}

TEST_CASE("convertToMPI IntegerDatatype" * doctest::test_suite("io")) {
  SUBCASE("int32 signed") {
    auto dt = std::make_shared<IntegerDatatype>(4, true);
    MPI_Datatype mpit = convertToMPI(dt);
    int size = 0;
    MPI_Type_size(mpit, &size);
    CHECK(size == 4);
  }

  SUBCASE("uint32 unsigned") {
    auto dt = std::make_shared<IntegerDatatype>(4, false);
    MPI_Datatype mpit = convertToMPI(dt);
    int size = 0;
    MPI_Type_size(mpit, &size);
    CHECK(size == 4);
  }

  SUBCASE("int8") {
    auto dt = std::make_shared<IntegerDatatype>(1, true);
    MPI_Datatype mpit = convertToMPI(dt);
    int size = 0;
    MPI_Type_size(mpit, &size);
    CHECK(size == 1);
  }

  SUBCASE("int64 unsigned") {
    auto dt = std::make_shared<IntegerDatatype>(8, false);
    MPI_Datatype mpit = convertToMPI(dt);
    int size = 0;
    MPI_Type_size(mpit, &size);
    CHECK(size == 8);
  }

  SUBCASE("int16 signed") {
    auto dt = std::make_shared<IntegerDatatype>(2, true);
    MPI_Datatype mpit = convertToMPI(dt);
    int size = 0;
    MPI_Type_size(mpit, &size);
    CHECK(size == 2);
  }
}

TEST_CASE("convertToMPI ArrayDatatype" * doctest::test_suite("io")) {
  auto base = std::make_shared<F64Datatype>();
  auto dt = std::make_shared<ArrayDatatype>(base, std::vector<std::size_t>{4, 3});
  MPI_Datatype mpit = convertToMPI(dt);
  int size = 0;
  MPI_Type_size(mpit, &size);
  CHECK(size == 4 * 3 * 8);
}

TEST_CASE("convertToMPI StructDatatype" * doctest::test_suite("io")) {
  auto f64 = std::make_shared<F64Datatype>();
  auto i32 = std::make_shared<IntegerDatatype>(4, true);
  std::vector<StructDatatype::MemberInfo> members = {
      {"a", 0, f64},
      {"b", 8, i32},
  };
  auto dt = std::make_shared<StructDatatype>(members, 16);
  MPI_Datatype mpit = convertToMPI(dt);

  // The struct type after create_resized should report the full padded size
  int size = 0;
  MPI_Type_size(mpit, &size);
  // MPI_Type_size returns the actual data size (without padding)
  CHECK(size == 8 + 4);

  MPI_Aint lb = 0;
  MPI_Aint extent = 0;
  MPI_Type_get_extent(mpit, &lb, &extent);
  CHECK(extent == 16);
}

TEST_CASE("convertToMPI OpaqueDatatype" * doctest::test_suite("io")) {
  auto dt = std::make_shared<OpaqueDatatype>(7);
  MPI_Datatype mpit = convertToMPI(dt);
  int size = 0;
  MPI_Type_size(mpit, &size);
  CHECK(size == 7);
}

TEST_CASE("convertToMPI StringDatatype" * doctest::test_suite("io")) {
  auto dt = std::make_shared<StringDatatype>(32);
  MPI_Datatype mpit = convertToMPI(dt);
  int size = 0;
  MPI_Type_size(mpit, &size);
  CHECK(size == 32);
}

TEST_CASE("convertToMPI caching" * doctest::test_suite("io")) {
  // Calling convertToMPI twice with the same datatype should return the same handle
  auto dt1 = std::make_shared<F32Datatype>();
  auto dt2 = std::make_shared<F32Datatype>();
  MPI_Datatype t1 = convertToMPI(dt1);
  MPI_Datatype t2 = convertToMPI(dt2);
  CHECK(t1 == t2);
}

} // namespace seissol::unit_test
