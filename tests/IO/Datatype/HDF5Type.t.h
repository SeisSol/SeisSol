// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/HDF5Type.h"

#include <hdf5.h>
#include <memory>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::io::datatype;

TEST_CASE("convertToHdf5 F32" * doctest::test_suite("io")) {
  auto dt = std::make_shared<F32Datatype>();
  hid_t h5t = convertToHdf5(dt);
  CHECK(H5Tget_class(h5t) == H5T_FLOAT);
  CHECK(H5Tget_size(h5t) == 4);
}

TEST_CASE("convertToHdf5 F64" * doctest::test_suite("io")) {
  auto dt = std::make_shared<F64Datatype>();
  hid_t h5t = convertToHdf5(dt);
  CHECK(H5Tget_class(h5t) == H5T_FLOAT);
  CHECK(H5Tget_size(h5t) == 8);
}

TEST_CASE("convertToHdf5 Integer" * doctest::test_suite("io")) {
  SUBCASE("int32 signed") {
    auto dt = std::make_shared<IntegerDatatype>(4, true);
    hid_t h5t = convertToHdf5(dt);
    CHECK(H5Tget_class(h5t) == H5T_INTEGER);
    CHECK(H5Tget_size(h5t) == 4);
    CHECK(H5Tget_sign(h5t) == H5T_SGN_2);
  }

  SUBCASE("uint32 unsigned") {
    auto dt = std::make_shared<IntegerDatatype>(4, false);
    hid_t h5t = convertToHdf5(dt);
    CHECK(H5Tget_class(h5t) == H5T_INTEGER);
    CHECK(H5Tget_size(h5t) == 4);
    CHECK(H5Tget_sign(h5t) == H5T_SGN_NONE);
  }

  SUBCASE("int8 signed") {
    auto dt = std::make_shared<IntegerDatatype>(1, true);
    hid_t h5t = convertToHdf5(dt);
    CHECK(H5Tget_class(h5t) == H5T_INTEGER);
    CHECK(H5Tget_size(h5t) == 1);
  }

  SUBCASE("int64 unsigned") {
    auto dt = std::make_shared<IntegerDatatype>(8, false);
    hid_t h5t = convertToHdf5(dt);
    CHECK(H5Tget_class(h5t) == H5T_INTEGER);
    CHECK(H5Tget_size(h5t) == 8);
    CHECK(H5Tget_sign(h5t) == H5T_SGN_NONE);
  }
}

TEST_CASE("convertToHdf5 ArrayDatatype" * doctest::test_suite("io")) {
  auto base = std::make_shared<F64Datatype>();
  auto dt = std::make_shared<ArrayDatatype>(base, std::vector<std::size_t>{3, 5});
  hid_t h5t = convertToHdf5(dt);
  CHECK(H5Tget_class(h5t) == H5T_ARRAY);
  CHECK(H5Tget_size(h5t) == 3 * 5 * 8);
  int ndims = H5Tget_array_ndims(h5t);
  CHECK(ndims == 2);
  std::vector<hsize_t> dims(ndims);
  H5Tget_array_dims(h5t, dims.data());
  CHECK(dims[0] == 3);
  CHECK(dims[1] == 5);
  H5Tclose(h5t);
}

TEST_CASE("convertToHdf5 StructDatatype" * doctest::test_suite("io")) {
  auto f64 = std::make_shared<F64Datatype>();
  auto f32 = std::make_shared<F32Datatype>();
  std::vector<StructDatatype::MemberInfo> members = {
      {"x", 0, f64},
      {"y", 8, f32},
  };
  auto dt = std::make_shared<StructDatatype>(members, 16);
  hid_t h5t = convertToHdf5(dt);
  CHECK(H5Tget_class(h5t) == H5T_COMPOUND);
  CHECK(H5Tget_nmembers(h5t) == 2);
  CHECK(H5Tget_size(h5t) == 16);

  // Check member names
  char* name0 = H5Tget_member_name(h5t, 0);
  CHECK(std::string(name0) == "x");
  H5free_memory(name0);

  char* name1 = H5Tget_member_name(h5t, 1);
  CHECK(std::string(name1) == "y");
  H5free_memory(name1);

  // Check offsets
  CHECK(H5Tget_member_offset(h5t, 0) == 0);
  CHECK(H5Tget_member_offset(h5t, 1) == 8);

  H5Tclose(h5t);
}

TEST_CASE("convertToHdf5 OpaqueDatatype" * doctest::test_suite("io")) {
  auto dt = std::make_shared<OpaqueDatatype>(42);
  hid_t h5t = convertToHdf5(dt);
  CHECK(H5Tget_class(h5t) == H5T_OPAQUE);
  CHECK(H5Tget_size(h5t) == 42);
  H5Tclose(h5t);
}

TEST_CASE("convertToHdf5 StringDatatype" * doctest::test_suite("io")) {
  auto dt = std::make_shared<StringDatatype>(64);
  hid_t h5t = convertToHdf5(dt);
  CHECK(H5Tget_class(h5t) == H5T_STRING);
  CHECK(H5Tget_size(h5t) == 64);
  H5Tclose(h5t);
}

} // namespace seissol::unit_test
