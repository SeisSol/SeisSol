// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_GEOMETRY_PUMLBOUNDARY_T_H_
#define SEISSOL_TESTS_GEOMETRY_PUMLBOUNDARY_T_H_

// The boundary condition tags of a mesh file may come in any integer type. PUMGen files them as
// signed integers, so that in i32 and i64, a tag of 128 resp. 32768 or more on the last face of a
// cell makes the value of the cell negative -- which HDF5 clips to 0 where it is read as an
// unsigned type. The files here are mesh.h5 with such a tag on the last face of every cell.

#include <doctest.h>

#include "Common/Constants.h"
#include "Geometry/PUMLReader.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Parallel/MPI.h"
#include "TestHelper.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <hdf5.h>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace pumlboundarytest {

using seissol::initializer::parameters::BoundaryFormat;

/// the tag of a face of a cell in the files here: the last face has one of the upper half of the
/// range of its format
inline std::uint32_t tagOf(BoundaryFormat format, std::size_t cell, std::size_t face) {
  const bool last = face + 1 == Cell::NumFaces;
  if (format == BoundaryFormat::I32) {
    return static_cast<std::uint32_t>(last ? 0xff - cell : (7 * cell + face) % 0x100);
  }
  if (format == BoundaryFormat::I64) {
    return static_cast<std::uint32_t>(last ? 0xffff - cell : (1031 * cell + face) % 0x10000);
  }
  return static_cast<std::uint32_t>(last ? 0x7fffffff - cell : 100000 * cell + face);
}

/// the values of /boundary, as their bits: one per cell that packs its tags, or one per face
inline std::vector<std::uint64_t> valuesOf(BoundaryFormat format, std::size_t cells) {
  if (format == BoundaryFormat::I32x4) {
    std::vector<std::uint64_t> values(cells * Cell::NumFaces);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        values[cell * Cell::NumFaces + face] = tagOf(format, cell, face);
      }
    }
    return values;
  }
  const std::size_t bits = format == BoundaryFormat::I32 ? 8 : 16;
  std::vector<std::uint64_t> values(cells);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      values[cell] |= static_cast<std::uint64_t>(tagOf(format, cell, face)) << (bits * face);
    }
  }
  return values;
}

/// mesh.h5 with the tags of tagOf in /boundary, filed as `fileType` (a little-endian integer type)
inline void writeMesh(const std::string& file, BoundaryFormat format, hid_t fileType) {
  if (seissol::Mpi::mpi.rank() == 0) {
    const hid_t source = H5Fopen(tpath("Testing/mesh.h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    const hid_t target = H5Fcreate(file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    for (const char* name : {"connect", "geometry", "group"}) {
      H5Ocopy(source, name, target, name, H5P_DEFAULT, H5P_DEFAULT);
    }

    const hid_t connect = H5Dopen2(source, "connect", H5P_DEFAULT);
    const hid_t connectSpace = H5Dget_space(connect);
    std::vector<hsize_t> connectDims(2);
    H5Sget_simple_extent_dims(connectSpace, connectDims.data(), nullptr);
    H5Sclose(connectSpace);
    H5Dclose(connect);
    const auto cells = static_cast<std::size_t>(connectDims[0]);

    // the bytes of the values in the file type, so that HDF5 writes them without converting them
    const auto values = valuesOf(format, cells);
    const auto width = H5Tget_size(fileType);
    std::vector<unsigned char> bytes(values.size() * width);
    for (std::size_t i = 0; i < values.size(); ++i) {
      for (std::size_t byte = 0; byte < width; ++byte) {
        bytes[i * width + byte] = static_cast<unsigned char>(values[i] >> (8 * byte));
      }
    }

    const bool perFace = format == BoundaryFormat::I32x4;
    const std::vector<hsize_t> dims{cells, Cell::NumFaces};
    const hid_t space = H5Screate_simple(perFace ? 2 : 1, dims.data(), nullptr);
    const hid_t dataset =
        H5Dcreate2(target, "boundary", fileType, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, fileType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bytes.data());
    H5Dclose(dataset);
    H5Sclose(space);
    H5Fclose(target);
    H5Fclose(source);
  }
  MPI_Barrier(seissol::Mpi::mpi.comm());
}

/// whether every face of every local cell of a mesh has the tag of tagOf
inline bool hasTheTags(const seissol::geometry::PumlMesh& mesh,
                       const seissol::geometry::PumlCellArrays& arrays) {
  const auto tags = seissol::geometry::boundaryTagsOf(mesh, arrays);
  const auto cellIds = mesh.data(arrays.cellIdsAsInFile);
  bool same = true;
  for (std::size_t cell = 0; cell < cellIds.size(); ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      same &= tags.tag(cell, face) == tagOf(arrays.boundaryFormat, cellIds[cell], face);
    }
  }
  return same;
}

} // namespace pumlboundarytest

TEST_CASE("PUMLReader reads the boundary condition tags whatever integer type the file has" *
          doctest::test_suite("geometry")) {
  using pumlboundarytest::BoundaryFormat;
  struct Case {
    std::string name;
    BoundaryFormat format;
    hid_t fileType;
  };
  // the first of every format as PUMGen writes it
  const std::vector<Case> cases = {
      {.name = "i32-int32", .format = BoundaryFormat::I32, .fileType = H5T_STD_I32LE},
      {.name = "i32-uint32", .format = BoundaryFormat::I32, .fileType = H5T_STD_U32LE},
      {.name = "i32-int64", .format = BoundaryFormat::I32, .fileType = H5T_STD_I64LE},
      {.name = "i64-int64", .format = BoundaryFormat::I64, .fileType = H5T_STD_I64LE},
      {.name = "i64-uint64", .format = BoundaryFormat::I64, .fileType = H5T_STD_U64LE},
      {.name = "i32x4-int32", .format = BoundaryFormat::I32x4, .fileType = H5T_STD_I32LE},
      {.name = "i32x4-uint32", .format = BoundaryFormat::I32x4, .fileType = H5T_STD_U32LE}};

  for (const auto& testCase : cases) {
    CAPTURE(testCase.name);
    const auto file = tpath("Testing/boundary-" + testCase.name + ".h5");
    pumlboundarytest::writeMesh(file, testCase.format, testCase.fileType);

    std::cout.setstate(std::ios_base::failbit);
    seissol::geometry::PumlMesh mesh;
    mesh.setComm(seissol::Mpi::mpi.comm());
    const auto arrays = seissol::geometry::PUMLReader::read(mesh, file, testCase.format);
    std::cout.clear();
    CHECK(pumlboundarytest::hasTheTags(mesh, arrays));

    // and the clone that a topology of its own starts from, through its own handles
    const auto clone = mesh.clone();
    CHECK(pumlboundarytest::hasTheTags(clone, arrays.rebasedTo(clone)));

    MPI_Barrier(seissol::Mpi::mpi.comm());
    if (seissol::Mpi::mpi.rank() == 0) {
      static_cast<void>(std::remove(file.c_str()));
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_GEOMETRY_PUMLBOUNDARY_T_H_
