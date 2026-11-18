// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/MPI.h"
#include "Reader/Datafield/CellData.h"
#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <filesystem>
#include <hdf5.h>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace seissol::unit_test::celldata {

using seissol::reader::datafield::CellDataStore;
using seissol::reader::datafield::GridDesc;
using seissol::reader::datafield::GridKind;

/// An HDF5 file in a directory of its own, removed again when done.
class CellFile {
  public:
  CellFile()
      : dir_(std::filesystem::temp_directory_path() /
             ("seissol-celldata-" + std::to_string(std::random_device{}()))) {
    std::filesystem::create_directories(dir_);
  }
  ~CellFile() {
    std::error_code error;
    std::filesystem::remove_all(dir_, error);
  }
  CellFile(const CellFile&) = delete;
  CellFile& operator=(const CellFile&) = delete;
  CellFile(CellFile&&) = delete;
  CellFile& operator=(CellFile&&) = delete;

  [[nodiscard]] std::string path() const { return (dir_ / "cells.h5").string(); }

  /// Writes the dataset `name` (in its groups, which it makes) of `rows` rows of `width` values,
  /// the value of row r and column c being 10 r + c; of one value per row without a second axis
  /// where `width` is 0. In binary32 where `single`.
  void write(const std::string& name, hsize_t rows, hsize_t width, bool single = false) const {
    const hid_t file = H5Fcreate(path().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t links = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(links, 1);
    const hsize_t columns = width == 0 ? 1 : width;
    const std::vector<hsize_t> dims =
        width == 0 ? std::vector<hsize_t>{rows} : std::vector<hsize_t>{rows, width};
    const hid_t space = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    std::vector<double> values(rows * columns);
    for (hsize_t r = 0; r < rows; ++r) {
      for (hsize_t c = 0; c < columns; ++c) {
        values[r * columns + c] = 10.0 * static_cast<double>(r) + static_cast<double>(c);
      }
    }
    const hid_t dataset = H5Dcreate2(file,
                                     name.c_str(),
                                     single ? H5T_IEEE_F32LE : H5T_IEEE_F64LE,
                                     space,
                                     links,
                                     H5P_DEFAULT,
                                     H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
    H5Dclose(dataset);
    H5Sclose(space);
    H5Pclose(links);
    H5Fclose(file);
  }

  /// The declaration `d` of the variable `variable` of the file, of `components`.
  [[nodiscard]] GridDesc desc(const std::string& variable,
                              const std::vector<std::string>& components) const {
    GridDesc desc;
    desc.kind = GridKind::Cells;
    desc.name = "d";
    desc.path = path();
    desc.variable = variable;
    desc.componentNames = components;
    return desc;
  }

  private:
  std::filesystem::path dir_;
};

/// The message of what `load` throws, empty where it throws nothing.
template <typename F>
std::string failureOf(F&& load) {
  try {
    load();
  } catch (const std::runtime_error& error) {
    return error.what();
  }
  return "";
}

TEST_CASE("Datafield cells: the rows of the cells of a rank, by their ids in the mesh file") {
  const CellFile file;
  file.write("material", 5, 2);
  CellDataStore store;
  const GridDesc desc = file.desc("material", {"rho", "vs"});
  CHECK(!store.loaded(desc));
  // the cells of this rank in another order than the file's, one of them twice
  store.load(desc, {4, 0, 2, 1, 3, 2}, seissol::Mpi::mpi.comm());
  REQUIRE(store.loaded(desc));
  for (std::size_t id = 0; id < 5; ++id) {
    const double* row = store.row(desc, id);
    CHECK(row[0] == 10.0 * static_cast<double>(id));
    CHECK(row[1] == 10.0 * static_cast<double>(id) + 1.0);
  }
  // a later read replaces the rows: the cells after the partition those before it
  store.load(desc, {3, 1, 4, 0, 2}, seissol::Mpi::mpi.comm());
  CHECK(store.row(desc, 3)[1] == 31.0);
}

TEST_CASE("Datafield cells: a value per cell in a group, in binary32") {
  const CellFile file;
  file.write("fields/rho", 3, 0, true);
  CellDataStore store;
  const GridDesc desc = file.desc("fields/rho", {"rho"});
  store.load(desc, {0, 1, 2}, seissol::Mpi::mpi.comm());
  CHECK(store.row(desc, 2)[0] == 20.0);
  CHECK(failureOf([&] {
          static_cast<void>(store.row(file.desc("fields/vs", {"vs"}), 0));
        }).find("was not read") != std::string::npos);
}

TEST_CASE("Datafield cells: a dataset that does not fit the mesh or the declaration is an error") {
  const CellFile file;
  file.write("material", 5, 2);
  CellDataStore store;
  const auto comm = seissol::Mpi::mpi.comm();
  CHECK(failureOf([&] {
          store.load(file.desc("material", {"rho", "vs", "vp"}), {0, 1, 2, 3, 4}, comm);
        }).find("a cell has 2 values in the file, but the declaration names 3 components") !=
        std::string::npos);
  CHECK(failureOf([&] {
          store.load(file.desc("material", {"rho", "vs"}), {0, 1, 2, 3}, comm);
        }).find("it has 5 rows, but the mesh has 4 cells") != std::string::npos);
  CHECK(failureOf([&] {
          store.load(file.desc("other", {"rho"}), {0, 1, 2, 3, 4}, comm);
        }).find("the file has no other") != std::string::npos);
  GridDesc missing = file.desc("material", {"rho", "vs"});
  missing.path += ".gone";
  CHECK(failureOf([&] {
          store.load(missing, {0, 1, 2, 3, 4}, comm);
        }).find("there is no such file") != std::string::npos);
}

} // namespace seissol::unit_test::celldata
