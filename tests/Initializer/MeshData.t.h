// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Equations/elastic/Model/Datastructures.h"
#include "Initializer/ParameterDB.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <hdf5.h>
#include <random>
#include <string>
#include <system_error>
#include <vector>

namespace seissol::unit_test::mesh_material {

/// A directory of its own, removed again when done.
class Directory {
  public:
  Directory()
      : path_(std::filesystem::temp_directory_path() /
              ("seissol-meshdata-" + std::to_string(std::random_device{}()))) {
    std::filesystem::create_directories(path_);
  }
  ~Directory() {
    std::error_code error;
    std::filesystem::remove_all(path_, error);
  }
  Directory(const Directory&) = delete;
  Directory& operator=(const Directory&) = delete;
  Directory(Directory&&) = delete;
  Directory& operator=(Directory&&) = delete;

  [[nodiscard]] std::string file(const std::string& name) const { return (path_ / name).string(); }

  private:
  std::filesystem::path path_;
};

/// The material of the row `row` of the file: rho, mu and lambda.
std::array<double, 3> materialOfRow(std::size_t row) {
  const auto r = static_cast<double>(row);
  return {1000.0 + r, (r + 1.0) * 1.0e9, (r + 2.0) * 1.0e9};
}

/// The dataset `material` of `rows` cells, a row of rho, mu and lambda each (materialOfRow()).
void writeMaterial(const std::string& path, std::size_t rows) {
  const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  const std::array<hsize_t, 2> dims{rows, 3};
  const hid_t space = H5Screate_simple(2, dims.data(), nullptr);
  std::vector<double> values;
  for (std::size_t row = 0; row < rows; ++row) {
    const auto material = materialOfRow(row);
    values.insert(values.end(), material.begin(), material.end());
  }
  const hid_t dataset =
      H5Dcreate2(file, "material", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
  H5Dclose(dataset);
  H5Sclose(space);
  H5Fclose(file);
}

TEST_CASE("Materials: data on the mesh gives the material of every cell" *
          doctest::test_suite("initializer")) {
  const Directory directory;
  writeMaterial(directory.file("cells.h5"), 3);
  const std::string script = directory.file("material.sderiv");
  std::ofstream(script)
      << "dataset m (kind = \"cells\", path = \"cells.h5\", variable = \"material\",\n"
         "           components = (\"rho\", \"mu\", \"lambda\"))\n"
         "out def rho = m.rho\n"
         "out def mu = m.mu\n"
         "out def lambda = select(lt(z, 1e9), m.lambda, 0)\n";

  // three cells of this rank, the rows 2, 0 and 1 of the file
  const std::vector<std::array<std::array<double, 3>, 4>> corners{
      {{{0.6, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.6, 1.0, 0.0}, {0.6, 0.0, 1.0}}},
      {{{0.4, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.6, 1.0, 0.0}, {0.6, 0.0, 1.0}}},
      {{{0.0, 0.0, 0.0}, {0.4, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}}};
  auto cells = initializer::CellToVertexArray::fromVectors(corners, {1, 1, 1});
  const std::vector<std::uint64_t> ids{2, 0, 1};
  cells.elementIds = [ids](std::size_t cell) { return ids[cell]; };
  initializer::loadCellData({script}, cells);

  // at the barycenters
  std::vector<model::ElasticMaterial> materials;
  initializer::MaterialParameterDB<model::ElasticMaterial> database;
  database.setMaterialVector(&materials);
  database.evaluateModel(script, initializer::ElementBarycenterGenerator(cells));
  REQUIRE(materials.size() == 3);
  for (std::size_t cell = 0; cell < 3; ++cell) {
    const auto expected = materialOfRow(ids[cell]);
    CHECK(materials[cell].rho == expected[0]);
    CHECK(materials[cell].mu == expected[1]);
    CHECK(materials[cell].lambda == expected[2]);
  }

  // at the points of the quadrature of a cell, its value at every one: the average is that value
  std::vector<model::ElasticMaterial> averaged;
  initializer::MaterialParameterDB<model::ElasticMaterial> averaging;
  averaging.setMaterialVector(&averaged);
  averaging.evaluateModel(script, initializer::ElementAverageGenerator(cells, 3));
  REQUIRE(averaged.size() == 3);
  for (std::size_t cell = 0; cell < 3; ++cell) {
    const auto expected = materialOfRow(ids[cell]);
    CHECK(averaged[cell].rho == doctest::Approx(expected[0]).epsilon(1e-14));
    CHECK(averaged[cell].mu == doctest::Approx(expected[1]).epsilon(1e-14));
    CHECK(averaged[cell].lambda == doctest::Approx(expected[2]).epsilon(1e-14));
  }
}

} // namespace seissol::unit_test::mesh_material
