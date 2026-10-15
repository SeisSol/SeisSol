// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/ParameterDB.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace seissol::unit_test::output_regions {

using initializer::OutputRegions;

/// A model written to a file of its own, removed again when done.
class ModelFile {
  public:
  ModelFile(const std::string& name, const std::string& content)
      : path_(std::filesystem::temp_directory_path() /
              ("seissol-outputregions-" + std::to_string(std::random_device{}()) + "-" + name)) {
    std::ofstream(path_) << content;
  }
  ~ModelFile() { std::filesystem::remove(path_); }
  ModelFile(const ModelFile&) = delete;
  ModelFile& operator=(const ModelFile&) = delete;
  ModelFile(ModelFile&&) = delete;
  ModelFile& operator=(ModelFile&&) = delete;

  [[nodiscard]] std::string path() const { return path_.string(); }

  private:
  std::filesystem::path path_;
};

/// Three tetrahedra along x: the first right of x = 0.5, the second across it, the third left of
/// it; in the groups 1, 2, 2, of elastic materials with the densities 1, 2, 3.
struct Cells {
  std::vector<std::array<std::array<double, 3>, 4>> corners{
      {{{0.6, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.6, 1.0, 0.0}, {0.6, 0.0, 1.0}}},
      {{{0.4, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.6, 1.0, 0.0}, {0.6, 0.0, 1.0}}},
      {{{0.0, 0.0, 0.0}, {0.4, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}}};
  std::vector<int> groups{1, 2, 2};
  std::vector<double> densities{1.0, 2.0, 3.0};

  [[nodiscard]] OutputRegions::Parameters parameters() const {
    OutputRegions::Parameters parameters;
    parameters.names = {"rho", "mu", "lambda"};
    parameters.values = [this](const std::string& name) {
      return name == "rho"
                 ? densities
                 : std::vector<double>(densities.size(), std::numeric_limits<double>::quiet_NaN());
    };
    return parameters;
  }

  [[nodiscard]] std::vector<bool> select(const OutputRegions& regions,
                                         const std::string& name) const {
    return regions.select(
        name,
        corners.size(),
        4,
        [this](std::size_t cell, std::size_t corner) { return corners[cell][corner]; },
        [this](std::size_t cell) { return groups[cell]; },
        parameters());
  }
};

TEST_CASE("OutputRegions: an element is in the region where the model is positive at a corner" *
          doctest::test_suite("initializer")) {
  const ModelFile file("region.sderiv",
                       "out def wavefield = 0.5 - x\n"
                       "out def surface = select(eq(group, 2.0), 1.0, 0.0)\n");
  const OutputRegions regions("sderiv:" + file.path());
  CHECK(regions.restricts(OutputRegions::WaveField));
  CHECK(regions.restricts(OutputRegions::Surface));

  const Cells cells;
  CHECK(cells.select(regions, OutputRegions::WaveField) == std::vector<bool>{false, true, true});
  CHECK(cells.select(regions, OutputRegions::Surface) == std::vector<bool>{false, true, true});
}

TEST_CASE("OutputRegions: a model reads the parameters of the material of a cell" *
          doctest::test_suite("initializer")) {
  const ModelFile file("region.sderiv",
                       "out def wavefield = rho - 1.5\n"
                       "out def surface = select(eq(group, 2.0), 2.5 - rho, 0.0)\n");
  const OutputRegions regions("sderiv:" + file.path());
  const Cells cells;
  CHECK(cells.select(regions, OutputRegions::WaveField) == std::vector<bool>{false, true, true});
  CHECK(cells.select(regions, OutputRegions::Surface) == std::vector<bool>{false, true, false});
}

TEST_CASE("OutputRegions: an output the model does not name is not restricted" *
          doctest::test_suite("initializer")) {
  const ModelFile file("region.lua",
                       "local M = {}\n"
                       "M.output_parameters = {\"wavefield\"}\n"
                       "M.input_parameters = {\"x\", \"y\", \"z\"}\n"
                       "function M.evaluate(fields, x, y, z)\n"
                       "  if x < 0.5 then\n"
                       "    return 1.0\n"
                       "  end\n"
                       "  return 0.0\n"
                       "end\n"
                       "return M\n");
  const OutputRegions regions("lua:" + file.path());
  CHECK(regions.restricts(OutputRegions::WaveField));
  CHECK(!regions.restricts(OutputRegions::Surface));

  const Cells cells;
  CHECK(cells.select(regions, OutputRegions::WaveField) == std::vector<bool>{false, true, true});
  CHECK(cells.select(regions, OutputRegions::Surface) == std::vector<bool>{true, true, true});
}

TEST_CASE("OutputRegions: without a model, nothing is restricted" *
          doctest::test_suite("initializer")) {
  const OutputRegions regions;
  CHECK(!regions.restricts(OutputRegions::WaveField));
  const Cells cells;
  CHECK(cells.select(regions, OutputRegions::WaveField) == std::vector<bool>{true, true, true});
  // and nothing to select from
  CHECK(regions
            .select(
                OutputRegions::WaveField,
                0,
                4,
                [](std::size_t, std::size_t) { return std::array<double, 3>{}; },
                [](std::size_t) { return 0; },
                cells.parameters())
            .empty());
}

} // namespace seissol::unit_test::output_regions
