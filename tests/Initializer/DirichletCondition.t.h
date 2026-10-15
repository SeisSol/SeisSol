// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/ParameterDB.h"
#include "Model/ScriptParameters.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace seissol::unit_test::dirichlet_condition {

/// A model written to a file of its own, removed again when done.
class ModelFile {
  public:
  ModelFile(const std::string& name, const std::string& content)
      : path_(std::filesystem::temp_directory_path() /
              ("seissol-dirichlet-" + std::to_string(std::random_device{}()) + "-" + name)) {
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

TEST_CASE("DirichletCondition: a script reads the group and the material of the cell" *
          doctest::test_suite("initializer")) {
  using Cfg = Config;
  using MaterialT = model::MaterialOf<Cfg>;
  const auto names = model::scriptParameterNames<MaterialT>();
  const std::string& first = names.front();
  const std::string& last = names.back();
  const std::string& velocity = MaterialT::Quantities[MaterialT::VelocityOffset];

  const ModelFile file("boundary.sderiv",
                       "out def const_" + velocity + " = 0.5 * group + " + first + "\n" +
                           "out def map_" + velocity + "_" + velocity + " = -" + last + "\n");
  const initializer::DirichletCondition condition("sderiv:" + file.path());

  std::vector<double> parameters(names.size());
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    parameters[i] = 1.0 + static_cast<double>(i);
  }
  const std::array<double, 3> barycenter{0.1, 0.2, 0.3};
  std::array<Real<Cfg>, tensor::dirichletMapGlobal<Cfg>::size()> map{};
  std::array<Real<Cfg>, tensor::dirichletOffsetGlobal<Cfg>::size()> offset{};

  // the same model for two faces of other cells
  for (const int group : {3, 7}) {
    CAPTURE(group);
    const auto frame = condition.query<Cfg>(
        barycenter.data(), group, parameters.data(), map.data(), offset.data());
    CHECK(frame == initializer::BoundaryFrame::Global);
    auto mapView = init::dirichletMapGlobal<Cfg>::view::create(map.data());
    auto offsetView = init::dirichletOffsetGlobal<Cfg>::view::create(offset.data());
    const auto v = MaterialT::VelocityOffset;
    CHECK(multisim::multisimWrap<Cfg>(offsetView, 0, v) ==
          static_cast<Real<Cfg>>(0.5 * group + parameters.front()));
    CHECK(mapView(v, v) == static_cast<Real<Cfg>>(-parameters.back()));
    // what the model does not give: the identity, and no offset
    CHECK(mapView(0, 0) == 1);
    CHECK(multisim::multisimWrap<Cfg>(offsetView, 0, 0) == 0);
    std::for_each(parameters.begin(), parameters.end(), [](double& value) { value *= 2.0; });
  }
}

} // namespace seissol::unit_test::dirichlet_condition
