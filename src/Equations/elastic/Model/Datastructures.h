// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_DATASTRUCTURES_H_

#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::model {
class ElasticLocalData;
class ElasticNeighborData;

struct ElasticMaterial : Material {
  static constexpr std::size_t NumQuantities = 9;
  static constexpr std::size_t NumElasticQuantities = 9;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Elastic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "elastic";
  static inline const std::array<std::string, NumQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

  using LocalSpecificData = ElasticLocalData;
  using NeighborSpecificData = ElasticNeighborData;

  double lambda;
  double mu;

  [[nodiscard]] double getLambdaBar() const override { return lambda; }

  [[nodiscard]] double getMuBar() const override { return mu; }

  ElasticMaterial() = default;
  ElasticMaterial(const std::vector<double>& materialValues)
      : Material(materialValues), lambda(materialValues.at(2)), mu(materialValues.at(1)) {}

  ~ElasticMaterial() override = default;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {

    auto stiffnessTensorView =
        seissol_general::init::stiffnessTensor::view::create(fullTensor.data());
    stiffnessTensorView.setZero();
    stiffnessTensorView(0, 0, 0, 0) = lambda + 2 * mu;
    stiffnessTensorView(0, 0, 1, 1) = lambda;
    stiffnessTensorView(0, 0, 2, 2) = lambda;
    stiffnessTensorView(0, 1, 0, 1) = mu;
    stiffnessTensorView(0, 1, 1, 0) = mu;
    stiffnessTensorView(0, 2, 0, 2) = mu;
    stiffnessTensorView(0, 2, 2, 0) = mu;
    stiffnessTensorView(1, 0, 0, 1) = mu;
    stiffnessTensorView(1, 0, 1, 0) = mu;
    stiffnessTensorView(1, 1, 0, 0) = lambda;
    stiffnessTensorView(1, 1, 1, 1) = lambda + 2 * mu;
    stiffnessTensorView(1, 1, 2, 2) = lambda;
    stiffnessTensorView(1, 2, 1, 2) = mu;
    stiffnessTensorView(1, 2, 2, 1) = mu;
    stiffnessTensorView(2, 0, 0, 2) = mu;
    stiffnessTensorView(2, 0, 2, 0) = mu;
    stiffnessTensorView(2, 1, 2, 1) = mu;
    stiffnessTensorView(2, 2, 0, 0) = lambda;
    stiffnessTensorView(2, 2, 1, 1) = lambda;
    stiffnessTensorView(2, 2, 2, 2) = lambda + 2 * mu;
  }

  [[nodiscard]] double getMaxWaveSpeed() const override { return getPWaveSpeed(); }

  [[nodiscard]] double getPWaveSpeed() const override { return std::sqrt((lambda + 2 * mu) / rho); }

  [[nodiscard]] double getSWaveSpeed() const override { return std::sqrt(mu / rho); }

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_DATASTRUCTURES_H_
