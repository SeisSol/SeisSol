// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf
// SPDX-FileContributor: Jinwen Pan

#ifndef SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_DATASTRUCTURES_H_

#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::model {
class AcousticLocalData;
class AcousticNeighborData;

struct AcousticMaterial : public Material {
  static constexpr std::size_t NumQuantities = 4;
  static constexpr std::size_t NumElasticQuantities = 4;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Acoustic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "acoustic";
  // The stress-velocity formulation of the elastic model is reused.
  // By definition, the normal stress and pressure are negatives of each other.
  static inline const std::array<std::string, NumQuantities> Quantities = {"-p", "v1", "v2", "v3"};

  using LocalSpecificData = AcousticLocalData;
  using NeighborSpecificData = AcousticNeighborData;

  double lambda;

  [[nodiscard]] double getLambdaBar() const override { return lambda; }

  [[nodiscard]] double getMuBar() const override { return 0.0; }

  AcousticMaterial() = default;
  AcousticMaterial(const std::vector<double>& materialValues)
      : Material(materialValues), lambda(materialValues.at(1)) {}

  ~AcousticMaterial() override = default;

  // The stiffness tensor of the elastic model is reused.
  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {

    auto stiffnessTensorView =
        seissol_general::init::stiffnessTensor::view::create(fullTensor.data());
    stiffnessTensorView.setZero();
    stiffnessTensorView(0, 0, 0, 0) = lambda;
    stiffnessTensorView(0, 0, 1, 1) = lambda;
    stiffnessTensorView(0, 0, 2, 2) = lambda;
    stiffnessTensorView(1, 1, 0, 0) = lambda;
    stiffnessTensorView(1, 1, 1, 1) = lambda;
    stiffnessTensorView(1, 1, 2, 2) = lambda;
    stiffnessTensorView(2, 2, 0, 0) = lambda;
    stiffnessTensorView(2, 2, 1, 1) = lambda;
    stiffnessTensorView(2, 2, 2, 2) = lambda;
  }

  [[nodiscard]] double getMaxWaveSpeed() const override { return getPWaveSpeed(); }

  [[nodiscard]] double getPWaveSpeed() const override { return std::sqrt(lambda / rho); }

  [[nodiscard]] double getSWaveSpeed() const override { return 0.0; }

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_DATASTRUCTURES_H_
