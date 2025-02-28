// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_DATASTRUCTURES_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include <array>
#include <cstddef>
#include <string>

namespace seissol::model {
class AnisotropicLocalData;
class AnisotropicNeighborData;

struct AnisotropicMaterial : public Material {
  static constexpr std::size_t NumQuantities = 9;
  static constexpr std::size_t NumElasticQuantities = 9;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Anisotropic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "anisotropic";
  static inline const std::array<std::string, NumQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

  using LocalSpecificData = AnisotropicLocalData;
  using NeighborSpecificData = AnisotropicNeighborData;

  double c11;
  double c12;
  double c13;
  double c14;
  double c15;
  double c16;
  double c22;
  double c23;
  double c24;
  double c25;
  double c26;
  double c33;
  double c34;
  double c35;
  double c36;
  double c44;
  double c45;
  double c46;
  double c55;
  double c56;
  double c66;

  [[nodiscard]] double getLambdaBar() const override;

  [[nodiscard]] double getMuBar() const override;

  AnisotropicMaterial();

  explicit AnisotropicMaterial(const ElasticMaterial& m);

  AnisotropicMaterial(const std::vector<double>& materialValues);

  ~AnisotropicMaterial() override;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override;

  // calculate maximal wave speed
  [[nodiscard]] double getMaxWaveSpeed() const override;

  // calculate P-wave speed based on averaged material parameters
  [[nodiscard]] double getPWaveSpeed() const override;

  // calculate S-wave speed based on averaged material parameters
  [[nodiscard]] double getSWaveSpeed() const override;

  [[nodiscard]] MaterialType getMaterialType() const override;
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_DATASTRUCTURES_H_
