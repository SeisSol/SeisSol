// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_DATASTRUCTURES_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/LinearCK/Solver.h"
#include "Model/CommonDatastructures.h"

#include <array>
#include <cstddef>
#include <string>

namespace seissol::model {
struct AnisotropicMaterial : public Material {
  static constexpr std::size_t NumQuantities = 9;
  static constexpr std::size_t NumElasticQuantities = 9;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t TractionQuantities = 6;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Anisotropic;
  static inline const std::string Text = "anisotropic";
  static inline const std::array<std::string, NumQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};
  static constexpr std::size_t Parameters = 21 + Material::Parameters;

  static constexpr bool SupportsDR = false;
  static constexpr bool SupportsLTS = true;

  template <typename Cfg>
  using LocalSpecificData = std::monostate;

  template <typename Cfg>
  using NeighborSpecificData = std::monostate;

  template <typename Cfg>
  using Solver = kernels::solver::linearck::Solver;

  double c11{};
  double c12{};
  double c13{};
  double c14{};
  double c15{};
  double c16{};
  double c22{};
  double c23{};
  double c24{};
  double c25{};
  double c26{};
  double c33{};
  double c34{};
  double c35{};
  double c36{};
  double c44{};
  double c45{};
  double c46{};
  double c55{};
  double c56{};
  double c66{};

  static const std::unordered_map<std::string, double AnisotropicMaterial::*> ParameterMap;

  [[nodiscard]] double getLambdaBar() const override;

  [[nodiscard]] double getMuBar() const override;

  AnisotropicMaterial();

  explicit AnisotropicMaterial(const ElasticMaterial& m);

  explicit AnisotropicMaterial(const std::vector<double>& materialValues);

  ~AnisotropicMaterial() override;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override;

  // calculate maximal wave speed
  [[nodiscard]] double getMaxWaveSpeed() const override;

  // calculate P-wave speed based on averaged material parameters
  [[nodiscard]] double getPWaveSpeed() const override;

  // calculate S-wave speed based on averaged material parameters
  [[nodiscard]] double getSWaveSpeed() const override;

  [[nodiscard]] MaterialType getMaterialType() const override;

  void setLameParameters(double mu, double lambda) override;
};

inline const std::unordered_map<std::string, double AnisotropicMaterial::*>
    AnisotropicMaterial::ParameterMap{
        {"rho", &AnisotropicMaterial::rho}, {"c11", &AnisotropicMaterial::c11},
        {"c12", &AnisotropicMaterial::c12}, {"c13", &AnisotropicMaterial::c13},
        {"c14", &AnisotropicMaterial::c14}, {"c15", &AnisotropicMaterial::c15},
        {"c16", &AnisotropicMaterial::c16}, {"c22", &AnisotropicMaterial::c22},
        {"c23", &AnisotropicMaterial::c23}, {"c24", &AnisotropicMaterial::c24},
        {"c25", &AnisotropicMaterial::c25}, {"c26", &AnisotropicMaterial::c26},
        {"c33", &AnisotropicMaterial::c33}, {"c34", &AnisotropicMaterial::c34},
        {"c35", &AnisotropicMaterial::c35}, {"c36", &AnisotropicMaterial::c36},
        {"c44", &AnisotropicMaterial::c44}, {"c45", &AnisotropicMaterial::c45},
        {"c46", &AnisotropicMaterial::c46}, {"c55", &AnisotropicMaterial::c55},
        {"c56", &AnisotropicMaterial::c56}, {"c66", &AnisotropicMaterial::c66},
    };

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_DATASTRUCTURES_H_
