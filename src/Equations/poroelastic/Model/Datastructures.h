// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_DATASTRUCTURES_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "Model/CommonDatastructures.h"
#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::model {
class PoroelasticLocalData;
class PoroelasticNeighborData;

struct PoroElasticMaterial : public ElasticMaterial {
  static constexpr std::size_t NumQuantities = 13;
  static constexpr std::size_t NumElasticQuantities = 13;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Poroelastic;
  static constexpr LocalSolver Solver = LocalSolver::SpaceTimePredictorPoroelastic;
  static inline const std::string Text = "poroelastic";
  static inline const std::array<std::string, NumQuantities> Quantities{"s_xx",
                                                                        "s_yy",
                                                                        "s_zz",
                                                                        "s_xy",
                                                                        "s_yz",
                                                                        "s_xz",
                                                                        "v1",
                                                                        "v2",
                                                                        "v3",
                                                                        "p",
                                                                        "v1_f",
                                                                        "v2_f",
                                                                        "v3_f"};

  using LocalSpecificData = PoroelasticLocalData;
  using NeighborSpecificData = PoroelasticNeighborData;

  double bulkSolid;
  double porosity;
  double permeability;
  double tortuosity;
  double bulkFluid;
  double rhoFluid;
  double viscosity;

  PoroElasticMaterial() = default;

  PoroElasticMaterial(const std::vector<double>& materialValues)
      : bulkSolid(materialValues.at(0)), porosity(materialValues.at(4)),
        permeability(materialValues.at(5)), tortuosity(materialValues.at(6)),
        bulkFluid(materialValues.at(7)), rhoFluid(materialValues.at(8)),
        viscosity(materialValues.at(9)) {
    // those are shifted right now by 1 compared to everywhere else;
    // hence we cannot use the base class initializer here
    rho = materialValues.at(1);
    lambda = materialValues.at(2);
    mu = materialValues.at(3);
  }
  ~PoroElasticMaterial() override = default;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {
    const std::vector<double> elasticMaterialVals{this->rho, this->mu, this->lambda};
    const ElasticMaterial em(elasticMaterialVals);
    em.getFullStiffnessTensor(fullTensor);
  }

  [[nodiscard]] double getMaxWaveSpeed() const override { return getPWaveSpeed(); }

  // only declare it here and define in a separate datastructures.cpp
  // to circumvent problems with circular includes
  [[nodiscard]] double getPWaveSpeed() const override;

  [[nodiscard]] double getSWaveSpeed() const override { return std::sqrt(mu / rho); }

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_DATASTRUCTURES_H_
