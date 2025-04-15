// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf
// SPDX-FileContributor: Zihua Niu

#ifndef SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DATASTRUCTURES_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::model {
class DamageLocalData;
class DamageNeighborData;

struct DamageMaterial : Material {
  static constexpr std::size_t NumQuantities = 10;
  static constexpr std::size_t NumElasticQuantities = 10;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Damage;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "damage";
  static inline const std::array<std::string, NumQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3", "alpha"};

  using LocalSpecificData = DamageLocalData;
  using NeighborSpecificData = DamageNeighborData;

  double lambda;
  double mu;
  double gammaR;
  double Cd;      // damage evolution coefficient
  double lambdaE; // effective
  double muE;     // effective
  // store the initial strain for predictor step
  double epsInit_xx;
  double epsInit_yy;
  double epsInit_zz;
  double epsInit_xy;
  double epsInit_yz;
  double epsInit_xz;
  // track the total strain for predictor step
  double epsTot_xx;
  double epsTot_yy;
  double epsTot_zz;
  double epsTot_xy;
  double epsTot_yz;
  double epsTot_xz;

  [[nodiscard]] double getLambdaBar() const override { return lambda; }

  [[nodiscard]] double getMuBar() const override { return mu; }

  DamageMaterial() = default;

  void assignTotalStrain() override {
    this->muE = this->mu;
    this->lambdaE = this->lambda;
    this->epsTot_xx = this->epsInit_xx;
    this->epsTot_yy = this->epsInit_yy;
    this->epsTot_zz = this->epsInit_zz;
    this->epsTot_xy = this->epsInit_xy;
    this->epsTot_yz = this->epsInit_yz;
    this->epsTot_xz = this->epsInit_xz;
  }

  // This initialization is not used for material initialization with easi
  DamageMaterial(const std::vector<double>& materialValues)
      : Material(materialValues), lambda(materialValues.at(2)), mu(materialValues.at(1)) {
    this->gammaR = materialValues.at(3);
    this->Cd = materialValues.at(4);
    this->epsInit_xx = materialValues.at(5);
    this->epsInit_yy = materialValues.at(6);
    this->epsInit_zz = materialValues.at(7);
    this->epsInit_xy = materialValues.at(8);
    this->epsInit_yz = materialValues.at(9);
    this->epsInit_xz = materialValues.at(10);
    // Initialize the total strain also as initial strain
    // This is not useful at the init phase, how to init these values later?
    this->muE = materialValues.at(1);
    this->lambdaE = materialValues.at(2);
    this->epsTot_xx = materialValues.at(5);
    this->epsTot_yy = materialValues.at(6);
    this->epsTot_zz = materialValues.at(7);
    this->epsTot_xy = materialValues.at(8);
    this->epsTot_yz = materialValues.at(9);
    this->epsTot_xz = materialValues.at(10);
  }

  ~DamageMaterial() override = default;

  // Use the same elasticity tensor calculation as Elastic case
  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {
    const std::vector<double> elasticMaterialVals{this->rho, this->mu, this->lambda};
    const ElasticMaterial em(elasticMaterialVals);
    em.getFullStiffnessTensor(fullTensor);
  }

  [[nodiscard]] double getMaxWaveSpeed() const override { return getPWaveSpeed(); }

  [[nodiscard]] double getPWaveSpeed() const override { return std::sqrt((lambda + 2 * mu) / rho); }

  [[nodiscard]] double getSWaveSpeed() const override { return std::sqrt(mu / rho); }

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DATASTRUCTURES_H_
