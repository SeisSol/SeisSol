#ifndef MODEL_POROELASTIC_DATASTRUCTURES_H_
#define MODEL_POROELASTIC_DATASTRUCTURES_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "Model/CommonDatastructures.h"
#include <array>
#include <cassert>
#include <cstddef>
#include <string>

namespace seissol::model {
struct PoroElasticMaterial : ElasticMaterial {
  static constexpr std::size_t NumQuantities = 13;
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

  double bulkSolid;
  double porosity;
  double permeability;
  double tortuosity;
  double bulkFluid;
  double rhoFluid;
  double viscosity;

  PoroElasticMaterial() = default;

  PoroElasticMaterial(const double* materialValues, int numMaterialValues) {
    assert(numMaterialValues == 10);

    this->bulkSolid = materialValues[0];
    this->rho = materialValues[1];
    this->lambda = materialValues[2];
    this->mu = materialValues[3];
    this->porosity = materialValues[4];
    this->permeability = materialValues[5];
    this->tortuosity = materialValues[6];
    this->bulkFluid = materialValues[7];
    this->rhoFluid = materialValues[8];
    this->viscosity = materialValues[9];
  };
  ~PoroElasticMaterial() override = default;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {
    double elasticMaterialVals[] = {this->rho, this->mu, this->lambda};
    ElasticMaterial em(elasticMaterialVals, 3);
    em.getFullStiffnessTensor(fullTensor);
  }

  double getMaxWaveSpeed() const override { return getPWaveSpeed(); }

  // only declare it here and define in a separate datastructures.cpp
  // to circumvent problems with circular includes
  double getPWaveSpeed() const override;

  double getSWaveSpeed() const override { return std::sqrt(mu / rho); }

  MaterialType getMaterialType() const override { return Type; }
};
} // namespace seissol::model

#endif
