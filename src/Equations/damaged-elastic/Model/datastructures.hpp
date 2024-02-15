#ifndef MODEL_DAMAGED_DATASTRUCTURES_H_
#define MODEL_DAMAGED_DATASTRUCTURES_H_

#include <Model/common_datastructures.hpp>
#include <cmath>
#include <generated_code/kernel.h>
#include <generated_code/init.h>

#include <iostream>

namespace seissol {
namespace model {
struct DamagedElasticMaterial : Material {
  double lambda0;
  double mu0;
  double gammaR;
  double xi0;
  double lambda; // effective
  double mu;     // effective
  double gamma;  // effective
  double epsxx_alpha;
  double epsyy_alpha;
  double epszz_alpha;
  double epsxy_alpha;
  double epsyz_alpha;
  double epszx_alpha;
  double Cd;

  DamagedElasticMaterial(){};
  DamagedElasticMaterial(double* materialValues, int numMaterialValues) {
    // assert(numMaterialValues == 3);

    this->rho = materialValues[0];
    this->mu0 = materialValues[1];
    this->lambda0 = materialValues[2];
    this->gammaR = materialValues[3];
    this->xi0 = materialValues[4];
    this->mu = materialValues[5];
    this->lambda = materialValues[6];
    this->gamma = materialValues[7];
    this->epsxx_alpha = materialValues[8];
    this->epsyy_alpha = materialValues[9];
    this->epszz_alpha = materialValues[10];
    this->epsxy_alpha = materialValues[11];
    this->epsyz_alpha = materialValues[12];
    this->epszx_alpha = materialValues[13];
    this->Cd = materialValues[14];
  }

  virtual ~DamagedElasticMaterial(){};

  void getFullStiffnessTensor(std::array<real, 81>& fullTensor) const final {

    auto stiffnessTensorView = init::stiffnessTensor::view::create(fullTensor.data());
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

  double getMaxWaveSpeed() const final { return getPWaveSpeed(); }

  double getPWaveSpeed() const final { return std::sqrt((lambda + 2 * mu) / rho); }

  double getSWaveSpeed() const final { return std::sqrt(mu / rho); }

  MaterialType getMaterialType() const override { return MaterialType::damaged; }
};
} // namespace model
} // namespace seissol

#endif
