#include "datastructures.hpp"
#include "PoroelasticSetup.h"
#include "Numerical_aux/Eigenvalues.h"

using namespace seissol::model;
PoroElasticMaterial::PoroElasticMaterial( double* i_materialVal, int i_numMaterialVals)
{ 
  assert(i_numMaterialVals == 10);

  this->bulk_solid = i_materialVal[0];
  this->rho = i_materialVal[1]; 
  this->lambda = i_materialVal[2];    
  this->mu = i_materialVal[3];
  this->porosity = i_materialVal[4]; 
  this->permeability = i_materialVal[5];
  this->tortuosity = i_materialVal[6];
  this->bulk_fluid = i_materialVal[7];
  this->rho_fluid = i_materialVal[8];
  this->viscosity = i_materialVal[9];  
}

void PoroElasticMaterial::getFullStiffnessTensor(std::array<real, 81>& fullTensor) const {
  double elasticMaterialVals[] = {this->rho, this->mu, this->lambda};
  ElasticMaterial em(elasticMaterialVals, 3);
  em.getFullStiffnessTensor(fullTensor);
}

double PoroElasticMaterial::getMaxWaveSpeed() const {
  return getPWaveSpeed();
}

double PoroElasticMaterial::getPWaveSpeed() const {
#ifdef HAS_ARMADILLO
  eigenvalues::Eigenpair<std::complex<double>, 13> eigendecomposition;
  std::array<std::complex<double>, 169> AT_values{};
  auto AT = yateto::DenseTensorView<2,std::complex<double>>(AT_values.data(), {13, 13});
  seissol::model::getTransposedCoefficientMatrix(*this, 0, AT);
  seissol::eigenvalues::computeEigenvaluesWithArmadillo(AT_values, eigendecomposition);
  double max_ev = std::numeric_limits<double>::lowest();
  for (int i = 0; i < 13; i++) {
    max_ev = eigendecomposition.values.at(i).real() > max_ev ? eigendecomposition.values.at(i).real() : max_ev;
  }
  return max_ev;
#else
  return std::sqrt(lambda + 2*mu) / rho;
#endif
}

double PoroElasticMaterial::getSWaveSpeed() const {
  return std::sqrt(mu / rho);
}

MaterialType PoroElasticMaterial::getMaterialType() const {
  return MaterialType::poroelastic;
}
