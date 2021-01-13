#include "datastructures.hpp"
#include "PoroelasticSetup.h"

using namespace seissol::model;
PoroElasticMaterial::PoroElasticMaterial( double* i_materialVal, int i_numMaterialVals)
{ 
#ifdef USE_POROELASTIC
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
#endif
}

void PoroElasticMaterial::getFullStiffnessTensor(std::array<real, 81>& fullTensor) const {
#ifdef USE_POROELASTIC
  double elasticMaterialVals[] = {this->rho, this->mu, this->lambda};
  ElasticMaterial em(elasticMaterialVals, 3);
  em.getFullStiffnessTensor(fullTensor);
#endif
}

double PoroElasticMaterial::getMaxWaveSpeed() const {
  return getPWaveSpeed();
}

double PoroElasticMaterial::getPWaveSpeed() const {
#ifdef USE_POROELASTIC
  Eigen::Matrix<double, 13, 13> AT = Eigen::Matrix<double, 13, 13>::Zero();
  seissol::model::getTransposedCoefficientMatrix(*this, 0, AT);
  Eigen::ComplexEigenSolver<Eigen::Matrix<double, 13, 13>> ces;
  ces.compute(AT);
  const auto eigenvalues = ces.eigenvalues();
  double max_ev = std::numeric_limits<double>::lowest();
  for (int i = 0; i < 13; i++) {
    max_ev = eigenvalues[i].real() > max_ev ? eigenvalues[i].real() : max_ev;
  }

  return max_ev;
#else
  return 0;
#endif
}

double PoroElasticMaterial::getSWaveSpeed() const {
#ifdef USE_POROELASTIC
  return std::sqrt(mu / rho);
#else
  return 0;
#endif
}

MaterialType PoroElasticMaterial::getMaterialType() const {
  return MaterialType::poroelastic;
}
