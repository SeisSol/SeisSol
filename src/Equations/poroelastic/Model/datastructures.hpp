#ifndef MODEL_POROELASTIC_DATASTRUCTURES_H_
#define MODEL_POROELASTIC_DATASTRUCTURES_H_

#include <cassert>
#include "Model/common_datastructures.hpp"
#include "Equations/elastic/Model/datastructures.hpp"

namespace seissol {
  namespace model {
    struct PoroElasticMaterial : Material {
      double bulk_solid;
      double lambda;
      double mu;
      double porosity;
      double permeability;
      double tortuosity;
      double bulk_fluid;
      double rho_fluid;
      double viscosity;

      PoroElasticMaterial() {}

      PoroElasticMaterial( double* i_materialVal, int i_numMaterialVals)
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
      };
      virtual ~PoroElasticMaterial() {};

      void getFullStiffnessTensor(std::array<real, 81>& fullTensor) const final 
      {
        double elasticMaterialVals[] = {this->rho, this->mu, this->lambda};
        ElasticMaterial em(elasticMaterialVals, 3);
        em.getFullStiffnessTensor(fullTensor);
      };

      double getMaxWaveSpeed() const final 
      {
        return getPWaveSpeed();
      }

#ifndef USE_POROELASTIC
      //Return an estimate
      double getPWaveSpeed() const final 
      {
        return std::sqrt(lambda + 2*mu) / rho;
      };
#else
      //only declare it here and define in a separate datastructures.cpp
      //to circumvent problems with circular includes
      double getPWaveSpeed() const final;
#endif

      double getSWaveSpeed() const final
      {
        return std::sqrt(mu / rho);
      };

      MaterialType getMaterialType() const override {
        return MaterialType::poroelastic;
      };
    };
  }
}


#endif
