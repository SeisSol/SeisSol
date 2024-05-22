#ifndef MODEL_POROELASTIC_DATASTRUCTURES_H_
#define MODEL_POROELASTIC_DATASTRUCTURES_H_

#include <cassert>
#include "Model/common_datastructures.hpp"
#include "Equations/elastic/Model/datastructures.hpp"

namespace seissol {
  namespace model {
    struct PoroElasticMaterial : Material {
      double bulkSolid;
      double lambda;
      double mu;
      double porosity;
      double permeability;
      double tortuosity;
      double bulkFluid;
      double rhoFluid;
      double viscosity;

      double getLambda() const override {
        return lambda;
      }

      double getMu() const override {
        return mu;
      }

      PoroElasticMaterial() {}

      PoroElasticMaterial( double* materialValues, int numMaterialValues)
      { 
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

      //only declare it here and define in a separate datastructures.cpp
      //to circumvent problems with circular includes
      double getPWaveSpeed() const final;

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
