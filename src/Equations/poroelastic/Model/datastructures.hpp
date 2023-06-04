#ifndef MODEL_POROELASTIC_DATASTRUCTURES_H_
#define MODEL_POROELASTIC_DATASTRUCTURES_H_

#include <cassert>
#include "Model/common_datastructures.hpp"
#include "Equations/elastic/Model/datastructures.hpp"

namespace seissol {
  namespace model {
    struct PoroElasticMaterial : ElasticMaterial {
      static constexpr std::size_t NumberOfQuantities = 13;
      static constexpr std::size_t NumberPerMechanism = 0;
      static constexpr std::size_t Mechanisms = 0;
      static constexpr MaterialType Type = MaterialType::poroelastic;

      double bulkSolid;
      // double lambda; // given by elasticity
      // double mu; // given by elasticity
      double porosity;
      double permeability;
      double tortuosity;
      double bulkFluid;
      double rhoFluid;
      double viscosity;

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
      virtual ~PoroElasticMaterial() {}

      void getFullStiffnessTensor(std::array<real, 81>& fullTensor) const 
      {
        double elasticMaterialVals[] = {this->rho, this->mu, this->lambda};
        ElasticMaterial em(elasticMaterialVals, 3);
        em.getFullStiffnessTensor(fullTensor);
      }

      double getMaxWaveSpeed() const 
      {
        return getPWaveSpeed();
      }

      //only declare it here and define in a separate datastructures.cpp
      //to circumvent problems with circular includes
      double getPWaveSpeed() const;

      double getSWaveSpeed() const
      {
        return std::sqrt(mu / rho);
      }

      MaterialType getMaterialType() const {
        return Type;
      }
    };
  }
}


#endif
