#ifndef MODEL_POROELASTIC_DATASTRUCTURES_H_
#define MODEL_POROELASTIC_DATASTRUCTURES_H_

#include "Model/common_datastructures.hpp"

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
      PoroElasticMaterial( double* i_materialVal, int i_numMaterialVals);
      virtual ~PoroElasticMaterial() {};

      void getFullStiffnessTensor(std::array<real, 81>& fullTensor) const final;
      double getMaxWaveSpeed() const final;
      double getPWaveSpeed() const final;
      double getSWaveSpeed() const final;
      MaterialType getMaterialType() const override;
    };
  }
}

#endif
