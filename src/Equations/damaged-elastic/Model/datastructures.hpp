#ifndef MODEL_DAMAGED_DATASTRUCTURES_H_
#define MODEL_DAMAGED_DATASTRUCTURES_H_

#include <Model/common_datastructures.hpp>
#include <cmath>
#include <generated_code/kernel.h>
#include <generated_code/init.h>

namespace seissol {
  namespace model {
    struct DamagedElasticMaterial : Material {
      double lambda;
      double mu;
      double sigmaxx_alpha;
      double sigmaxy_alpha;
      double sigmaxz_alpha;
      double sigmayx_alpha;
      double sigmayy_alpha;
      double sigmayz_alpha;
      double sigmazx_alpha;
      double sigmazy_alpha;
      double sigmazz_alpha;

      DamagedElasticMaterial() {};
      DamagedElasticMaterial(double* materialValues, int numMaterialValues)
      {
        assert(numMaterialValues == 12);

        this->rho = materialValues[0];
        this->mu = materialValues[1];
        this->lambda = materialValues[2];
        this->sigmaxx_alpha = materialValues[3];
        this->sigmaxy_alpha = materialValues[4];
        this->sigmaxz_alpha = materialValues[5];
        this->sigmayx_alpha = materialValues[6];
        this->sigmayy_alpha = materialValues[7];
        this->sigmayz_alpha = materialValues[8];
        this->sigmazx_alpha = materialValues[9];
        this->sigmazy_alpha = materialValues[10];
        this->sigmazz_alpha = materialValues[11];
      }

      virtual ~DamagedElasticMaterial() {};

      void getFullStiffnessTensor(std::array<real, 81>& fullTensor) const final {

        auto stiffnessTensorView = init::stiffnessTensor::view::create(fullTensor.data());
        stiffnessTensorView.setZero();
        stiffnessTensorView(0,0,0,0) = lambda + 2*mu;
        stiffnessTensorView(0,0,1,1) = lambda;
        stiffnessTensorView(0,0,2,2) = lambda;
        stiffnessTensorView(0,1,0,1) = mu;
        stiffnessTensorView(0,1,1,0) = mu;
        stiffnessTensorView(0,2,0,2) = mu;
        stiffnessTensorView(0,2,2,0) = mu;
        stiffnessTensorView(1,0,0,1) = mu;
        stiffnessTensorView(1,0,1,0) = mu;
        stiffnessTensorView(1,1,0,0) = lambda;
        stiffnessTensorView(1,1,1,1) = lambda + 2*mu;
        stiffnessTensorView(1,1,2,2) = lambda;
        stiffnessTensorView(1,2,1,2) = mu;
        stiffnessTensorView(1,2,2,1) = mu;
        stiffnessTensorView(2,0,0,2) = mu;
        stiffnessTensorView(2,0,2,0) = mu;
        stiffnessTensorView(2,1,2,1) = mu;
        stiffnessTensorView(2,2,0,0) = lambda;
        stiffnessTensorView(2,2,1,1) = lambda;
        stiffnessTensorView(2,2,2,2) = lambda + 2*mu;
      }

      double getMaxWaveSpeed() const final {
        return getPWaveSpeed();
      }

      double getPWaveSpeed() const final {
        return std::sqrt((lambda + 2*mu) / rho);
      }

      double getSWaveSpeed() const final {
        return std::sqrt(mu / rho);
      }

      MaterialType getMaterialType() const override{
        return MaterialType::damaged;
      }
    };
  }
}

#endif
