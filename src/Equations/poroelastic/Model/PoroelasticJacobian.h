#ifndef MODEL_POROELASTIC_JACOBIAN_H_
#define MODEL_POROELASTIC_JACOBIAN_H_

#include <armadillo>
#include <Eigen/Eigen>

namespace seissol {
  namespace model {
    struct PoroElasticMaterial;

    template<typename T>
    inline void setToZero(T& AT) {
      AT.setZero();
    };

    template<>
    inline void setToZero(arma::Mat<std::complex<double> >::fixed<13, 13>& AT) {
      for (size_t row = 0; row < AT.n_rows; row++) {
        for (size_t col = 0; col < AT.n_cols; col++) {
          AT(row, col) = 0;
        }
      }
    };

    template<typename T>
    inline void getTransposedCoefficientMatrix( PoroElasticMaterial const& material,
        unsigned dim,
        T& AT)
    {
      setToZero<T>(AT);

      Eigen::Matrix<double, 6, 1> alpha;
      alpha << 1 - (3*material.lambda + 2*material.mu) / (3*material.bulk_solid),
        1 - (3*material.lambda + 2*material.mu) / (3*material.bulk_solid),
        1 - (3*material.lambda + 2*material.mu) / (3*material.bulk_solid),
        -0.0,
        -0.0,
        -0.0;

      Eigen::Matrix<double, 6, 6> c;
      c << material.lambda + 2*material.mu , material.lambda , material.lambda , 0 , 0 , 0,
         material.lambda , material.lambda + 2*material.mu , material.lambda , 0 , 0 , 0,
         material.lambda , material.lambda , material.lambda + 2*material.mu , 0 , 0 , 0,
         0 , 0 , 0 , material.mu , 0 , 0,
         0 , 0 , 0 , 0 , material.mu , 0,
         0 , 0 , 0 , 0 , 0 , material.mu;

      double K_bar = material.lambda + 2*material.mu/3;
      double M = material.bulk_solid / (1 - material.porosity - K_bar / material.bulk_solid + material.porosity * material.bulk_solid / material.bulk_fluid);
      double m =  material.rho_fluid * material.tortuosity / material.porosity;

      Eigen::Matrix<double, 6, 6> c_bar = c + M * alpha * alpha.transpose();

      double rho_bar = (1 - material.porosity) * material.rho + material.porosity * material.rho_fluid;
      double rho_1 = rho_bar - material.rho_fluid*material.rho_fluid / m;
      double rho_2 = material.rho_fluid - m * rho_bar/material.rho_fluid;
      double beta_1 = material.rho_fluid / m;
      double beta_2 = rho_bar / material.rho_fluid;

      switch(dim){
        case 0:
          AT(0,6)  = -1 / rho_1;
          AT(0,10) = -1 / rho_2;
          AT(3,7)  = -1 / rho_1;
          AT(3,11) = -1 / rho_2;
          AT(5,8)  = -1 / rho_1;
          AT(5,12) = -1 / rho_2;

          AT(6, 0) = -c_bar(0, 0);
          AT(6, 1) = -c_bar(1, 0);
          AT(6, 2) = -c_bar(2, 0);
          AT(6, 3) = -c_bar(5, 0);
          AT(6, 4) = -c_bar(3, 0);
          AT(6, 5) = -c_bar(4, 0);
          AT(6, 9) = M * alpha(0);

          AT(7, 0) = -c_bar(0, 5);
          AT(7, 1) = -c_bar(1, 5);
          AT(7, 2) = -c_bar(2, 5);
          AT(7, 3) = -c_bar(5, 5);
          AT(7, 4) = -c_bar(3, 5);
          AT(7, 5) = -c_bar(4, 5);
          AT(7, 9) = M * alpha(5);

          AT(8, 0) = -c_bar(0, 4);
          AT(8, 1) = -c_bar(1, 4);
          AT(8, 2) = -c_bar(2, 4);
          AT(8, 3) = -c_bar(5, 4);
          AT(8, 4) = -c_bar(3, 4);
          AT(8, 5) = -c_bar(4, 4);
          AT(8, 9) = M * alpha(4);

          AT(9,6)  = - beta_1 / rho_1;
          AT(9,10) = - beta_2 / rho_2;

          AT(10,0) = - M*alpha(0);
          AT(10,1) = - M*alpha(1);
          AT(10,2) = - M*alpha(2);
          AT(10,3) = - M*alpha(5);
          AT(10,4) = - M*alpha(3);
          AT(10,5) = - M*alpha(4);
          AT(10,9) = M;
          break;
        case 1:
          AT(1,7)  = -1 / rho_1;
          AT(1,11) = -1 / rho_2;
          AT(3,6)  = -1 / rho_1;
          AT(3,10) = -1 / rho_2;
          AT(4,8)  = -1 / rho_1;
          AT(4,12) = -1 / rho_2;

          AT(6, 0) = -c_bar(0, 5);
          AT(6, 1) = -c_bar(1, 5);
          AT(6, 2) = -c_bar(2, 5);
          AT(6, 3) = -c_bar(5, 5);
          AT(6, 4) = -c_bar(3, 5);
          AT(6, 5) = -c_bar(4, 5);
          AT(6, 9) = M * alpha(5);

          AT(7, 0) = -c_bar(0, 1);
          AT(7, 1) = -c_bar(1, 1);
          AT(7, 2) = -c_bar(2, 1);
          AT(7, 3) = -c_bar(5, 1);
          AT(7, 4) = -c_bar(3, 1);
          AT(7, 5) = -c_bar(4, 1);
          AT(7, 9) = M * alpha(1);

          AT(8, 0) = -c_bar(0, 3);
          AT(8, 1) = -c_bar(1, 3);
          AT(8, 2) = -c_bar(2, 3);
          AT(8, 3) = -c_bar(5, 3);
          AT(8, 4) = -c_bar(3, 3);
          AT(8, 5) = -c_bar(4, 3);
          AT(8, 9) = M * alpha(3);

          AT(9,7)  = - beta_1 / rho_1;
          AT(9,11) = - beta_2 / rho_2;

          AT(11,0) = - M*alpha(0);
          AT(11,1) = - M*alpha(1);
          AT(11,2) = - M*alpha(2);
          AT(11,3) = - M*alpha(5);
          AT(11,4) = - M*alpha(3);
          AT(11,5) = - M*alpha(4);
          AT(11,9) = M;
          break;
        case 2:
          AT(2,8)  = -1 / rho_1;
          AT(2,12) = -1 / rho_2;
          AT(4,7)  = -1 / rho_1;
          AT(4,11) = -1 / rho_2;
          AT(5,6)  = -1 / rho_1;
          AT(5,10) = -1 / rho_2;

          AT(6, 0) = -c_bar(0, 4);
          AT(6, 1) = -c_bar(1, 4);
          AT(6, 2) = -c_bar(2, 4);
          AT(6, 3) = -c_bar(5, 4);
          AT(6, 4) = -c_bar(3, 4);
          AT(6, 5) = -c_bar(4, 4);
          AT(6, 9) = M * alpha(4);

          AT(7, 0) = -c_bar(0, 3);
          AT(7, 1) = -c_bar(1, 3);
          AT(7, 2) = -c_bar(2, 3);
          AT(7, 3) = -c_bar(5, 3);
          AT(7, 4) = -c_bar(3, 3);
          AT(7, 5) = -c_bar(4, 3);
          AT(7, 9) = M * alpha(3);

          AT(8, 0) = -c_bar(0, 2);
          AT(8, 1) = -c_bar(1, 2);
          AT(8, 2) = -c_bar(2, 2);
          AT(8, 3) = -c_bar(5, 2);
          AT(8, 4) = -c_bar(3, 2);
          AT(8, 5) = -c_bar(4, 2);
          AT(8, 9) = M * alpha(2);

          AT(9,8)  = - beta_1 / rho_1;
          AT(9,12) = - beta_2 / rho_2;

          AT(12,0) = - M*alpha(0);
          AT(12,1) = - M*alpha(1);
          AT(12,2) = - M*alpha(2);
          AT(12,3) = - M*alpha(5);
          AT(12,4) = - M*alpha(3);
          AT(12,5) = - M*alpha(4);
          AT(12,9) = M;
          break;

        default:
          break;
      }
    }
  }
}

#endif
