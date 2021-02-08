#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <iomanip>

#include <yateto/TensorView.h>

#include "Model/common.hpp"
#include "Kernels/common.hpp"
#include "Numerical_aux/Transformation.h"
#include "generated_code/init.h"

namespace seissol {
  namespace model {
#ifdef USE_POROELASTIC
    struct additionalPoroelasticParameters {
      Eigen::Matrix<double, 6, 1> alpha;
      double K_bar; 
      double M; 
      double m; 
      Eigen::Matrix<double, 6, 6> c_bar; 
      double rho_bar; 
      double rho_1; 
      double rho_2; 
      double beta_1; 
      double beta_2; 
    };

    inline additionalPoroelasticParameters getAdditionalParameters(const PoroElasticMaterial& material) {
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

      return {alpha, K_bar, M, m, c_bar, rho_bar, rho_1, rho_2, beta_1, beta_2};
    }

    template<typename T>
    inline void setToZero(T& AT) {
      AT.setZero();
    }

    template<>
    inline void setToZero(arma::Mat<std::complex<double> >::fixed<13, 13>& AT) {
      for (size_t row = 0; row < AT.n_rows; row++) {
        for (size_t col = 0; col < AT.n_cols; col++) {
          AT(row, col) = 0;
        }
      }
    }

    template<typename T>
    inline void getTransposedCoefficientMatrix( PoroElasticMaterial const& material,
        unsigned dim,
        T& AT)
    {
      setToZero<T>(AT);
      const additionalPoroelasticParameters params = getAdditionalParameters(material);
      switch(dim){
        case 0:
          AT(0,6)  = -1 / params.rho_1;
          AT(0,10) = -1 / params.rho_2;
          AT(3,7)  = -1 / params.rho_1;
          AT(3,11) = -1 / params.rho_2;
          AT(5,8)  = -1 / params.rho_1;
          AT(5,12) = -1 / params.rho_2;

          AT(6, 0) = -params.c_bar(0, 0);
          AT(6, 1) = -params.c_bar(1, 0); AT(6, 2) = -params.c_bar(2, 0);
          AT(6, 3) = -params.c_bar(5, 0);
          AT(6, 4) = -params.c_bar(3, 0);
          AT(6, 5) = -params.c_bar(4, 0);
          AT(6, 9) = params.M * params.alpha(0);

          AT(7, 0) = -params.c_bar(0, 5);
          AT(7, 1) = -params.c_bar(1, 5);
          AT(7, 2) = -params.c_bar(2, 5);
          AT(7, 3) = -params.c_bar(5, 5);
          AT(7, 4) = -params.c_bar(3, 5);
          AT(7, 5) = -params.c_bar(4, 5);
          AT(7, 9) = params.M * params.alpha(5);

          AT(8, 0) = -params.c_bar(0, 4);
          AT(8, 1) = -params.c_bar(1, 4);
          AT(8, 2) = -params.c_bar(2, 4);
          AT(8, 3) = -params.c_bar(5, 4);
          AT(8, 4) = -params.c_bar(3, 4);
          AT(8, 5) = -params.c_bar(4, 4);
          AT(8, 9) = params.M * params.alpha(4);

          AT(9,6)  = - params.beta_1 / params.rho_1;
          AT(9,10) = - params.beta_2 / params.rho_2;

          AT(10,0) = - params.M*params.alpha(0);
          AT(10,1) = - params.M*params.alpha(1);
          AT(10,2) = - params.M*params.alpha(2);
          AT(10,3) = - params.M*params.alpha(5);
          AT(10,4) = - params.M*params.alpha(3);
          AT(10,5) = - params.M*params.alpha(4);
          AT(10,9) = params.M;
          break;
        case 1:
          AT(1,7)  = -1 / params.rho_1;
          AT(1,11) = -1 / params.rho_2;
          AT(3,6)  = -1 / params.rho_1;
          AT(3,10) = -1 / params.rho_2;
          AT(4,8)  = -1 / params.rho_1;
          AT(4,12) = -1 / params.rho_2;

          AT(6, 0) = -params.c_bar(0, 5);
          AT(6, 1) = -params.c_bar(1, 5);
          AT(6, 2) = -params.c_bar(2, 5);
          AT(6, 3) = -params.c_bar(5, 5);
          AT(6, 4) = -params.c_bar(3, 5);
          AT(6, 5) = -params.c_bar(4, 5);
          AT(6, 9) = params.M * params.alpha(5);

          AT(7, 0) = -params.c_bar(0, 1);
          AT(7, 1) = -params.c_bar(1, 1);
          AT(7, 2) = -params.c_bar(2, 1);
          AT(7, 3) = -params.c_bar(5, 1);
          AT(7, 4) = -params.c_bar(3, 1);
          AT(7, 5) = -params.c_bar(4, 1);
          AT(7, 9) = params.M * params.alpha(1);

          AT(8, 0) = -params.c_bar(0, 3);
          AT(8, 1) = -params.c_bar(1, 3);
          AT(8, 2) = -params.c_bar(2, 3);
          AT(8, 3) = -params.c_bar(5, 3);
          AT(8, 4) = -params.c_bar(3, 3);
          AT(8, 5) = -params.c_bar(4, 3);
          AT(8, 9) = params.M * params.alpha(3);

          AT(9,7)  = - params.beta_1 / params.rho_1;
          AT(9,11) = - params.beta_2 / params.rho_2;

          AT(11,0) = - params.M*params.alpha(0);
          AT(11,1) = - params.M*params.alpha(1);
          AT(11,2) = - params.M*params.alpha(2);
          AT(11,3) = - params.M*params.alpha(5);
          AT(11,4) = - params.M*params.alpha(3);
          AT(11,5) = - params.M*params.alpha(4);
          AT(11,9) = params.M;
          break;
        case 2:
          AT(2,8)  = -1 / params.rho_1;
          AT(2,12) = -1 / params.rho_2;
          AT(4,7)  = -1 / params.rho_1;
          AT(4,11) = -1 / params.rho_2;
          AT(5,6)  = -1 / params.rho_1;
          AT(5,10) = -1 / params.rho_2;

          AT(6, 0) = -params.c_bar(0, 4);
          AT(6, 1) = -params.c_bar(1, 4);
          AT(6, 2) = -params.c_bar(2, 4);
          AT(6, 3) = -params.c_bar(5, 4);
          AT(6, 4) = -params.c_bar(3, 4);
          AT(6, 5) = -params.c_bar(4, 4);
          AT(6, 9) = params.M * params.alpha(4);

          AT(7, 0) = -params.c_bar(0, 3);
          AT(7, 1) = -params.c_bar(1, 3);
          AT(7, 2) = -params.c_bar(2, 3);
          AT(7, 3) = -params.c_bar(5, 3);
          AT(7, 4) = -params.c_bar(3, 3);
          AT(7, 5) = -params.c_bar(4, 3);
          AT(7, 9) = params.M * params.alpha(3);

          AT(8, 0) = -params.c_bar(0, 2);
          AT(8, 1) = -params.c_bar(1, 2);
          AT(8, 2) = -params.c_bar(2, 2);
          AT(8, 3) = -params.c_bar(5, 2);
          AT(8, 4) = -params.c_bar(3, 2);
          AT(8, 5) = -params.c_bar(4, 2);
          AT(8, 9) = params.M * params.alpha(2);

          AT(9,8)  = - params.beta_1 / params.rho_1;
          AT(9,12) = - params.beta_2 / params.rho_2;

          AT(12,0) = - params.M*params.alpha(0);
          AT(12,1) = - params.M*params.alpha(1);
          AT(12,2) = - params.M*params.alpha(2);
          AT(12,3) = - params.M*params.alpha(5);
          AT(12,4) = - params.M*params.alpha(3);
          AT(12,5) = - params.M*params.alpha(4);
          AT(12,9) = params.M;
          break;

        default:
          break;
      }
    }

    template<typename T>
    inline void getTransposedSourceCoefficientTensor( PoroElasticMaterial const& material, T& ET) {
      const additionalPoroelasticParameters params = getAdditionalParameters(material);
      const double e_1 = params.beta_1 * material.viscosity / (params.rho_1 * material.permeability);
      const double e_2 = params.beta_2 * material.viscosity / (params.rho_2 * material.permeability);

      ET.setZero();
      ET(10,6) = e_1;
      ET(11,7) = e_1;
      ET(12,8) = e_1;

      ET(10,10) = e_2;
      ET(11,11) = e_2;
      ET(12,12) = e_2;
    }

    template<>
    inline void getTransposedGodunovState( PoroElasticMaterial const&        local,
        PoroElasticMaterial const&        neighbor,
        FaceType                          faceType,
        init::QgodLocal::view::type&      QgodLocal,
        init::QgodNeighbor::view::type&   QgodNeighbor )
    {
      constexpr auto tolerance = 1.0e-10;

      using CMatrix = typename arma::Mat<std::complex<double>>::template fixed<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES>;
      using Matrix = typename arma::Mat<double>;
      using Vector = typename arma::Col<std::complex<double>>::template fixed<NUMBER_OF_QUANTITIES>;
      struct eigenDecomposition { Vector eigenvalues; CMatrix eigenvectors; };
      auto getEigenDecomposition = [&tolerance](PoroElasticMaterial const& material) {
        CMatrix coeff(arma::fill::zeros);
        getTransposedCoefficientMatrix(material, 0, coeff);

	CMatrix arma_eigenvectors(arma::fill::zeros);
	Vector arma_eigenvalues(arma::fill::zeros);
	arma::eig_gen(arma_eigenvalues, arma_eigenvectors, coeff.t());

#ifndef NDEBUG
        //check number of eigenvalues
        int ev_neg = 0;
        int ev_pos = 0;
        for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
          if (arma_eigenvalues(i).real() < -tolerance) {
            ++ev_neg;
          } else if (arma_eigenvalues(i).real() > tolerance) {
            ++ev_pos;
          }
        }
        assert(ev_neg == 4);
        assert(ev_pos == 4);

        //check whether eigensolver is good enough
        const CMatrix matrix_mult = coeff.t() * arma_eigenvectors;
	CMatrix eigenvalue_matrix(arma::fill::zeros);
        for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
          eigenvalue_matrix(i,i) = arma_eigenvalues(i);
        }
        const CMatrix vector_mult = arma_eigenvectors * eigenvalue_matrix;
        const CMatrix diff = matrix_mult - vector_mult;
        const double norm = arma::norm(diff);
        
        assert(norm < tolerance);
#endif
        return eigenDecomposition({arma_eigenvalues, arma_eigenvectors});
      };
      auto eigen_local = getEigenDecomposition(local);
      auto eigen_neighbor = getEigenDecomposition(neighbor); 	

      CMatrix chi_minus(arma::fill::zeros);
      CMatrix chi_plus(arma::fill::zeros);
      for(int i = 0; i < 13; i++) {
        if(eigen_local.eigenvalues(i).real() < -tolerance) {
          chi_minus(i,i) = 1.0;
        }
        if(eigen_local.eigenvalues(i).real() > tolerance) {
          chi_plus(i,i) = 1.0;
        }
      }

      //Matrix R = eigen_local.eigenvectors() * chi_minus + eigen_neighbor.eigenvectors() * chi_plus;
      CMatrix R = eigen_local.eigenvectors;
      if (faceType == FaceType::freeSurface) {
        Matrix R_real = arma::real(R);
        getTransposedFreeSurfaceGodunovState(false, QgodLocal, QgodNeighbor, R_real);
      } else {
	CMatrix R_inv = inv(R);
        CMatrix godunov_minus = R * chi_minus * R_inv;
        CMatrix godunov_plus =  R * chi_plus * R_inv;

        for (unsigned i = 0; i < QgodLocal.shape(1); ++i) {
          for (unsigned j = 0; j < QgodLocal.shape(0); ++j) {
            QgodLocal(j,i) = godunov_plus(i,j).real();
            QgodNeighbor(j,i) = godunov_minus(i,j).real();
#ifndef NDEBUG
	    assert(std::abs(godunov_plus(j,i).imag()) < tolerance);
	    assert(std::abs(godunov_minus(j,i).imag()) < tolerance);
#endif
          }
        }
      }
    }

    inline void calcZinv( yateto::DenseTensorView<2, real, unsigned> &Zinv, 
        yateto::DenseTensorView<2, real, unsigned> &sourceMatrix, 
        size_t quantity,
        real timeStepWidth) {
      using Matrix = Eigen::Matrix<real, CONVERGENCE_ORDER, CONVERGENCE_ORDER>;
      using Vector = Eigen::Matrix<real, CONVERGENCE_ORDER, 1>;

      Matrix Z(init::Z::Values);
      if(quantity >= 10) {
        for(int j = 0; j < CONVERGENCE_ORDER; j++) {
          Z(j,j) = Z(j,j) - timeStepWidth * sourceMatrix(quantity, quantity);
        }
      }

      auto solver = Z.colPivHouseholderQr();
      for(int col = 0; col < CONVERGENCE_ORDER; col++) {
        Vector rhs = Vector::Zero();
        rhs(col) = 1.0;
        auto Zinv_col = solver.solve(rhs);
        for(int row = 0; row < CONVERGENCE_ORDER; row++) {
          Zinv(row,col) = Zinv_col(row);
        }
      }
    }

    //constexpr for loop since we need to instatiate the view templates
     template<size_t i_start, size_t i_end>
     struct for_loop {
       for_loop(PoroelasticLocalData* ld,
           yateto::DenseTensorView<2, real, unsigned> &sourceMatrix, 
           real timeStepWidth) {
         auto Zinv = init::Zinv::view<i_start>::create(ld->Zinv[i_start]); 
         calcZinv(Zinv, sourceMatrix, i_start, timeStepWidth);
         if constexpr(i_start < i_end-1) {
           for_loop<i_start+1, i_end>(ld, sourceMatrix, timeStepWidth);
         }
       };
     };

    inline void initializeSpecificLocalData( PoroElasticMaterial const& material,
        real timeStepWidth,
        PoroelasticLocalData* localData )
    {
      auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
      sourceMatrix.setZero();
      getTransposedSourceCoefficientTensor(material, sourceMatrix);

      for_loop<0, NUMBER_OF_QUANTITIES>(localData, sourceMatrix, timeStepWidth);

      localData->typicalTimeStepWidth = timeStepWidth;
    }
#endif
  }
}
