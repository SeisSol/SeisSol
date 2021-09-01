#ifndef MODEL_POROELASTICSETUP_H_
#define MODEL_POROELASTICSETUP_H_

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

#include <cassert>

#include <yateto/TensorView.h>

#include "Model/common.hpp"
#include "Kernels/common.hpp"
#include "Numerical_aux/Transformation.h"
#include "Numerical_aux/Eigenvalues.h"
#include "generated_code/init.h"

namespace seissol {
  namespace model {
    struct AdditionalPoroelasticParameters {
      Eigen::Matrix<double, 6, 1> alpha;
      double KBar; 
      double M; 
      double m; 
      Eigen::Matrix<double, 6, 6> cBar; 
      double rhoBar; 
      double rho1; 
      double rho2; 
      double beta1; 
      double beta2; 
    };

    inline AdditionalPoroelasticParameters getAdditionalParameters(const PoroElasticMaterial& material) {
      Eigen::Matrix<double, 6, 1> alpha;
      alpha << 1 - (3*material.lambda + 2*material.mu) / (3*material.bulkSolid),
        1 - (3*material.lambda + 2*material.mu) / (3*material.bulkSolid),
        1 - (3*material.lambda + 2*material.mu) / (3*material.bulkSolid),
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

      double KBar = material.lambda + 2*material.mu/3;
      double M = material.bulkSolid / (1 - material.porosity - KBar / material.bulkSolid + material.porosity * material.bulkSolid / material.bulkFluid);
      double m =  material.rhoFluid * material.tortuosity / material.porosity;

      Eigen::Matrix<double, 6, 6> cBar = c + M * alpha * alpha.transpose();

      double rhoBar = (1 - material.porosity) * material.rho + material.porosity * material.rhoFluid;
      double rho1 = rhoBar - material.rhoFluid*material.rhoFluid / m;
      double rho2 = material.rhoFluid - m * rhoBar/material.rhoFluid;
      double beta1 = material.rhoFluid / m;
      double beta2 = rhoBar / material.rhoFluid;

      return {alpha, KBar, M, m, cBar, rhoBar, rho1, rho2, beta1, beta2};
    }

    template<typename T>
    inline void setToZero(T& AT) {
      AT.setZero();
    }

    template<>
    inline void setToZero(arma::Mat<std::complex<double>>::fixed<13,13>& AT) {
      AT.zeros();
    }

    template<typename T>
    inline void getTransposedCoefficientMatrix( PoroElasticMaterial const& material,
                                                unsigned dim,
                                                T& AT)
    {
      setToZero<T>(AT);
      const AdditionalPoroelasticParameters params = getAdditionalParameters(material);
      switch(dim){
        case 0:
          AT(0,6)  = -1 / params.rho1;
          AT(0,10) = -1 / params.rho2;
          AT(3,7)  = -1 / params.rho1;
          AT(3,11) = -1 / params.rho2;
          AT(5,8)  = -1 / params.rho1;
          AT(5,12) = -1 / params.rho2;

          AT(6, 0) = -params.cBar(0, 0);
          AT(6, 1) = -params.cBar(1, 0);
          AT(6, 2) = -params.cBar(2, 0);
          AT(6, 3) = -params.cBar(5, 0);
          AT(6, 4) = -params.cBar(3, 0);
          AT(6, 5) = -params.cBar(4, 0);
          AT(6, 9) = params.M * params.alpha(0);

          AT(7, 0) = -params.cBar(0, 5);
          AT(7, 1) = -params.cBar(1, 5);
          AT(7, 2) = -params.cBar(2, 5);
          AT(7, 3) = -params.cBar(5, 5);
          AT(7, 4) = -params.cBar(3, 5);
          AT(7, 5) = -params.cBar(4, 5);
          AT(7, 9) = params.M * params.alpha(5);

          AT(8, 0) = -params.cBar(0, 4);
          AT(8, 1) = -params.cBar(1, 4);
          AT(8, 2) = -params.cBar(2, 4);
          AT(8, 3) = -params.cBar(5, 4);
          AT(8, 4) = -params.cBar(3, 4);
          AT(8, 5) = -params.cBar(4, 4);
          AT(8, 9) = params.M * params.alpha(4);

          AT(9,6)  = - params.beta1 / params.rho1;
          AT(9,10) = - params.beta2 / params.rho2;

          AT(10,0) = - params.M*params.alpha(0);
          AT(10,1) = - params.M*params.alpha(1);
          AT(10,2) = - params.M*params.alpha(2);
          AT(10,3) = - params.M*params.alpha(5);
          AT(10,4) = - params.M*params.alpha(3);
          AT(10,5) = - params.M*params.alpha(4);
          AT(10,9) = params.M;
          break;
        case 1:
          AT(1,7)  = -1 / params.rho1;
          AT(1,11) = -1 / params.rho2;
          AT(3,6)  = -1 / params.rho1;
          AT(3,10) = -1 / params.rho2;
          AT(4,8)  = -1 / params.rho1;
          AT(4,12) = -1 / params.rho2;

          AT(6, 0) = -params.cBar(0, 5);
          AT(6, 1) = -params.cBar(1, 5);
          AT(6, 2) = -params.cBar(2, 5);
          AT(6, 3) = -params.cBar(5, 5);
          AT(6, 4) = -params.cBar(3, 5);
          AT(6, 5) = -params.cBar(4, 5);
          AT(6, 9) = params.M * params.alpha(5);

          AT(7, 0) = -params.cBar(0, 1);
          AT(7, 1) = -params.cBar(1, 1);
          AT(7, 2) = -params.cBar(2, 1);
          AT(7, 3) = -params.cBar(5, 1);
          AT(7, 4) = -params.cBar(3, 1);
          AT(7, 5) = -params.cBar(4, 1);
          AT(7, 9) = params.M * params.alpha(1);

          AT(8, 0) = -params.cBar(0, 3);
          AT(8, 1) = -params.cBar(1, 3);
          AT(8, 2) = -params.cBar(2, 3);
          AT(8, 3) = -params.cBar(5, 3);
          AT(8, 4) = -params.cBar(3, 3);
          AT(8, 5) = -params.cBar(4, 3);
          AT(8, 9) = params.M * params.alpha(3);

          AT(9,7)  = - params.beta1 / params.rho1;
          AT(9,11) = - params.beta2 / params.rho2;

          AT(11,0) = - params.M*params.alpha(0);
          AT(11,1) = - params.M*params.alpha(1);
          AT(11,2) = - params.M*params.alpha(2);
          AT(11,3) = - params.M*params.alpha(5);
          AT(11,4) = - params.M*params.alpha(3);
          AT(11,5) = - params.M*params.alpha(4);
          AT(11,9) = params.M;
          break;
        case 2:
          AT(2,8)  = -1 / params.rho1;
          AT(2,12) = -1 / params.rho2;
          AT(4,7)  = -1 / params.rho1;
          AT(4,11) = -1 / params.rho2;
          AT(5,6)  = -1 / params.rho1;
          AT(5,10) = -1 / params.rho2;

          AT(6, 0) = -params.cBar(0, 4);
          AT(6, 1) = -params.cBar(1, 4);
          AT(6, 2) = -params.cBar(2, 4);
          AT(6, 3) = -params.cBar(5, 4);
          AT(6, 4) = -params.cBar(3, 4);
          AT(6, 5) = -params.cBar(4, 4);
          AT(6, 9) = params.M * params.alpha(4);

          AT(7, 0) = -params.cBar(0, 3);
          AT(7, 1) = -params.cBar(1, 3);
          AT(7, 2) = -params.cBar(2, 3);
          AT(7, 3) = -params.cBar(5, 3);
          AT(7, 4) = -params.cBar(3, 3);
          AT(7, 5) = -params.cBar(4, 3);
          AT(7, 9) = params.M * params.alpha(3);

          AT(8, 0) = -params.cBar(0, 2);
          AT(8, 1) = -params.cBar(1, 2);
          AT(8, 2) = -params.cBar(2, 2);
          AT(8, 3) = -params.cBar(5, 2);
          AT(8, 4) = -params.cBar(3, 2);
          AT(8, 5) = -params.cBar(4, 2);
          AT(8, 9) = params.M * params.alpha(2);

          AT(9,8)  = - params.beta1 / params.rho1;
          AT(9,12) = - params.beta2 / params.rho2;

          AT(12,0) = - params.M*params.alpha(0);
          AT(12,1) = - params.M*params.alpha(1);
          AT(12,2) = - params.M*params.alpha(2);
          AT(12,3) = - params.M*params.alpha(5);
          AT(12,4) = - params.M*params.alpha(3);
          AT(12,5) = - params.M*params.alpha(4);
          AT(12,9) = params.M;
          break;

        default:
          logError() << "Cannot create transposed coefficient matrix for dimension " << dim << ", has to be either 0, 1 or 2.";
      }
    }

    template<typename T>
    inline void getTransposedSourceCoefficientTensor( PoroElasticMaterial const& material, T& ET) {
      const AdditionalPoroelasticParameters params = getAdditionalParameters(material);
      const double e1 = params.beta1 * material.viscosity / (params.rho1 * material.permeability);
      const double e2 = params.beta2 * material.viscosity / (params.rho2 * material.permeability);

      ET.setZero();
      ET(10,6) = e1;
      ET(11,7) = e1;
      ET(12,8) = e1;

      ET(10,10) = e2;
      ET(11,11) = e2;
      ET(12,12) = e2;
    }

    template<typename T>
    void getTransposedFreeSurfaceGodunovState( MaterialType materialtype,
                                               T&        QgodLocal,
                                               T&        QgodNeighbor,
                                               arma::Mat<double>::fixed<13,13>& R)
    {
      if (materialtype != MaterialType::poroelastic) {
        logError() << "This specialization is for armadillo matrices. This is only used for poroelastic materials. You should never end up here.";
      }
    
      constexpr size_t relevantQuantities = NUMBER_OF_QUANTITIES - 6*NUMBER_OF_RELAXATION_MECHANISMS;
      for (size_t i = 0; i < relevantQuantities; i++) {
        for (size_t j = 0; j < relevantQuantities; j++) {
          QgodNeighbor(i,j) = std::numeric_limits<double>::signaling_NaN();
        }
      }
    
      QgodLocal.setZero();
      arma::uvec tractionIndices = {0,3,5,9};
      arma::uvec velocityIndices = {6,7,8,10,11,12};
      arma::uvec columnIndices = {5, 7, 9, 11};
      arma::mat R11 = R.submat(tractionIndices, columnIndices);
      arma::mat R21 = R.submat(velocityIndices, columnIndices);
      arma::mat S = (-(R21 * inv(R11))).eval();
      setBlocks(QgodLocal, S, tractionIndices, velocityIndices);
    }

    template<>
    inline void getTransposedGodunovState( PoroElasticMaterial const&        local,
                                           PoroElasticMaterial const&        neighbor,
                                           FaceType                          faceType,
                                           init::QgodLocal::view::type&      QgodLocal,
                                           init::QgodNeighbor::view::type&   QgodNeighbor )
    {
      //Will be used to check, whether numbers are (numerically) zero
      constexpr auto zeroThreshold = 1e-7;
      using CMatrix = typename arma::Mat<std::complex<double>>::template fixed<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES>;
      using Matrix = typename arma::Mat<double>::template fixed<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES>;
      using CVector = typename arma::Col<std::complex<double>>::template fixed<NUMBER_OF_QUANTITIES>;
      auto getEigenDecomposition = [&zeroThreshold](PoroElasticMaterial const& material) {
        CMatrix coeff(arma::fill::zeros);
        getTransposedCoefficientMatrix(material, 0, coeff);

        CMatrix armaEigenvectors(arma::fill::zeros);
        CVector armaEigenvalues(arma::fill::zeros);
        arma::eig_gen(armaEigenvalues, armaEigenvectors, coeff.t(), "balance");

#ifndef NDEBUG
        //check number of eigenvalues
        //also check that the imaginary parts are zero
        int evNeg = 0;
        int evPos = 0;
        for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
          assert(std::abs(armaEigenvalues(i).imag()) < zeroThreshold);
          if (armaEigenvalues(i).real() < -zeroThreshold) {
            ++evNeg;
          } else if (armaEigenvalues(i).real() > zeroThreshold) {
            ++evPos;
          }
        }
        assert(evNeg == 4);
        assert(evPos == 4);

        //check whether eigensolver is good enough
        const CMatrix matrixMult = coeff.t() * armaEigenvectors;
        CMatrix eigenvalueMatrix(arma::fill::zeros);
        for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
          eigenvalueMatrix(i,i) = armaEigenvalues(i);
        }
        const CMatrix vectorMult = armaEigenvectors * eigenvalueMatrix;
        const CMatrix diff = matrixMult - vectorMult;
        const double norm = arma::norm(diff);
        
        std::stringstream messageStream;
        messageStream << "Residual " << norm << " is larger than " << zeroThreshold << ": Eigensolver is not accurate enough";
        assert((messageStream.str().c_str(), norm < zeroThreshold));
#endif
        return std::pair<CVector, CMatrix>({armaEigenvalues, armaEigenvectors});
      };
      auto [localEigenvalues, localEigenvectors] = getEigenDecomposition(local);
      auto [neighborEigenvalues, neighborEigenvectors] = getEigenDecomposition(neighbor); 


      CMatrix chiMinus(arma::fill::zeros);
      CMatrix chiPlus(arma::fill::zeros);
      for(int i = 0; i < 13; i++) {
        if(localEigenvalues(i).real() < -zeroThreshold) {
          chiMinus(i,i) = 1.0;
        }
        if(localEigenvalues(i).real() > zeroThreshold) {
          chiPlus(i,i) = 1.0;
        }
      }

      CMatrix R = localEigenvectors * chiMinus + neighborEigenvectors * chiPlus;
      //set null space eigenvectors manually
      R(1,0) = 1.0;
      R(2,1) = 1.0;
      R(12,2) = 1.0;
      R(11,3) = 1.0;
      R(4,12) = 1.0;
      if (faceType == FaceType::freeSurface) {
        Matrix realR = arma::real(R);
        getTransposedFreeSurfaceGodunovState(MaterialType::poroelastic, QgodLocal, QgodNeighbor, realR);
      } else {
        CMatrix invR = inv(R);
        CMatrix godunovMinus = R * chiMinus * invR;
        CMatrix godunovPlus =  R * chiPlus * invR;

        for (unsigned i = 0; i < QgodLocal.shape(1); ++i) {
          for (unsigned j = 0; j < QgodLocal.shape(0); ++j) {
            QgodLocal(j,i) = godunovPlus(i,j).real();
            QgodNeighbor(j,i) = godunovMinus(i,j).real();
            assert(std::abs(godunovPlus(j,i).imag()) < zeroThreshold);
            assert(std::abs(godunovMinus(j,i).imag()) < zeroThreshold);
          }
        }
      }
    }

    template<typename Tview>
    inline void calcZinv( yateto::DenseTensorView<2, real, unsigned> &Zinv, 
        Tview &sourceMatrix, 
        size_t quantity,
        real timeStepWidth) {
      using Matrix = Eigen::Matrix<real, CONVERGENCE_ORDER, CONVERGENCE_ORDER>;
      using Vector = Eigen::Matrix<real, CONVERGENCE_ORDER, 1>;

      Matrix Z(init::Z::Values);
      //sourceMatrix[i,i] = 0 for i < 10
      //This is specific to poroelasticity, so change this for another equation
      //We need this check, because otherwise the lookup sourceMatrix(quantity, quantity) fails
      if(quantity >= 10) {
        Z -= timeStepWidth * sourceMatrix(quantity, quantity) * Matrix::Identity();
      }

      auto solver = Z.colPivHouseholderQr();
      for(int col = 0; col < CONVERGENCE_ORDER; col++) {
        Vector rhs = Vector::Zero();
        rhs(col) = 1.0;
        auto ZinvCol = solver.solve(rhs);
        for(int row = 0; row < CONVERGENCE_ORDER; row++) {
          //save as transposed
          Zinv(col,row) = ZinvCol(row);
        }
      }
    }

    //constexpr for loop since we need to instatiate the view templates
    template<size_t iStart, size_t iEnd, typename Tview>
    struct zInvInitializerForLoop {
      zInvInitializerForLoop(real ZinvData[NUMBER_OF_QUANTITIES][CONVERGENCE_ORDER*CONVERGENCE_ORDER],
          Tview &sourceMatrix, 
          real timeStepWidth) {
        auto Zinv = init::Zinv::view<iStart>::create(ZinvData[iStart]); 
        calcZinv(Zinv, sourceMatrix, iStart, timeStepWidth);
        if constexpr(iStart < iEnd-1) {
          zInvInitializerForLoop<iStart+1, iEnd, Tview>(ZinvData, sourceMatrix, timeStepWidth);
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

      zInvInitializerForLoop<0, NUMBER_OF_QUANTITIES, decltype(sourceMatrix)>(localData->Zinv, sourceMatrix, timeStepWidth);
      std::fill(localData->G, localData->G+NUMBER_OF_QUANTITIES, 0.0);
      localData->G[10] = sourceMatrix(10, 6);
      localData->G[11] = sourceMatrix(11, 7);
      localData->G[12] = sourceMatrix(12, 8);

      localData->typicalTimeStepWidth = timeStepWidth;
    }
  }
}
#endif
