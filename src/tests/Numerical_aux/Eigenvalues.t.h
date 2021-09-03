#include <cxxtest/TestSuite.h>

#include <Numerical_aux/Eigenvalues.h>
#include <Kernels/precision.hpp>


namespace seissol {
  namespace unit_test {
    class EigenvaluesTestSuite;
  }
}

class seissol::unit_test::EigenvaluesTestSuite : public CxxTest::TestSuite
{
  private:
    static constexpr double epsilon = 10*std::numeric_limits<real>::epsilon();
  public:  
    void testEigenvalues3()
    {
      const unsigned dim = 3;

      std::array<std::array<std::complex<double>, dim*dim>, 2> matrices = {{
        { 2.0,  -3.0,  -3.0, -2.0,  1.0, -2.0, -1.0,  1.0,   4.0},
        {-2.0,   3.0,   3.0,  2.0, -1.0,  2.0,  1.0, -1.0,  -4.0},
      }};
      for (auto& M : matrices) {
        eigenvalues::Eigenpair<std::complex<double>, dim> eigenpair{};
        eigenvalues::computeEigenvaluesWithEigen3(M, eigenpair);
        testResidual<dim>(M, eigenpair);
#ifdef USE_POROELASTIC
        eigenvalues::computeEigenvaluesWithLapack(M, eigenpair);
        testResidual<dim>(M, eigenpair);
#endif
      }

    }

    template<unsigned dim>
    void testResidual(std::array<std::complex<double>,dim*dim>& M, eigenvalues::Eigenpair<std::complex<double>, dim>& eigenpair) {
      //compute M*R
      std::array<std::complex<double>, dim*dim> M_R{};
      for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
          for (unsigned int k = 0; k < dim; k++) {
            M_R[i+dim*j] += M[i+dim*k] * eigenpair.vectors[k+dim*j];
          }
        }
      }
      //compute R*L
      std::array<std::complex<double>, dim*dim> R_L{};
      for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
          R_L[i+dim*j] = eigenpair.vectors[i+dim*j] * eigenpair.values[j];
        }
      }
      //compare residual
      for (unsigned i = 0; i < dim*dim; i++) {
        TS_ASSERT_DELTA(std::abs(M_R[i] - R_L[i]), 0, epsilon);
      }
    }
};
