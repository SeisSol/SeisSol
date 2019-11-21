#include <cmath>
#include <array>
#include <stdexcept>
#include <numeric>
#ifdef HAS_EIGEN3
#include <Eigen/Eigenvalues>
#endif

#include <Kernels/precision.hpp>
#include <Physics/InitialField.h>
#include <Model/Setup.h>
#include <Solver/Interoperability.h>

extern seissol::Interoperability e_interoperability;

seissol::physics::Planarwave::Planarwave(double phase)
  : m_setVar(2),
    m_kVec{3.14159265358979323846, 3.14159265358979323846, 3.14159265358979323846},
    m_phase(phase)
{
#ifdef HAS_EIGEN3
  const double rho = 1.0;
  const double mu = 1.0;
  const double lambda = 2.0;
  const double Qp = 20.0;
  const double Qs = 10.0;


  seissol::model::ViscoElasticMaterial material;
  e_interoperability.fitAttenuation(rho, mu, lambda, Qp, Qs, material);

  std::complex<double> planeWaveOperator[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
  model::getPlaneWaveOperator(material, m_kVec.data(), planeWaveOperator);

  using Matrix = Eigen::Matrix<std::complex<double>, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES, Eigen::ColMajor>;
  Matrix A(planeWaveOperator);
  Eigen::ComplexEigenSolver<Matrix> ces;
  ces.compute(A);

  auto eigenvalues = ces.eigenvalues();
  for (size_t i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    m_lambdaA[i] = eigenvalues(i,0);
  }

  std::vector<size_t> varField(NUMBER_OF_QUANTITIES);
  std::iota(varField.begin(), varField.end(), 0);

  std::sort(varField.begin(), varField.end(), [&eigenvalues](size_t a, size_t b) {
    return eigenvalues[a].real() < eigenvalues[b].real();
  });

  // Select S-wave in opposite direction (1) and P-wave along direction (last)
  std::array<size_t, 2> selectVars = {1, NUMBER_OF_QUANTITIES-1};
  assert(m_setVar == selectVars.size());

  for (auto& var : selectVars) {
    m_varField.push_back(var);
    m_ampField.push_back(1.0);
  }

  auto& eigenvectors = ces.eigenvectors();
  auto R = yateto::DenseTensorView<2,std::complex<double>>(m_eigenvectors, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});
  for (size_t j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
    for (size_t i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      R(i,j) = eigenvectors(i,j);
    }
  }
#else
  throw std::runtime_error("Eigen3 required for anelastic planarwave.");
#endif
}

void seissol::physics::Planarwave::evaluate(  double time,
                                              std::vector<std::array<double, 3>> const& points,
                                              yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const
{
  dofsQP.setZero();

  auto R = yateto::DenseTensorView<2,std::complex<double>>(
             const_cast<std::complex<double>*>(m_eigenvectors),
             {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES}
           );
  for (int v = 0; v < m_setVar; ++v) {
    const auto omega =  m_lambdaA[m_varField[v]];
    for (unsigned j = 0; j < dofsQP.shape(1); ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        dofsQP(i,j) += (R(j,m_varField[v]) * m_ampField[v] *
                        std::exp(std::complex<double>(0.0, 1.0) * (
                          omega * time - m_kVec[0]*points[i][0] - m_kVec[1]*points[i][1] - m_kVec[2]*points[i][2] + m_phase
                        ))).real();
      }
    }
  }
}
