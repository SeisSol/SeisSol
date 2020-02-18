#include <cmath>
#include <array>
#include <stdexcept>
#include <Eigen/Eigenvalues>

#include <Kernels/precision.hpp>
#include <Physics/InitialField.h>
#include <Model/Setup.h>
#include <Solver/Interoperability.h>
#include <Initializer/ParameterDB.h>

extern seissol::Interoperability e_interoperability;

seissol::physics::PoroelasticPlanarwave::PoroelasticPlanarwave(seissol::model::Material material, double phase)
  : m_setVar(3),
    m_kVec{3.14159265358979323846, 3.14159265358979323846, 3.14159265358979323846},
    m_phase(phase)
{

  std::cout << "material values: " << std::endl
  << "bulk_solid   = " << material.bulk_solid  << std::endl
  << "rho_solid    = " << material.rho  << std::endl
  << "lambda       = " << material.lambda  << std::endl
  << "mu           = " << material.mu  << std::endl
  << "porosity     = " << material.porosity  << std::endl
  << "permeability = " << material.permeability  << std::endl
  << "tortuosity   = " << material.tortuosity  << std::endl
  << "bulk_fluid   = " << material.bulk_fluid  << std::endl
  << "rho_fluid    = " << material.rho_fluid  << std::endl
  << "viscosity    = " << material.viscosity << std::endl;

  std::complex<real> planeWaveOperator[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
  model::getPlaneWaveOperator(material, m_kVec.data(), planeWaveOperator);

  using Matrix = Eigen::Matrix<std::complex<real>, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES, Eigen::ColMajor>;
  Matrix A(planeWaveOperator);
  Eigen::ComplexEigenSolver<Matrix> ces;
  ces.compute(A);

  auto eigenvalues = ces.eigenvalues();
  std::cout << eigenvalues << std::endl;
  for (size_t i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    m_lambdaA[i] = eigenvalues(i,0);
  }

  auto& eigenvectors = ces.eigenvectors();
  std::vector<size_t> varField(NUMBER_OF_QUANTITIES);
  std::iota(varField.begin(), varField.end(), 0);

  std::sort(varField.begin(), varField.end(), [&eigenvalues](size_t a, size_t b) {
      return eigenvalues[a].real() < eigenvalues[b].real();
      });

  // Select S-wave and slow P-wave in opposite direction (1, 3) and P-wave along direction (last)
  std::array<size_t, 3> selectVars = {1, 3, NUMBER_OF_QUANTITIES-1};
  assert(m_setVar == selectVars.size());
  for (int i = 0; i < m_setVar; i++) {
    m_varField.push_back(selectVars[i]);
    m_ampField.push_back(100.0);
  }

//  using Vector = Eigen::Matrix<std::complex<double>, NUMBER_OF_QUANTITIES, 1, Eigen::ColMajor>;
//  auto solver = eigenvectors.colPivHouseholderQr();
//  Vector rhs = Vector::Zero();
//  rhs(1) = 1.0;
//  auto ampField = solver.solve(rhs);
//  std::cout << ampField << std::endl;
//  for (int i = 0; i < m_setVar; i++) {
//    m_varField.push_back(i);
//    m_ampField.push_back(ampField(i));
//  }

  auto R = yateto::DenseTensorView<2,std::complex<double>>(m_eigenvectors, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});
  for (size_t j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
    for (size_t i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      R(i,j) = eigenvectors(i,j);
    }
  }
}

void seissol::physics::PoroelasticPlanarwave::evaluate(double time,
                                                       std::vector<std::array<double, 3>> const& points,
                                                       const CellMaterialData& materialData,
                                                       yateto::DenseTensorView<2,double,unsigned>& dofsQP) const
{
  dofsQP.setZero();

  auto R = yateto::DenseTensorView<2,std::complex<real>>(
             const_cast<std::complex<real>*>(m_eigenvectors),
             {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES}
           );
  for (int v = 0; v < m_setVar; ++v) {
    const auto omega =  m_lambdaA[m_varField[v]];
    for (unsigned j = 0; j < dofsQP.shape(1); ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        dofsQP(i,j) += (R(j,m_varField[v]) * m_ampField[v] *
                        std::exp(std::complex<real>(0.0, 1.0) * (
                          omega * time 
                          - m_kVec[0]*points[i][0] 
                          - m_kVec[1]*points[i][1] 
                          - m_kVec[2]*points[i][2] 
                          + m_phase
                          ))).real();
      }
    }
  }
}
