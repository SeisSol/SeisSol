#include <cmath>
#include <array>
#include <Eigen/Eigenvalues>

#include <Kernels/precision.hpp>
#include <Physics/InitialField.h>
#include <Model/Setup.h>
#include <yateto/TensorView.h>

seissol::physics::Planarwave::Planarwave(real phase)
  : m_setVar(9),
    m_kVec{6.283185307179587E-002, 6.283185307179587E-002, 6.283185307179587E-002},
    m_phase(phase)
{ 

  double materialVal[22] = { 
      1.0, //rho
    192.0, //c11
     66.0, //c12
     60.0, //c13
      0.0, //c14
      0.0, //c15
      0.0, //c16
    160.0, //c22
     56.0, //c23
      0.0, //c24
      0.0, //c25
      0.0, //c26
    272.0, //c33
      0.0, //c34
      0.0, //c35
      0.0, //c36
     60.0, //c44
      0.0, //c45
      0.0, //c46
     62.0, //c55
      0.0, //c56
     49.0, //c66
  };
  seissol::model::Material material;
  seissol::model::setMaterial(materialVal, 22, &material); 

  std::complex<real> planeWaveOperator[9*9];
  model::getPlaneWaveOperator(material, m_kVec.data(), planeWaveOperator);

  using Matrix = Eigen::Matrix<std::complex<real>, 9, 9, Eigen::ColMajor>;
  using Vector = Eigen::Matrix<std::complex<real>, 9, 1, Eigen::ColMajor>;
  Matrix A(planeWaveOperator);
  Eigen::ComplexEigenSolver<Matrix> ces;
  ces.compute(A);

  auto eigenvalues = ces.eigenvalues();
  for (size_t i = 0; i < 9; ++i) {
    m_lambdaA[i] = eigenvalues(i,0);
  }

  Vector ic;
  for (size_t j = 0; j < 9; ++j) {
    ic(j) = 1.0;
  }

  auto eigenvectors = ces.eigenvectors();
  Vector amp = eigenvectors.colPivHouseholderQr().solve(ic);
  for (int j = 0; j < m_setVar; ++j) {
    m_varField.push_back(j);
    m_ampField.push_back(amp(j));
  }

  auto R = yateto::DenseTensorView<2,std::complex<real>>(m_eigenvectors, {9, 9});
  for (size_t j = 0; j < 9; ++j) {
    for (size_t i = 0; i < 9; ++i) {
      R(i,j) = eigenvectors(i,j);
    }
  }
}

void seissol::physics::Planarwave::evaluate(  double time,
                                              std::vector<std::array<double, 3>> const& points,
                                              yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const
{
  dofsQP.setZero();

  auto ra = yateto::DenseTensorView<2,std::complex<real>>(const_cast<std::complex<real>*>(m_eigenvectors), {9, 9});
  for (int v = 0; v < m_setVar; ++v) {
    for (int j = 0; j < 9; ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        dofsQP(i,j) += ra(j,m_varField[v]).real() * m_ampField[v].real()
                       * std::sin(m_kVec[0]*points[i][0]+m_kVec[1]*points[i][1]+m_kVec[2]*points[i][2] - m_lambdaA[m_varField[v]].real() * time + m_phase);
      }
    }
  }
}
