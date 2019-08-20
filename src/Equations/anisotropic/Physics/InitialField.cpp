#include <cmath>
#include <array>
#include <Eigen/Eigenvalues>

#include <Kernels/precision.hpp>
#include <Equations/anisotropic/Physics/InitialField.h>
#include <Model/Setup.h>
#include <yateto/TensorView.h>

seissol::physics::Planarwave::Planarwave(real phase)
  : m_setVar(9),
    m_kVec1{M_PI, 0.0, 0.0},
    m_kVec2{0.0, M_PI, 0.0},
    m_kVec3{0.0, 0.0, M_PI},
    m_phase(phase)
{ 
  double materialVal[22] = { 
//use this for the anisotropic test case
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
//use this for the isotropic test case 
 //     1.0, //rho
 //     4.0, //c11
 //     2.0, //c12
 //     2.0, //c13
 //     0.0, //c14
 //     0.0, //c15
 //     0.0, //c16
 //     4.0, //c22
 //     2.0, //c23
 //     0.0, //c24
 //     0.0, //c25
 //     0.0, //c26
 //     4.0, //c33
 //     0.0, //c34
 //     0.0, //c35
 //     0.0, //c36
 //     1.0, //c44
 //     0.0, //c45
 //     0.0, //c46
 //     1.0, //c55
 //     0.0, //c56
 //     1.0, //c66
  };
  seissol::model::Material material;
  seissol::model::setMaterial(materialVal, 22, &material); 

  std::complex<real> planeWaveOperator_1[9*9];
  std::complex<real> planeWaveOperator_2[9*9];
  std::complex<real> planeWaveOperator_3[9*9];

  model::getPlaneWaveOperator(material, m_kVec1.data(), planeWaveOperator_1);
  model::getPlaneWaveOperator(material, m_kVec2.data(), planeWaveOperator_2);
  model::getPlaneWaveOperator(material, m_kVec3.data(), planeWaveOperator_3);

  using Matrix = Eigen::Matrix<std::complex<real>, 9, 9, Eigen::ColMajor>;
  Matrix op_1(planeWaveOperator_1);
  Matrix op_2(planeWaveOperator_2);
  Matrix op_3(planeWaveOperator_3);
  Eigen::ComplexEigenSolver<Matrix> ces1;
  Eigen::ComplexEigenSolver<Matrix> ces2;
  Eigen::ComplexEigenSolver<Matrix> ces3;
  ces1.compute(op_1);
  ces2.compute(op_2);
  ces3.compute(op_3);

  auto eigenvalues1 = ces1.eigenvalues();
  auto eigenvalues2 = ces2.eigenvalues();
  auto eigenvalues3 = ces3.eigenvalues();
  for (size_t i = 0; i < 9; ++i) {
    m_lambda1[i] = eigenvalues1(i,0);
    m_lambda2[i] = eigenvalues2(i,0);
    m_lambda3[i] = eigenvalues3(i,0);
  }
  logInfo() << "Wave frequencies: " << m_lambda1[0] << ", " << m_lambda2[0] << ", " << m_lambda3[0];

  auto eigenvectors1 = ces1.eigenvectors();
  auto eigenvectors2 = ces2.eigenvectors();
  auto eigenvectors3 = ces3.eigenvectors();

  auto R1 = yateto::DenseTensorView<2,std::complex<real>>(m_eigenvectors1, {9, 9});
  auto R2 = yateto::DenseTensorView<2,std::complex<real>>(m_eigenvectors2, {9, 9});
  auto R3 = yateto::DenseTensorView<2,std::complex<real>>(m_eigenvectors3, {9, 9});
  for (size_t j = 0; j < 9; ++j) {
    for (size_t i = 0; i < 9; ++i) {
      R1(i,j) = eigenvectors1(i,j);
      R2(i,j) = eigenvectors2(i,j);
      R3(i,j) = eigenvectors3(i,j);
    }
  }
}

void seissol::physics::Planarwave::evaluate(  double time,
                                              std::vector<std::array<double, 3>> const& points,
                                              yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const
{
  dofsQP.setZero();

  auto r1 = yateto::DenseTensorView<2,std::complex<real>>(const_cast<std::complex<real>*>(m_eigenvectors1), {9, 9});
  auto r2 = yateto::DenseTensorView<2,std::complex<real>>(const_cast<std::complex<real>*>(m_eigenvectors2), {9, 9});
  auto r3 = yateto::DenseTensorView<2,std::complex<real>>(const_cast<std::complex<real>*>(m_eigenvectors3), {9, 9});
  for (int v = 0; v < m_setVar; ++v) {
    for (int j = 0; j < 9; ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
          dofsQP(i,j) += r1(j,v).real() * 100.0
                       * std::sin(m_kVec1[0]*points[i][0]+m_kVec1[1]*points[i][1]+m_kVec1[2]*points[i][2] - m_lambda1[v].real() * time + m_phase)
                       + r2(j,v).real() * 100.0
                       * std::sin(m_kVec2[0]*points[i][0]+m_kVec2[1]*points[i][1]+m_kVec2[2]*points[i][2] - m_lambda2[v].real() * time + m_phase)
                       + r3(j,v).real() * 100.0
                       * std::sin(m_kVec3[0]*points[i][0]+m_kVec3[1]*points[i][1]+m_kVec3[2]*points[i][2] - m_lambda3[v].real() * time + m_phase);
      }
    }
  }
}
