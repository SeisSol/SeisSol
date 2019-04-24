#include <cmath>
#include <array>

#include <Kernels/precision.hpp>
#include <Physics/InitialField.h>
#include <yateto/TensorView.h>

seissol::physics::Planarwave::Planarwave(real phase)
  : m_setVar(2),
    m_varField{1,8},
    m_ampField{1.0, 1.0},
    m_kVec{6.283185307179587E-002, 6.283185307179587E-002, 6.283185307179587E-002},
    m_phase(phase)
{
  const auto rho0 = 1.0;
  const auto mu = 1.0;
  const auto lambda = 2.0;

  // Eigenvalues.
  const auto cp = std::sqrt((lambda+2*mu)/(rho0));
  const auto cs = std::sqrt((mu)/(rho0));

  double kVecNorm = 0.0;
  for (auto k: m_kVec) {
    kVecNorm += k * k;
  }
  kVecNorm = std::sqrt(kVecNorm);

  double omegaP = cp * kVecNorm;
  double omegaS = cs * kVecNorm;

  m_lambdaA = std::array<std::complex<real>, 9>{
    -omegaP, -omegaS, -omegaS, 0, 0, 0, omegaS, omegaS, omegaP
  };

  const auto n = std::array<real, 3>{0.577350269189626,
					 0.577350269189626,
					 0.577350269189626};

  auto ra = yateto::DenseTensorView<2,std::complex<real>>(m_eigenvectors, {9, 9});

  ra(0,0) = rho0*(-2*n[1]*n[1]*mu-2*n[2]*n[2]*mu+lambda+2*mu);
  ra(0,1) = -2*mu*n[1]*rho0*n[0]*n[0]*n[2];
  ra(0,2) = -2*n[2]*n[2]*mu*rho0*n[0]*n[0];
  ra(0,3) = -n[2]*n[2]*n[0];
  ra(0,4) = 0;
  ra(0,5) = -n[1]*n[0]*n[2];
  ra(0,6) = 2*n[2]*n[2]*mu*rho0*n[0]*n[0];
  ra(0,7) = 2*mu*n[1]*rho0*n[0]*n[0]*n[2];
  ra(0,8) = -rho0*(-2*n[1]*n[1]*mu-2*n[2]*n[2]*mu+lambda+2*mu);
    
  ra(1,0) = rho0*(2*n[1]*n[1]*mu+lambda);
  ra(1,1) = 2*mu*n[1]*rho0*n[0]*n[0]*n[2];
  ra(1,2) = 0;
  ra(1,3) = 0;
  ra(1,4) = -n[2]*n[2]/n[1]*n[0]*n[0];
  ra(1,5) = -1/n[1]*n[0]*n[0]*n[0]*n[2];
  ra(1,6) = 0;
  ra(1,7) = -2*mu*n[1]*rho0*n[0]*n[0]*n[2];
  ra(1,8) = -rho0*(2*n[1]*n[1]*mu+lambda);
    
  ra(2,0) = rho0*(2*n[2]*n[2]*mu+lambda);
  ra(2,1) = 0;
  ra(2,2) = 2*n[2]*n[2]*mu*rho0*n[0]*n[0];
  ra(2,3) = -n[0]*n[0]*n[0];
  ra(2,4) = -n[1]*n[0]*n[0];
  ra(2,5) = 0;
  ra(2,6) = -2*n[2]*n[2]*mu*rho0*n[0]*n[0];
  ra(2,7) = 0;
  ra(2,8) = -rho0*(2*n[2]*n[2]*mu+lambda);
    
  ra(3,0) = 2*n[1]*mu*n[0]*rho0;
  ra(3,1) = -mu*rho0*n[0]*(2*n[1]*n[1]+n[2]*n[2]-1)*n[2];
  ra(3,2) = -mu*n[1]*n[2]*n[2]*rho0*n[0];
  ra(3,3) = 0;
  ra(3,4) = 0;
  ra(3,5) = n[0]*n[0]*n[2];
  ra(3,6) = mu*n[1]*n[2]*n[2]*rho0*n[0];
  ra(3,7) = mu*rho0*n[0]*(2*n[1]*n[1]+n[2]*n[2]-1)*n[2];
  ra(3,8) = -2*n[1]*mu*n[0]*rho0;
    
  ra(4,0) = 2*mu*n[1]*rho0*n[2];
  ra(4,1) = n[2]*n[2]*mu*rho0*n[0]*n[0];
  ra(4,2) = mu*n[1]*rho0*n[0]*n[0]*n[2];
  ra(4,3) = 0;
  ra(4,4) = n[0]*n[0]*n[2];
  ra(4,5) = 0;
  ra(4,6) = -mu*n[1]*rho0*n[0]*n[0]*n[2];
  ra(4,7) = -n[2]*n[2]*mu*rho0*n[0]*n[0];
  ra(4,8) = -2*mu*n[1]*rho0*n[2];
    
  ra(5,0) = 2*mu*n[0]*rho0*n[2];
  ra(5,1) = -mu*n[1]*n[2]*n[2]*rho0*n[0];
  ra(5,2) = mu*rho0*n[0]*(-2*n[2]*n[2]-n[1]*n[1]+1)*n[2];
  ra(5,3) = n[0]*n[0]*n[2];
  ra(5,4) = 0;
  ra(5,5) = 0;
  ra(5,6) = -mu*rho0*n[0]*(-2*n[2]*n[2]-n[1]*n[1]+1)*n[2];
  ra(5,7) = mu*n[1]*n[2]*n[2]*rho0*n[0];
  ra(5,8) = -2*mu*n[0]*rho0*n[2];
    
  ra(6,0) = n[0]*std::sqrt(rho0*(lambda+2*mu));
  ra(6,1) = -n[1]*n[0]*n[2]*std::sqrt(rho0*mu);
  ra(6,2) = -n[2]*n[2]*n[0]*std::sqrt(rho0*mu);
  ra(6,3) = 0;
  ra(6,4) = 0;
  ra(6,5) = 0;
  ra(6,6) = -n[2]*n[2]*n[0]*std::sqrt(rho0*mu);
  ra(6,7) = -n[1]*n[0]*n[2]*std::sqrt(rho0*mu);
  ra(6,8) = n[0]*std::sqrt(rho0*(lambda+2*mu));
    
  ra(7,0) = n[1]*std::sqrt(rho0*(lambda+2*mu));
  ra(7,1) = n[0]*n[0]*n[2]*std::sqrt(rho0*mu);
  ra(7,2) = 0;
  ra(7,3) = 0;
  ra(7,4) = 0;
  ra(7,5) = 0;
  ra(7,6) = 0;
  ra(7,7) = n[0]*n[0]*n[2]*std::sqrt(rho0*mu);
  ra(7,8) = n[1]*std::sqrt(rho0*(lambda+2*mu));
    
  ra(8,0) = std::sqrt(rho0*(lambda+2*mu))*n[2];
  ra(8,1) = 0;
  ra(8,2) = n[0]*n[0]*n[2]*std::sqrt(rho0*mu);
  ra(8,3) = 0;
  ra(8,4) = 0;
  ra(8,5) = 0;
  ra(8,6) = n[0]*n[0]*n[2]*std::sqrt(rho0*mu);
  ra(8,7) = 0;
  ra(8,8) = std::sqrt(rho0*(lambda+2*mu))*n[2];
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
