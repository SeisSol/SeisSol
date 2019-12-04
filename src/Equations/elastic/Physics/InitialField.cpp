#include <cmath>
#include <array>

#include <Kernels/precision.hpp>
#include <Physics/InitialField.h>
#include <yateto/TensorView.h>

seissol::physics::Planarwave::Planarwave(double phase)
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

  m_lambdaA = std::array<std::complex<double>, 9>{
    -omegaP, -omegaS, -omegaS, 0, 0, 0, omegaS, omegaS, omegaP
  };

  const auto n = std::array<double, 3>{0.577350269189626,
					 0.577350269189626,
					 0.577350269189626};

  auto ra = yateto::DenseTensorView<2,std::complex<double>>(m_eigenvectors, {9, 9});

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

void seissol::physics::Planarwave::evaluate(double time,
                                            std::vector<std::array<double, 3>> const& points,
					    const CellMaterialData& materialData,
                                            yateto::DenseTensorView<2,real,unsigned>& dofsQP) const
{
  dofsQP.setZero();

  auto ra = yateto::DenseTensorView<2,std::complex<double>>(const_cast<std::complex<double>*>(m_eigenvectors), {9, 9});
  for (int v = 0; v < m_setVar; ++v) {
    for (int j = 0; j < 9; ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        dofsQP(i,j) += ra(j,m_varField[v]).real() * m_ampField[v].real()
                       * std::sin(m_kVec[0]*points[i][0]+m_kVec[1]*points[i][1]+m_kVec[2]*points[i][2] - m_lambdaA[m_varField[v]].real() * time + m_phase);
      }
    }
  }
}

void seissol::physics::ScholteWave::evaluate(double time,
					     std::vector<std::array<double, 3>> const& points,
					     const CellMaterialData& materialData,
					     yateto::DenseTensorView<2,real,unsigned>& dofsQp) const {
  const real omega = 2.0 * std::acos(-1);

  for (size_t i = 0; i < points.size(); ++i) {
    const auto& x = points[i];
    const bool isAcousticPart = std::abs(materialData.local.mu) < std::numeric_limits<real>::epsilon();
    const auto x_1 = x[0];
    const auto x_3 = x[2];
    const auto t = time;
    if (isAcousticPart) {
      dofsQp(i,0) = 0.35944997730200889*std::pow(omega, 2)*std::exp(-0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // sigma_xx
      dofsQp(i,1) = 0.35944997730200889*std::pow(omega, 2)*std::exp(-0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // sigma_yy
      dofsQp(i,2) = 0.35944997730200889*std::pow(omega, 2)*std::exp(-0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // sigma_zz
      dofsQp(i,3) = 0; // sigma_xy
      dofsQp(i,4) = 0; // sigma_yz
      dofsQp(i,5) = 0; // sigma_xz
      dofsQp(i,6) = -0.50555429848461109*std::pow(omega, 2)*std::exp(-0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // u
      dofsQp(i,7) = 0; // v
      dofsQp(i,8) = 0.35550086150929727*std::pow(omega, 2)*std::exp(-0.98901344820674908*omega*x_3)*std::cos(omega*t - 1.406466352506808*omega*x_1); // w
    } else {
      dofsQp(i,0) = -2.7820282741590652*std::pow(omega, 2)*std::exp(0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1) + 3.5151973269883681*std::pow(omega, 2)*std::exp(1.2825031256883821*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // sigma_xx
      dofsQp(i,1) = -6.6613381477509402e-16*std::pow(omega, 2)*std::exp(0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1) + 0.27315475753283058*std::pow(omega, 2)*std::exp(1.2825031256883821*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // sigma_yy
      dofsQp(i,2) = 2.7820282741590621*std::pow(omega, 2)*std::exp(0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1) - 2.4225782968570462*std::pow(omega, 2)*std::exp(1.2825031256883821*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // sigma_zz
      dofsQp(i,3) = 0; // sigma_xy
      dofsQp(i,4) = 0; // sigma_yz
      dofsQp(i,5) = -2.956295201467618*std::pow(omega, 2)*std::exp(0.98901344820674908*omega*x_3)*std::cos(omega*t - 1.406466352506808*omega*x_1) + 2.9562952014676029*std::pow(omega, 2)*std::exp(1.2825031256883821*omega*x_3)*std::cos(omega*t - 1.406466352506808*omega*x_1); // sigma_xz
      dofsQp(i,6) = 0.98901344820675241*std::pow(omega, 2)*std::exp(0.98901344820674908*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1) - 1.1525489264912381*std::pow(omega, 2)*std::exp(1.2825031256883821*omega*x_3)*std::sin(omega*t - 1.406466352506808*omega*x_1); // u
      dofsQp(i,7) = 0; // v
      dofsQp(i,8) = 1.406466352506812*std::pow(omega, 2)*std::exp(0.98901344820674908*omega*x_3)*std::cos(omega*t - 1.406466352506808*omega*x_1) - 1.050965490997515*std::pow(omega, 2)*std::exp(1.2825031256883821*omega*x_3)*std::cos(omega*t - 1.406466352506808*omega*x_1); // w
    }
  }
}

void seissol::physics::SnellsLaw::evaluate(double time,
					   std::vector<std::array<double, 3>> const& points,
					   const CellMaterialData& materialData,
					   yateto::DenseTensorView<2,real,unsigned>& dofsQp) const {
  const double pi = std::acos(-1);
  const double omega = 2.0 * pi;

  for (size_t i = 0; i < points.size(); ++i) {
    const auto &x = points[i];
    const bool isAcousticPart = std::abs(materialData.local.mu) < std::numeric_limits<real>::epsilon();

    const auto x_1 = x[0];
    const auto x_3 = x[2];
    const auto t = time;
    if (isAcousticPart) {
      dofsQp(i,0) = 1.0*omega*std::sin(omega*t - omega*(0.19866933079506119*x_1 + 0.98006657784124163*x_3)) + 0.48055591432167399*omega*std::sin(omega*t - omega*(0.19866933079506149*x_1 - 0.98006657784124152*x_3)); // sigma_xx
      dofsQp(i,1) = 1.0*omega*std::sin(omega*t - omega*(0.19866933079506119*x_1 + 0.98006657784124163*x_3)) + 0.48055591432167399*omega*std::sin(omega*t - omega*(0.19866933079506149*x_1 - 0.98006657784124152*x_3)); // sigma_yy
      dofsQp(i,2) = 1.0*omega*std::sin(omega*t - omega*(0.19866933079506119*x_1 + 0.98006657784124163*x_3)) + 0.48055591432167399*omega*std::sin(omega*t - omega*(0.19866933079506149*x_1 - 0.98006657784124152*x_3)); // sigma_zz
      dofsQp(i,3) = 0; // sigma_xy
      dofsQp(i,4) = 0; // sigma_yz
      dofsQp(i,5) = 0; // sigma_xz
      dofsQp(i,6) = -0.19866933079506119*omega*std::sin(omega*t - omega*(0.19866933079506119*x_1 + 0.98006657784124163*x_3)) - 0.095471721907895893*omega*std::sin(omega*t - omega*(0.19866933079506149*x_1 - 0.98006657784124152*x_3)); // u
      dofsQp(i,7) = 0; // v
      dofsQp(i,8) = -0.98006657784124163*omega*std::sin(omega*t - omega*(0.19866933079506119*x_1 + 0.98006657784124163*x_3)) + 0.47097679041061191*omega*std::sin(omega*t - omega*(0.19866933079506149*x_1 - 0.98006657784124152*x_3)); // w
    } else {
      dofsQp(i,0) = -0.59005639909185559*omega*std::sin(omega*t - 1.0/2.0*omega*(0.39733866159012299*x_1 + 0.91767204817721759*x_3)) + 0.55554011463785213*omega*std::sin(omega*t - 1.0/3.0*omega*(0.59600799238518454*x_1 + 0.8029785009656123*x_3)); // sigma_xx
      dofsQp(i,1) = 0.14460396298676709*omega*std::sin(omega*t - 1.0/3.0*omega*(0.59600799238518454*x_1 + 0.8029785009656123*x_3)); // sigma_yy
      dofsQp(i,2) = 0.59005639909185559*omega*std::sin(omega*t - 1.0/2.0*omega*(0.39733866159012299*x_1 + 0.91767204817721759*x_3)) + 0.89049951522981918*omega*std::sin(omega*t - 1.0/3.0*omega*(0.59600799238518454*x_1 + 0.8029785009656123*x_3)); // sigma_zz
      dofsQp(i,3) = 0; // sigma_xy
      dofsQp(i,4) = 0; // sigma_yz
      dofsQp(i,5) = -0.55363837274201066*omega*std::sin(omega*t - 1.0/2.0*omega*(0.39733866159012299*x_1 + 0.91767204817721759*x_3)) + 0.55363837274201*omega*std::sin(omega*t - 1.0/3.0*omega*(0.59600799238518454*x_1 + 0.8029785009656123*x_3)); // sigma_xz
      dofsQp(i,6) = 0.37125533967075403*omega*std::sin(omega*t - 1.0/2.0*omega*(0.39733866159012299*x_1 + 0.91767204817721759*x_3)) - 0.2585553530120539*omega*std::sin(omega*t - 1.0/3.0*omega*(0.59600799238518454*x_1 + 0.8029785009656123*x_3)); // u
      dofsQp(i,7) = 0; // v
      dofsQp(i,8) = -0.16074816713222639*omega*std::sin(omega*t - 1.0/2.0*omega*(0.39733866159012299*x_1 + 0.91767204817721759*x_3)) - 0.34834162029840349*omega*std::sin(omega*t - 1.0/3.0*omega*(0.59600799238518454*x_1 + 0.8029785009656123*x_3)); // w
    }
  }
}

void seissol::physics::Ocean::evaluate(double time,
                                       std::vector<std::array<double, 3>> const& points,
                                       const CellMaterialData& materialData,
                                       yateto::DenseTensorView<2,real,unsigned>& dofsQp) const {
  for (size_t i = 0; i < points.size(); ++i) {
    const auto x = points[i][0];
    const auto y = points[i][1];
    const auto z = points[i][2];
    const auto t = time;

    const double g = 9.81; // m/s
    const double pi = std::acos(-1);
    assert(materialData.local.mu == 0); // has to be acoustic
    const double rho = materialData.local.rho;

    const double k_x = pi/100; // 1/m
    const double k_y = pi/100; // 1/m
    constexpr double k_star = 0.0444284459948;
    constexpr double omega = 0.276857520383318;

    const auto B = g * k_star/(omega*omega);
    const auto pressure = -std::sin(k_x*x)*std::sin(k_y*y)*std::sin(omega*t)*(std::sinh(k_star*z)
        + B * std::cosh(k_star*z));

    dofsQp(i,0) = pressure;
    dofsQp(i,1) = pressure;
    dofsQp(i,2) = pressure;
    dofsQp(i,3) = 0.0;
    dofsQp(i,4) = 0.0;
    dofsQp(i,5) = 0.0;
    dofsQp(i,6) = (k_x/(omega*rho))*cos(k_x*x)*std::sin(k_y*y)*cos(omega*t)*(std::sinh(k_star*z)
        + B * std::cosh(k_star*z));
    dofsQp(i,7) = (k_y/(omega*rho))*std::sin(k_x*x)*cos(k_y*y)*cos(omega*t)*(std::sinh(k_star*z)
        + B * std::cosh(k_star*z));
    dofsQp(i,8) = (k_star/(omega*rho))*std::sin(k_x*x)*std::sin(k_y*y)*cos(omega*t)*(std::cosh(k_star*z)
        + B * std::sinh(k_star*z));
  }
}
