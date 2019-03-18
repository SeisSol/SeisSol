#include <cmath>
#include <array>

#include <Kernels/precision.hpp>
#include <yateto/TensorView.h>

extern "C" {
  void initial_field_planarwave(double time, double x, double y, double z, double* variables) {

    std::fill_n(variables, NUMBER_OF_QUANTITIES, 0.0);
    real raData[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] = {0.0};
      
    auto ra = yateto::DenseTensorView<2,real>(raData, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});

    // Constants
    const auto rho0 = 1.0;
    const auto mu = 1.0;
    const auto lambda = 2.0;
    const auto n = std::array<real, 3>{0.577350269189626,
					 0.577350269189626,
					 0.577350269189626};
    const auto setVar = 2;
    const auto varField = std::array<int, 2>{1, 8};
    const auto ampField = std::array<real, 2>{1.0, 1.0};
    const auto kVec = std::array<real, 3>{6.283185307179587E-002,
					    6.283185307179587E-002,
					    6.283185307179587E-002};

    // Eigenvalues.
    const auto cp = std::sqrt((lambda+2*mu)/(rho0));
    const auto cs = std::sqrt((mu)/(rho0));

    const auto lambdaA = std::array<double, NUMBER_OF_QUANTITIES>{
      -cp, -cs, -cs, 0, 0, 0, cs, cs, cp
    };

    // Eigenvectors.
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

    real kVecNorm = 0.0;
    for (auto k: kVec) {
      kVecNorm += k * k;
    }
    kVecNorm = std::sqrt(kVecNorm);
      
    for (int i = 0; i < setVar; ++i) {
      const auto omega =  lambdaA[varField[i]] * kVecNorm;
      for (int j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
	variables[j] += ra(j,varField[i]) *
	  ampField[i] * std::sin(kVec[0]*x+kVec[1]*y+kVec[2]*z - omega * time);
      }
    }
  }
}
