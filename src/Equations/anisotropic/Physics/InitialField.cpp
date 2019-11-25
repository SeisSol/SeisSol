#include <cmath>
#include <array>
#include <Eigen/Eigenvalues>

#include <Kernels/precision.hpp>
#include <Physics/InitialField.h>
#include <Model/Setup.h>
#include <Model/common.hpp>
#include <yateto/TensorView.h>


seissol::physics::AnisotropicPlanarwave::AnisotropicPlanarwave(real phase)
  : m_kVec1{M_PI, 0.0, 0.0},
    m_kVec2{0.0, M_PI, 0.0},
    m_kVec3{0.0, 0.0, M_PI},
    m_phase(phase)
{ 

  m_pw1 = seissol::physics::Planarwave(phase, m_kVec1);
  m_pw2 = seissol::physics::Planarwave(phase, m_kVec2);
  m_pw3 = seissol::physics::Planarwave(phase, m_kVec3);
}

void seissol::physics::AnisotropicPlanarwave::evaluate(  double time,
                                              std::vector<std::array<double, 3>> const& points,
                                              yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const
{
  dofsQP.setZero();
  
  m_pw1.evaluate(time, points, dofsQP);
  m_pw2.evaluate(time, points, dofsQP);
  m_pw3.evaluate(time, points, dofsQP);
}
