#include <cmath>
#include <array>
#include <Eigen/Eigenvalues>

#include <Kernels/precision.hpp>
#include <Physics/InitialField.h>
#include <Model/Setup.h>
#include <Model/common.hpp>
#include <yateto/TensorView.h>


seissol::physics::AnisotropicPlanarwave::AnisotropicPlanarwave(real phase)
  : m_kVec({{{M_PI, 0.0, 0.0},
             {0.0, M_PI, 0.0},
             {0.0, 0.0, M_PI}}}),
    m_phase(phase)
{ 

  for (int i = 0; i < 3; i++) {
    m_pw.at(i) = seissol::physics::Planarwave(phase, m_kVec.at(i));
  }
}

void seissol::physics::AnisotropicPlanarwave::evaluate( double time,
                                                        std::vector<std::array<double, 3>> const& points,
                                                        const CellMaterialData& materialData,
                                                        yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const
{
  dofsQP.setZero();
 
  real dofsPW_data[tensor::dofsQP::size()];
  yateto::DenseTensorView<2,real,unsigned> dofsPW = init::dofsQP::view::create(dofsPW_data);

  for (int pw = 0; pw < 3; pw++) {
    m_pw.at(pw).evaluate(time, points, materialData, dofsPW);
    for (unsigned j = 0; j < dofsQP.shape(1); ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        dofsQP(i,j) += dofsPW(i,j);
      }
    }
  }
}
