/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Local kernel of SeisSol.
 **/

#include "Kernels/Local.h"

#ifndef NDEBUG
#pragma message "compiling local kernel with assertions"
#endif

#include <yateto.h>


#include <array>
#include <cassert>
#include <stdint.h>
#include <Eigen/Dense>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop

#include <Kernels/common.hpp>
GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

void seissol::kernels::Local::setGlobalData(GlobalData const* global) {
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert( ((uintptr_t)global->stiffnessMatrices(stiffness)) % ALIGNMENT == 0 );
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed(flux)) % ALIGNMENT == 0 );
    assert( ((uintptr_t)global->changeOfBasisMatrices(flux)) % ALIGNMENT == 0 );
  }
#endif

  m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;

  m_nodalLfKrnlPrototype.project2nFaceTo3m = global->project2nFaceTo3m;

  m_projectKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
  m_projectRotatedKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
}

void seissol::kernels::Local::setDatReader( seissol::sourceterm::DAT* dat ) {
  m_dat = dat;
}

void seissol::kernels::Local::computeIntegral(real i_timeIntegratedDegreesOfFreedom[tensor::I::size()],
                                              LocalData& data,
                                              LocalTmp& tmp,
                                              // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                                              const CellMaterialData* materialData,
                                              CellBoundaryMapping const (*cellBoundaryMapping)[4],
                                              double time,
                                              double timeStepWidth) {
  assert(reinterpret_cast<uintptr_t>(i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(data.dofs) % ALIGNMENT == 0);

  kernel::volume volKrnl = m_volumeKernelPrototype;
  volKrnl.Q = data.dofs;
  volKrnl.I = i_timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = data.localIntegration.starMatrices[i];
  }

  // Optional source term
  set_ET(volKrnl, get_ptr_sourceMatrix(data.localIntegration.specific));

  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = data.dofs;
  lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs + tensor::Q::size();
  
  volKrnl.execute();

  for (int face = 0; face < 4; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {
      lfKrnl.AplusT = data.localIntegration.nApNm1[face];
      lfKrnl.execute(face);
    }

    alignas(ALIGNMENT) real dofsFaceBoundaryNodal[tensor::INodal::size()];
    auto nodalLfKrnl = m_nodalLfKrnlPrototype;
    nodalLfKrnl.Q = data.dofs;
    nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
    nodalLfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + tensor::I::size();
    nodalLfKrnl._prefetch.Q = data.dofs + tensor::Q::size();
    nodalLfKrnl.AminusT = data.neighboringIntegration.nAmNm1[face];

    // Include some boundary conditions here.
    switch (data.cellInformation.faceTypes[face]) {
    case FaceType::freeSurfaceGravity:
      {
      assert(cellBoundaryMapping != nullptr);
      assert(materialData != nullptr);
      auto* displ = tmp.nodalAvgDisplacements[face];
      auto displacement = init::INodalDisplacement::view::create(displ);
        auto applyFreeSurfaceBc = [&displacement, &materialData](
            const real*, // nodes are unused
            init::INodal::view::type& boundaryDofs) {
          for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
            const double rho = materialData->local.rho;
            const double g = 9.81; // [m/s^2]
            const double pressureAtBnd = -1 * rho * g * displacement(i);

            boundaryDofs(i,0) = 2 * pressureAtBnd - boundaryDofs(i,0);
            boundaryDofs(i,1) = 2 * pressureAtBnd - boundaryDofs(i,1);
            boundaryDofs(i,2) = 2 * pressureAtBnd - boundaryDofs(i,2);
          }
      };

      dirichletBoundary.evaluate(i_timeIntegratedDegreesOfFreedom,
                                 face,
                                 (*cellBoundaryMapping)[face],
                                 m_projectKrnlPrototype,
                                 applyFreeSurfaceBc,
                                 dofsFaceBoundaryNodal);

      nodalLfKrnl.execute(face);
      break;
      }
    case FaceType::dirichlet:
      {
      assert(cellBoundaryMapping != nullptr);
      auto* easiBoundaryMap = (*cellBoundaryMapping)[face].easiBoundaryMap;
      auto* easiBoundaryConstant = (*cellBoundaryMapping)[face].easiBoundaryConstant;
      assert(easiBoundaryConstant != nullptr);
      assert(easiBoundaryMap != nullptr);
      auto applyEasiBoundary = [easiBoundaryMap, easiBoundaryConstant](
          const real* nodes,
          init::INodal::view::type& boundaryDofs) {
        seissol::kernel::createEasiBoundaryGhostCells easiBoundaryKernel;
        easiBoundaryKernel.easiBoundaryMap = easiBoundaryMap;
        easiBoundaryKernel.easiBoundaryConstant = easiBoundaryConstant;
        easiBoundaryKernel.easiIdentMap = init::easiIdentMap::Values;
        easiBoundaryKernel.INodal = boundaryDofs.data();
        easiBoundaryKernel.execute();
      };

      // Compute boundary in [n, t_1, t_2] basis
      dirichletBoundary.evaluate(i_timeIntegratedDegreesOfFreedom,
				 face,
				 (*cellBoundaryMapping)[face],
				 m_projectRotatedKrnlPrototype,
				 applyEasiBoundary,
				 dofsFaceBoundaryNodal);

      // We do not need to rotate the boundary data back to the [x,y,z] basis
      // as we set the Tinv matrix to the identity matrix in the flux solver
      // See init. in CellLocalMatrices.initializeCellLocalMatrices!

      nodalLfKrnl.execute(face);
      break;
      }
      case FaceType::analytical:
      {
      assert(cellBoundaryMapping != nullptr);
      auto applyAnalyticalSolution = [materialData, this](const real* nodes,
                                                    double time,
                                                    init::INodal::view::type&,
                                                    init::INodal::view::type& boundaryDofs) {
          auto nodesVec = std::vector<std::array<double, 3>>{};
          int offset = 0;
          for (unsigned int i = 0; i < tensor::INodal::Shape[0]; ++i) {
            auto curNode = std::array<double, 3>{};
            curNode[0] = nodes[offset++];
            curNode[1] = nodes[offset++];
            curNode[2] = nodes[offset++];
            nodesVec.push_back(curNode);
          }
          assert(initConds != nullptr);
          // TODO(Lukas) Support multiple init. conds?
          assert(initConds->size() == 1);
          (*initConds)[0]->evaluate(time, nodesVec, *materialData, boundaryDofs);
      };

      dirichletBoundary.evaluateTimeDependent(i_timeIntegratedDegreesOfFreedom,
                                              face,
                                              (*cellBoundaryMapping)[face],
                                              m_projectKrnlPrototype,
                                              applyAnalyticalSolution,
                                              dofsFaceBoundaryNodal,
                                              time,
                                              timeStepWidth);
      nodalLfKrnl.execute(face);
      break;
      }
      case FaceType::velocityInlet:
      {
        assert(cellBoundaryMapping != nullptr);
        // Note: Everything happens in [n, t_1, t_2] basis
        // Comments of Dirichlet bc are also valid for this!
         auto applyInlet = [materialData, cellBoundaryMapping, face, this](const real* nodes,
                                               double time,
                                               init::INodal::view::type& boundaryDofsInterior,
                                               init::INodal::view::type& boundaryDofs) {
          int offset = 0;
          for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {

            // Only able to measure the pressure field during a finite time-interval [0, T]
            //  T = ENDTIME as defined in paramters.par, and thereby equal to the last time entry
            // of the receivers.
            const auto endtime = m_dat->endtime;

            // x, y, z: points on the surface
            const auto x = nodes[offset+0];
            const auto y = nodes[offset+1];
            const auto z = nodes[offset+2];
            offset += 3;
            
            Eigen::Vector3d position(x, y, z);

            auto H = [](double t) -> double {
              return t > 0 ? 1.0 : 0.0;
            };
            

            // T is a 9x9 rotation matrix
            // transforms from normal frame to global cartesian frame.
            auto T_inv = init::T::view::create((*cellBoundaryMapping)[face].TinvData);

            Eigen::MatrixXd T_inv_matrix(m_dat->q_dim, m_dat->q_dim);
            
            for (unsigned int row = 0; row < m_dat->q_dim; ++row) {
              for (unsigned int col = 0; col < m_dat->q_dim; ++col){
                  T_inv_matrix(row, col) = T_inv(row, col);
              }
            }

            Eigen::VectorXd q_cartesian = m_dat->getQ(position, endtime - time);

            Eigen::VectorXd q_normal;
            q_normal = T_inv_matrix * q_cartesian;


            for (unsigned int j = 0; j < m_dat->q_dim; ++j) {
              q_normal(j) = q_normal(j) * H(endtime - time);
            }


            boundaryDofs(i,0) = 2 * q_normal(0) - boundaryDofsInterior(i,0);
            boundaryDofs(i,1) = 2 * q_normal(1) - boundaryDofsInterior(i,1);
            boundaryDofs(i,2) = 2 * q_normal(2) - boundaryDofsInterior(i,2);
            boundaryDofs(i,3) = 2 * q_normal(3) - boundaryDofsInterior(i,3);
            boundaryDofs(i,4) = 2 * q_normal(4) - boundaryDofsInterior(i,4);
            boundaryDofs(i,5) = 2 * q_normal(5) - boundaryDofsInterior(i,5);
            // Either the stress tensor is specified at the boundary, or
            // the velocities:
            // boundaryDofs(i,6) = 2 * q_normal(6) - boundaryDofsInterior(i,6);
            // boundaryDofs(i,7) = 2 * q_normal(7) - boundaryDofsInterior(i,7);
            // boundaryDofs(i,8) = 2 * q_normal(8) - boundaryDofsInterior(i,8);


          }
        };

        dirichletBoundary.evaluateTimeDependent(i_timeIntegratedDegreesOfFreedom,
                                                face,
                                                (*cellBoundaryMapping)[face],
                                                m_projectRotatedKrnlPrototype,
                                                //m_projectKrnlPrototype,
                                                applyInlet,
                                                dofsFaceBoundaryNodal,
                                                time,
                                                timeStepWidth);
        nodalLfKrnl.execute(face);
        break;
      }
    default:
      // No boundary condition.
      break;
    }
  }
}

void seissol::kernels::Local::flopsIntegral(FaceType const i_faceTypes[4],
                                            unsigned int &o_nonZeroFlops,
                                            unsigned int &o_hardwareFlops)
{
  o_nonZeroFlops = seissol::kernel::volume::NonZeroFlops;
  o_hardwareFlops = seissol::kernel::volume::HardwareFlops;

  for( unsigned int face = 0; face < 4; ++face ) {
    // Local flux is executed for all faces that are not dynamic rupture.
    // For those cells, the flux is taken into account during the neighbor kernel.
    if (i_faceTypes[face] != FaceType::dynamicRupture) {
      o_nonZeroFlops += seissol::kernel::localFlux::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFlux::hardwareFlops(face);
    }

    // Take boundary condition flops into account.
    // Note that this only includes the flops of the kernels but not of the
    // boundary condition implementation.
    // The (probably incorrect) assumption is that they are negligible.
    switch (i_faceTypes[face]) {
    case FaceType::freeSurfaceGravity:
      o_nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
	seissol::kernel::projectToNodalBoundary::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
	seissol::kernel::projectToNodalBoundary::hardwareFlops(face);
      break;
    case FaceType::dirichlet:
      o_nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
	seissol::kernel::projectToNodalBoundaryRotated::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
	seissol::kernel::projectToNodalBoundary::hardwareFlops(face);
      break;
    case FaceType::analytical:
      o_nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
	CONVERGENCE_ORDER * seissol::kernel::updateINodal::NonZeroFlops;
      o_hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
	CONVERGENCE_ORDER * seissol::kernel::updateINodal::HardwareFlops;
      break;
    default:
      break;
    }
  }
}

unsigned seissol::kernels::Local::bytesIntegral()
{
  unsigned reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size();
  
  return reals * sizeof(real);
}
