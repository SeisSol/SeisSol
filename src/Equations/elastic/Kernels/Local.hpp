/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * Volume kernel of SeisSol.
 **/

#ifndef WAVEPROP_KERNEL_LOCAL_CK_H_
#define WAVEPROP_KERNEL_LOCAL_CK_H_

#include <yateto.h>

#include <array>
#include <cassert>
#include <stdint.h>
#include "GravitationalFreeSurfaceBC.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop
#include <memory>
#include "Common/configtensor.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop
#include "Physics/InitialField.h"
#include "Equations/Local.hpp"

#ifdef ACL_DEVICE
#include <device.h>
#endif
template <typename Config>
struct GlobalData;
namespace seissol::waveprop::kernel::local {
template <typename Config>
class Local<
    Config,
    std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevski ||
                     Config::MaterialT::Solver ==
                         seissol::model::LocalSolver::SpaceTimePredictorPoroelastic>> {
  using RealT = typename Config::RealT;

  protected:
  static void checkGlobalData(GlobalData<Config> const* global, size_t alignment) {
#ifndef NDEBUG
    for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
      assert(((uintptr_t)global->stiffnessMatrices(stiffness)) % alignment == 0);
    }
    for (unsigned flux = 0; flux < 4; ++flux) {
      assert(((uintptr_t)global->localChangeOfBasisMatricesTransposed(flux)) % alignment == 0);
      assert(((uintptr_t)global->changeOfBasisMatrices(flux)) % alignment == 0);
    }
#endif
  }

  typename Yateto<Config>::Kernel::volume m_volumeKernelPrototype;
  typename Yateto<Config>::Kernel::localFlux m_localFluxKernelPrototype;
  typename Yateto<Config>::Kernel::localFluxNodal m_nodalLfKrnlPrototype;

  typename Yateto<Config>::Kernel::projectToNodalBoundary m_projectKrnlPrototype;
  typename Yateto<Config>::Kernel::projectToNodalBoundaryRotated m_projectRotatedKrnlPrototype;

  seissol::kernels::DirichletBoundary<Config> dirichletBoundary;

#ifdef ACL_DEVICE
  typename Yateto<Config>::Kernel::gpu_volume deviceVolumeKernelPrototype;
  typename Yateto<Config>::Kernel::gpu_localFlux deviceLocalFluxKernelPrototype;
  typename Yateto<Config>::Kernel::gpu_localFluxNodal deviceNodalLfKrnlPrototype;
  typename Yateto<Config>::Kernel::gpu_projectToNodalBoundaryRotated
      deviceProjectRotatedKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  const std::vector<std::unique_ptr<physics::InitialField>>* initConds;

  public:
  virtual void setInitConds(decltype(initConds) initConds) { this->initConds = initConds; }

  physics::InitialField* getInitCond(size_t index) {
    const auto& condition = this->initConds->at(index);
    return condition.get();
  }

  void setHostGlobalData(GlobalData<Config> const* global) {
    checkGlobalData(global, Alignment);
    m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
    m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
    m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;

    m_nodalLfKrnlPrototype.project2nFaceTo3m = global->project2nFaceTo3m;

    m_projectKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
    m_projectRotatedKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
  }

  void setGlobalData(const CompoundGlobalData<Config>& global) {
    setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
    assert(global.onDevice != nullptr);
    const auto deviceAlignment = device.api->getGlobMemAlignment();
    checkGlobalData(global.onDevice, deviceAlignment);

    deviceVolumeKernelPrototype.kDivM = global.onDevice->stiffnessMatrices;
    deviceLocalFluxKernelPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
    deviceLocalFluxKernelPrototype.fMrT = global.onDevice->localChangeOfBasisMatricesTransposed;
    deviceNodalLfKrnlPrototype.project2nFaceTo3m = global.onDevice->project2nFaceTo3m;

    deviceProjectRotatedKrnlPrototype.V3mTo2nFace = global.onDevice->V3mTo2nFace;
#endif
  }

  template <typename LocalDataType>
  struct ApplyAnalyticalSolution {
    ApplyAnalyticalSolution(seissol::physics::InitialField* initCondition, LocalDataType& data)
        : initCondition(initCondition), localData(data) {}

    void operator()(const RealT* nodes,
                    double time,
                    typename seissol::Yateto<Config>::Init::INodal::view::type& boundaryDofs) {

      auto nodesVec = std::vector<std::array<double, 3>>{};
      int offset = 0;
      for (unsigned int i = 0; i < seissol::Yateto<Config>::Tensor::INodal::Shape[0]; ++i) {
        auto curNode = std::array<double, 3>{};
        curNode[0] = nodes[offset++];
        curNode[1] = nodes[offset++];
        curNode[2] = nodes[offset++];
        nodesVec.push_back(curNode);
      }

      assert(initCondition != nullptr);
      initCondition->evaluate(time, nodesVec, localData.materialData, boundaryDofs);
    }

private:
    seissol::physics::InitialField* initCondition{};
    LocalDataType& localData;
  };

  void computeIntegral(RealT i_timeIntegratedDegreesOfFreedom[Yateto<Config>::Tensor::I::size()],
                       LocalData<Config>& data,
                       LocalTmp<Config>& tmp,
                       // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                       const CellMaterialData* materialData,
                       CellBoundaryMapping const (*cellBoundaryMapping)[4],
                       double time,
                       double timeStepWidth) {
    assert(reinterpret_cast<uintptr_t>(i_timeIntegratedDegreesOfFreedom) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(data.dofs) % Alignment == 0);

    typename Yateto<Config>::Kernel::volume volKrnl = m_volumeKernelPrototype;
    volKrnl.Q = data.dofs;
    volKrnl.I = i_timeIntegratedDegreesOfFreedom;
    for (unsigned i = 0; i < yateto::numFamilyMembers<typename Yateto<Config>::Tensor::star>();
         ++i) {
      volKrnl.star(i) = data.localIntegration.starMatrices[i];
    }

    // Optional source term
    set_ET(volKrnl, get_ptr_sourceMatrix(data.localIntegration.specific));

    typename Yateto<Config>::Kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
    lfKrnl.Q = data.dofs;
    lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
    lfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + Yateto<Config>::Tensor::I::size();
    lfKrnl._prefetch.Q = data.dofs + Yateto<Config>::Tensor::Q::size();

    volKrnl.execute();

    for (int face = 0; face < 4; ++face) {
      // no element local contribution in the case of dynamic rupture boundary conditions
      if (data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {
        lfKrnl.AplusT = data.localIntegration.nApNm1[face];
        lfKrnl.execute(face);
      }

      alignas(Alignment) RealT dofsFaceBoundaryNodal[Yateto<Config>::Tensor::INodal::size()];
      auto nodalLfKrnl = m_nodalLfKrnlPrototype;
      nodalLfKrnl.Q = data.dofs;
      nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
      nodalLfKrnl._prefetch.I =
          i_timeIntegratedDegreesOfFreedom + Yateto<Config>::Tensor::I::size();
      nodalLfKrnl._prefetch.Q = data.dofs + Yateto<Config>::Tensor::Q::size();
      nodalLfKrnl.AminusT = data.neighboringIntegration.nAmNm1[face];

      // Include some boundary conditions here.
      switch (data.cellInformation.faceTypes[face]) {
      case FaceType::freeSurfaceGravity: {
        assert(cellBoundaryMapping != nullptr);
        assert(materialData != nullptr);
        auto* displ = tmp.nodalAvgDisplacements[face].data();
        auto displacement = Yateto<Config>::Init::averageNormalDisplacement::view::create(displ);
        auto applyFreeSurfaceBc =
            [&displacement,
             &materialData](const RealT*, // nodes are unused
                            typename Yateto<Config>::Init::INodal::view::type& boundaryDofs) {
              for (unsigned int i = 0; i < Yateto<Config>::Tensor::nodal::nodes2D::Shape[0]; ++i) {
                const double rho = materialData->local->rho;
                const double g = getGravitationalAcceleration(); // [m/s^2]
                const double pressureAtBnd = -1 * rho * g * displacement(i);

                boundaryDofs(i, 0) = 2 * pressureAtBnd - boundaryDofs(i, 0);
                boundaryDofs(i, 1) = 2 * pressureAtBnd - boundaryDofs(i, 1);
                boundaryDofs(i, 2) = 2 * pressureAtBnd - boundaryDofs(i, 2);
              }
            };

        dirichletBoundary.evaluate(i_timeIntegratedDegreesOfFreedom,
                                   face,
                                   (*cellBoundaryMapping)[face],
                                   m_projectRotatedKrnlPrototype,
                                   applyFreeSurfaceBc,
                                   dofsFaceBoundaryNodal);

        nodalLfKrnl.execute(face);
        break;
      }
      case FaceType::dirichlet: {
        assert(cellBoundaryMapping != nullptr);
        auto* easiBoundaryMap = (*cellBoundaryMapping)[face].easiBoundaryMap;
        auto* easiBoundaryConstant = (*cellBoundaryMapping)[face].easiBoundaryConstant;
        assert(easiBoundaryConstant != nullptr);
        assert(easiBoundaryMap != nullptr);
        auto applyEasiBoundary = [easiBoundaryMap, easiBoundaryConstant](
                                     const RealT* nodes,
                                     typename Yateto<Config>::Init::INodal::view::type&
                                         boundaryDofs) {
          typename seissol::Yateto<Config>::Kernel::createEasiBoundaryGhostCells easiBoundaryKernel;
          easiBoundaryKernel.easiBoundaryMap = easiBoundaryMap;
          easiBoundaryKernel.easiBoundaryConstant = easiBoundaryConstant;
          easiBoundaryKernel.easiIdentMap = Yateto<Config>::Init::easiIdentMap::Values;
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
      case FaceType::analytical: {
        assert(cellBoundaryMapping != nullptr);
        assert(initConds != nullptr);
        assert(initConds->size() == 1);
        ApplyAnalyticalSolution applyAnalyticalSolution(this->getInitCond(0), data);

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
      default:
        // No boundary condition.
        break;
      }
    }
  }

  void computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                              ConditionalMaterialTable& materialTable,
                              ConditionalIndicesTable& indicesTable,
                              typename kernels::LocalData<Config>::Loader& loader,
                              LocalTmp<Config>& tmp,
                              double timeStepWidth) {
#ifdef ACL_DEVICE
    // Volume integral
    ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    typename Yateto<Config>::Kernel::gpu_volume volKrnl = deviceVolumeKernelPrototype;
    typename Yateto<Config>::Kernel::gpu_localFlux localFluxKrnl = deviceLocalFluxKernelPrototype;

    const auto maxTmpMem = yateto::getMaxTmpMemRequired(volKrnl, localFluxKrnl);

    RealT* tmpMem = nullptr;
    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];

      unsigned maxNumElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
      volKrnl.numElements = maxNumElements;

      // volume kernel always contains more elements than any local one
      tmpMem = (RealT*)(device.api->getStackMemory(maxTmpMem * maxNumElements));

      volKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      volKrnl.I =
          const_cast<const RealT**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

      unsigned starOffset = 0;
      for (size_t i = 0; i < yateto::numFamilyMembers<typename Yateto<Config>::Tensor::star>();
           ++i) {
        volKrnl.star(i) =
            const_cast<const RealT**>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
        volKrnl.extraOffset_star(i) = starOffset;
        starOffset += Yateto<Config>::Tensor::star::size(i);
      }
      volKrnl.linearAllocator.initialize(tmpMem);
      volKrnl.streamPtr = device.api->getDefaultStream();
      volKrnl.execute();
    }

    // Local Flux Integral
    for (unsigned face = 0; face < 4; ++face) {
      key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

      if (dataTable.find(key) != dataTable.end()) {
        auto& entry = dataTable[key];
        localFluxKrnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
        localFluxKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
        localFluxKrnl.I =
            const_cast<const RealT**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
        localFluxKrnl.AplusT =
            const_cast<const RealT**>(entry.get(inner_keys::Wp::Id::AplusT)->getDeviceDataPtr());
        localFluxKrnl.linearAllocator.initialize(tmpMem);
        localFluxKrnl.streamPtr = device.api->getDefaultStream();
        localFluxKrnl.execute(face);
      }

      ConditionalKey fsgKey(
          *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, face);
      if (dataTable.find(fsgKey) != dataTable.end()) {
        auto nodalAvgDisplacements =
            dataTable[fsgKey].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();
        auto rhos = materialTable[fsgKey].get(inner_keys::Material::Id::Rho)->getDeviceDataPtr();
        local_flux::aux::FreeSurfaceGravity freeSurfaceGravityBc;
        freeSurfaceGravityBc.g = getGravitationalAcceleration();
        freeSurfaceGravityBc.rhos = rhos;
        freeSurfaceGravityBc.displacementDataPtrs = nodalAvgDisplacements;
        dirichletBoundary.evaluateOnDevice(face,
                                           fsgKey,
                                           deviceProjectRotatedKrnlPrototype,
                                           deviceNodalLfKrnlPrototype,
                                           freeSurfaceGravityBc,
                                           dataTable,
                                           device);
      }

      ConditionalKey dirichletKey(
          *KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
      if (dataTable.find(dirichletKey) != dataTable.end()) {
        auto easiBoundaryMapPtrs =
            dataTable[dirichletKey].get(inner_keys::Wp::Id::EasiBoundaryMap)->getDeviceDataPtr();
        auto easiBoundaryConstantPtrs = dataTable[dirichletKey]
                                            .get(inner_keys::Wp::Id::EasiBoundaryConstant)
                                            ->getDeviceDataPtr();

        local_flux::aux::EasiBoundary easiBoundaryBc;
        easiBoundaryBc.easiBoundaryMapPtrs = easiBoundaryMapPtrs;
        easiBoundaryBc.easiBoundaryConstantPtrs = easiBoundaryConstantPtrs;

        dirichletBoundary.evaluateOnDevice(face,
                                           dirichletKey,
                                           deviceProjectRotatedKrnlPrototype,
                                           deviceNodalLfKrnlPrototype,
                                           easiBoundaryBc,
                                           dataTable,
                                           device);
      }
    }
    if (tmpMem != nullptr) {
      device.api->popStackMemory();
    }
#else
    assert(false && "no implementation provided");
#endif
  }

  void evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                      ConditionalIndicesTable& indicesTable,
                                      kernels::LocalData<Config>::Loader& loader,
                                      double time,
                                      double timeStepWidth) {

#ifdef ACL_DEVICE
    for (unsigned face = 0; face < 4; ++face) {
      ConditionalKey analyticalKey(
          *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
      if (indicesTable.find(analyticalKey) != indicesTable.end()) {
        auto idofsPtrs = dataTable[analyticalKey].get(inner_keys::Wp::Id::Idofs)->getHostData();

        auto cellIds =
            indicesTable[analyticalKey].get(inner_keys::Indices::Id::Cells)->getHostData();
        const size_t numElements = cellIds.size();

        for (unsigned index{0}; index < numElements; ++index) {
          auto cellId = cellIds[index];
          auto data = loader.entry(cellId);

          alignas(Alignment) RealT dofsFaceBoundaryNodal[Yateto<Config>::Tensor::INodal::size()];

          assert(initConds != nullptr);
          assert(initConds->size() == 1);
          ApplyAnalyticalSolution applyAnalyticalSolution(this->getInitCond(0), data);

          dirichletBoundary.evaluateTimeDependent(idofsPtrs[index],
                                                  face,
                                                  data.boundaryMapping[face],
                                                  m_projectKrnlPrototype,
                                                  applyAnalyticalSolution,
                                                  dofsFaceBoundaryNodal,
                                                  time,
                                                  timeStepWidth);

          auto nodalLfKrnl = this->m_nodalLfKrnlPrototype;
          nodalLfKrnl.Q = data.dofs;
          nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
          nodalLfKrnl.AminusT = data.neighboringIntegration.nAmNm1[face];
          nodalLfKrnl.execute(face);
        }
      }
    }
#else
    assert(false && "no implementation provided");
#endif // ACL_DEVICE
  }

  void flopsIntegral(FaceType const i_faceTypes[4],
                     unsigned int& o_nonZeroFlops,
                     unsigned int& o_hardwareFlops) {
    o_nonZeroFlops = seissol::Yateto<Config>::Kernel::volume::NonZeroFlops;
    o_hardwareFlops = seissol::Yateto<Config>::Kernel::volume::HardwareFlops;

    for (unsigned int face = 0; face < 4; ++face) {
      // Local flux is executed for all faces that are not dynamic rupture.
      // For those cells, the flux is taken into account during the neighbor kernel.
      if (i_faceTypes[face] != FaceType::dynamicRupture) {
        o_nonZeroFlops += seissol::Yateto<Config>::Kernel::localFlux::nonZeroFlops(face);
        o_hardwareFlops += seissol::Yateto<Config>::Kernel::localFlux::hardwareFlops(face);
      }

      // Take boundary condition flops into account.
      // Note that this only includes the flops of the kernels but not of the
      // boundary condition implementation.
      // The (probably incorrect) assumption is that they are negligible.
      switch (i_faceTypes[face]) {
      case FaceType::freeSurfaceGravity:
        o_nonZeroFlops +=
            seissol::Yateto<Config>::Kernel::localFluxNodal::nonZeroFlops(face) +
            seissol::Yateto<Config>::Kernel::projectToNodalBoundary::nonZeroFlops(face);
        o_hardwareFlops +=
            seissol::Yateto<Config>::Kernel::localFluxNodal::hardwareFlops(face) +
            seissol::Yateto<Config>::Kernel::projectToNodalBoundary::hardwareFlops(face);
        break;
      case FaceType::dirichlet:
        o_nonZeroFlops +=
            seissol::Yateto<Config>::Kernel::localFluxNodal::nonZeroFlops(face) +
            seissol::Yateto<Config>::Kernel::projectToNodalBoundaryRotated::nonZeroFlops(face);
        o_hardwareFlops +=
            seissol::Yateto<Config>::Kernel::localFluxNodal::hardwareFlops(face) +
            seissol::Yateto<Config>::Kernel::projectToNodalBoundary::hardwareFlops(face);
        break;
      case FaceType::analytical:
        o_nonZeroFlops +=
            seissol::Yateto<Config>::Kernel::localFluxNodal::nonZeroFlops(face) +
            Config::ConvergenceOrder * seissol::Yateto<Config>::Kernel::updateINodal::NonZeroFlops;
        o_hardwareFlops +=
            seissol::Yateto<Config>::Kernel::localFluxNodal::hardwareFlops(face) +
            Config::ConvergenceOrder * seissol::Yateto<Config>::Kernel::updateINodal::HardwareFlops;
        break;
      default:
        break;
      }
    }
  }

  unsigned bytesIntegral() {
    unsigned reals = 0;

    // star matrices load
    reals += yateto::computeFamilySize<typename Yateto<Config>::Kernel::star>();
    // flux solvers
    reals += 4 * Yateto<Config>::Kernel::AplusT::size();

    // DOFs write
    reals += Yateto<Config>::Kernel::Q::size();

    return reals * sizeof(RealT);
  }
};

} // namespace seissol::waveprop::kernel::local
#endif
