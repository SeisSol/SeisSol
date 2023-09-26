/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
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
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Setup of SeisSol's cell local matrices.
 **/

#include "CellLocalMatrices.h"

#include <cassert>

#include <Initializer/ParameterDB.h>
#include "Initializer/MemoryManager.h"
#include <Numerical_aux/Transformation.h>
#include <Equations/Setup.h>
#include <Model/common.hpp>
#include <Geometry/MeshTools.h>
#include <generated_code/tensor.h>
#include <generated_code/kernel.h>
#include <utils/logger.h>
#include "Common/cellconfigconv.hpp"
#ifdef ACL_DEVICE
#include <device.h>
#endif

template <typename Config>
static void setStarMatrix( typename Config::RealT* i_AT,
                    typename Config::RealT* i_BT,
                    typename Config::RealT* i_CT,
                    typename Config::RealT  i_grad[3],
                    typename Config::RealT* o_starMatrix )
{
  for (unsigned idx = 0; idx < seissol::Yateto<Config>::Tensor::star::size(0); ++idx) {
    o_starMatrix[idx] = i_grad[0] * i_AT[idx];
  }

  for (unsigned idx = 0; idx < seissol::Yateto<Config>::Tensor::star::size(1); ++idx) {
    o_starMatrix[idx] += i_grad[1] * i_BT[idx];
  }

  for (unsigned idx = 0; idx < seissol::Yateto<Config>::Tensor::star::size(2); ++idx) {
    o_starMatrix[idx] += i_grad[2] * i_CT[idx];
  }
}

void seissol::initializers::initializeCellLocalMatrices( seissol::geometry::MeshReader const&      i_meshReader,
                                                         ClusterLTSForest&      forest,
                                                         TimeStepping const&    timeStepping )
{
  std::vector<Element> const& elements = i_meshReader.getElements();
  std::vector<Vertex> const& vertices = i_meshReader.getVertices();

  assert(seissol::Yateto<Config>::Tensor::AplusT::Shape[0] == seissol::Yateto<Config>::Tensor::AminusT::Shape[0]);
  assert(seissol::Yateto<Config>::Tensor::AplusT::Shape[1] == seissol::Yateto<Config>::Tensor::AminusT::Shape[1]);

  forest.visit([&](auto&& ltsview) {
    using Config = typename std::decay_t<decltype(ltsview.lts)>::ConfigT;
    using RealT = typename Config::RealT;
    using MaterialT = typename Config::MaterialT;

  assert(LayerMask(Ghost) == ltsview.lts.material.mask);
  assert(LayerMask(Ghost) == ltsview.lts.localIntegration.mask);
  assert(LayerMask(Ghost) == ltsview.lts.neighboringIntegration.mask);

  for (LTSTree::leaf_iterator it = ltsview.tree.beginLeaf(LayerMask(Ghost)); it != ltsview.tree.endLeaf(); ++it) {
    CellMaterialData*           material                = it->var(ltsview.lts.material);
    MaterialT*           materialData                = it->var(ltsview.lts.materialData);
    LocalIntegrationData<Config>*       localIntegration        = it->var(ltsview.lts.localIntegration);
    NeighboringIntegrationData<Config>* neighboringIntegration  = it->var(ltsview.lts.neighboringIntegration);
    CellLocalInformation*       cellInformation         = it->var(ltsview.lts.cellInformation);
    SecondaryCellLocalInformation*       secondaryCellInformation         = it->var(ltsview.lts.secondaryCellInformation);

#ifdef _OPENMP
  #pragma omp parallel
    {
#endif
    RealT ATData[Yateto<Config>::Tensor::star::size(0)];
    RealT ATtildeData[Yateto<Config>::Tensor::star::size(0)];
    RealT BTData[Yateto<Config>::Tensor::star::size(1)];
    RealT CTData[Yateto<Config>::Tensor::star::size(2)];
    auto AT = Yateto<Config>::Init::star::template view<0>::create(ATData);
    // AT with elastic parameters in local coordinate system, used for flux kernel
    auto ATtilde = Yateto<Config>::Init::star::template view<0>::create(ATtildeData);
    auto BT = Yateto<Config>::Init::star::template view<0>::create(BTData);
    auto CT = Yateto<Config>::Init::star::template view<0>::create(CTData);

    RealT TData[seissol::Yateto<Config>::Tensor::T::size()];
    RealT TinvData[seissol::Yateto<Config>::Tensor::Tinv::size()];
    auto T = Yateto<Config>::Init::T::view::create(TData);
    auto Tinv = Yateto<Config>::Init::Tinv::view::create(TinvData);

    RealT QgodLocalData[Yateto<Config>::Tensor::QgodLocal::size()];
    RealT QgodNeighborData[Yateto<Config>::Tensor::QgodNeighbor::size()];
    auto QgodLocal = Yateto<Config>::Init::QgodLocal::view::create(QgodLocalData);
    auto QgodNeighbor = Yateto<Config>::Init::QgodNeighbor::view::create(QgodNeighborData);
    
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
      unsigned clusterId = secondaryCellInformation[cell].clusterId;
      auto timeStepWidth = timeStepping.globalCflTimeStepWidths[clusterId];
      unsigned meshId = secondaryCellInformation[cell].meshId;

      RealT x[4];
      RealT y[4];
      RealT z[4];
      RealT gradXi[3];
      RealT gradEta[3];
      RealT gradZeta[3];

      // Iterate over all 4 vertices of the tetrahedron
      for (unsigned vertex = 0; vertex < 4; ++vertex) {
        VrtxCoords const& coords = vertices[ elements[meshId].vertices[vertex] ].coords;
        x[vertex] = coords[0];
        y[vertex] = coords[1];
        z[vertex] = coords[2];
      }

      seissol::transformations::tetrahedronGlobalToReferenceJacobian( x, y, z, gradXi, gradEta, gradZeta );

      seissol::model::getTransposedCoefficientMatrix( materialData[cell], 0, AT );
      seissol::model::getTransposedCoefficientMatrix( materialData[cell], 1, BT );
      seissol::model::getTransposedCoefficientMatrix( materialData[cell], 2, CT );
      setStarMatrix<Config>(ATData, BTData, CTData, gradXi, localIntegration[cell].starMatrices[0]);
      setStarMatrix<Config>(ATData, BTData, CTData, gradEta, localIntegration[cell].starMatrices[1]);
      setStarMatrix<Config>(ATData, BTData, CTData, gradZeta, localIntegration[cell].starMatrices[2]);

      double volume = MeshTools::volume(elements[meshId], vertices);

      for (unsigned side = 0; side < 4; ++side) {
        VrtxCoords normal;
        VrtxCoords tangent1;
        VrtxCoords tangent2;
        MeshTools::normalAndTangents(elements[meshId], side, vertices, normal, tangent1, tangent2);
        double surface = MeshTools::surface(normal);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);

        RealT NLocalData[6*6];
        seissol::model::getBondMatrix(normal, tangent1, tangent2, NLocalData);

        std::visit([&](auto&& config) {
          using NeighborConfig = std::decay_t<decltype(config)>;
          using NeighborMaterialT = typename NeighborConfig::MaterialT;
          const auto* neighborMaterialData = material[cell].neighbor[side];
          const auto* neighborMaterialCasted = dynamic_cast<const NeighborMaterialT*>(neighborMaterialData);
          auto neighborMaterialConverged = MaterialT(*neighborMaterialCasted);
          // TODO(David): implement others here...
          seissol::model::getTransposedGodunovState(  seissol::model::getRotatedMaterialCoefficients(NLocalData, materialData[cell]),
                                                      seissol::model::getRotatedMaterialCoefficients(NLocalData, neighborMaterialConverged),
                                                      cellInformation[cell].faceTypes[side],
                                                      QgodLocal,
                                                      QgodNeighbor );
        }, ConfigInstances[cellInformation.neighborConfigIds[side]]);
        seissol::model::getTransposedCoefficientMatrix( seissol::model::getRotatedMaterialCoefficients(NLocalData, materialData[cell]), 0, ATtilde );

        // Calculate transposed T instead
        seissol::model::getFaceRotationMatrix<MaterialT>(normal, tangent1, tangent2, T, Tinv);

        // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
        // must be subtracted.
        RealT fluxScale = -2.0 * surface / (6.0 * volume);

        typename Yateto<Config>::Kernel::computeFluxSolverLocal localKrnl;
        localKrnl.fluxScale = fluxScale;
        localKrnl.AplusT = localIntegration[cell].nApNm1[side];
        localKrnl.QgodLocal = QgodLocalData;
        localKrnl.T = TData;
        localKrnl.Tinv = TinvData;
        localKrnl.star(0) = ATtildeData;
        localKrnl.execute();
        
        typename Yateto<Config>::Kernel::computeFluxSolverNeighbor neighKrnl;
        neighKrnl.fluxScale = fluxScale;
        neighKrnl.AminusT = neighboringIntegration[cell].nAmNm1[side];
        neighKrnl.QgodNeighbor = QgodNeighborData;
        neighKrnl.T = TData;
        neighKrnl.Tinv = TinvData;
        neighKrnl.star(0) = ATtildeData;
        if (cellInformation[cell].faceTypes[side] == FaceType::dirichlet ||
            cellInformation[cell].faceTypes[side] == FaceType::freeSurfaceGravity) {
          // Already rotated!
          neighKrnl.Tinv = Yateto<Config>::Init::identityT::Values;
        }
        neighKrnl.execute();
      }

      seissol::model::initializeSpecificLocalData(  materialData[cell],
                                                    timeStepWidth,
                                                    &localIntegration[cell].specific );

      seissol::model::initializeSpecificNeighborData( materialData[cell],
                                                      &neighboringIntegration[cell].specific );

    }
#ifdef _OPENMP
    }
#endif
  }

  });
}

void surfaceAreaAndVolume(  seissol::geometry::MeshReader const&      i_meshReader,
                            unsigned               meshId,
                            unsigned               side,
                            double*                surfaceArea,
                            double*                volume )
{
  std::vector<Vertex> const& vertices = i_meshReader.getVertices();
  std::vector<Element> const& elements = i_meshReader.getElements();

  VrtxCoords normal;
  VrtxCoords tangent1;
  VrtxCoords tangent2;
  MeshTools::normalAndTangents(elements[meshId], side, vertices, normal, tangent1, tangent2);

  *volume = MeshTools::volume(elements[meshId], vertices);
  *surfaceArea = MeshTools::surface(normal);
}

void seissol::initializers::initializeBoundaryMappings(const seissol::geometry::MeshReader& i_meshReader,
                                                       const EasiBoundary* easiBoundary,
                                                       ClusterLTSForest&      forest) {
  std::vector<Element> const& elements = i_meshReader.getElements();
  std::vector<Vertex> const& vertices = i_meshReader.getVertices();

  forest.visit([&](auto&& ltsview) {
    using Config = typename std::decay_t<decltype(ltsview.lts)>::ConfigT;
    using RealT = typename Config::RealT;
    using MaterialT = typename Config::MaterialT;
  for (LTSTree::leaf_iterator it = ltsview.tree.beginLeaf(LayerMask(Ghost)); it != ltsview.tree.endLeaf(); ++it) {
    auto* cellInformation = it->var(ltsview.lts.cellInformation);
    auto* secondaryCellInformation = it->var(ltsview.lts.secondaryCellInformation);
    auto* boundary = it->var(ltsview.lts.boundaryMapping);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
      const auto& element = elements[secondaryCellInformation[cell].meshId];
      double const* coords[4];
      for (unsigned v = 0; v < 4; ++v) {
        coords[v] = vertices[ element.vertices[ v ] ].coords;
      }
      for (unsigned side = 0; side < 4; ++side) {
        if (cellInformation[cell].faceTypes[side] != FaceType::freeSurfaceGravity
            && cellInformation[cell].faceTypes[side] != FaceType::dirichlet
            && cellInformation[cell].faceTypes[side] != FaceType::analytical) {
          continue;
        }
        // Compute nodal points in global coordinates for each side.
        RealT nodesReferenceData[Yateto<Config>::Tensor::nodal::nodes2D::Size];
        std::copy_n(Yateto<Config>::Init::nodal::nodes2D::Values,
                    Yateto<Config>::Tensor::nodal::nodes2D::Size,
                    nodesReferenceData);
        auto nodesReference = Yateto<Config>::Init::nodal::nodes2D::view::create(nodesReferenceData);
        auto nodes = boundary[cell][side].nodes;
        assert(nodes != nullptr);
        auto offset = 0;
        for (unsigned int i = 0; i < Yateto<Config>::Tensor::nodal::nodes2D::Shape[0]; ++i) {
          double nodeReference[2];
          nodeReference[0] = nodesReference(i,0);
          nodeReference[1] = nodesReference(i,1);
          // Compute the global coordinates for the nodal points.
          double xiEtaZeta[3], xyz[3];
          seissol::transformations::chiTau2XiEtaZeta(side,
                                                     nodeReference,
                                                     xiEtaZeta);
          seissol::transformations::tetrahedronReferenceToGlobal(coords[0],
                                                                 coords[1],
                                                                 coords[2],
                                                                 coords[3],
                                                                 xiEtaZeta,
                                                                 xyz);
          nodes[offset++] = xyz[0];
          nodes[offset++] = xyz[1];
          nodes[offset++] = xyz[2];
        }

        // Compute map that rotates to normal aligned coordinate system.
        RealT* TData = boundary[cell][side].TData;
        RealT* TinvData = boundary[cell][side].TinvData;
        assert(TData != nullptr);
        assert(TinvData != nullptr);
        auto T = Yateto<Config>::Init::T::view::create(TData);
        auto Tinv = Yateto<Config>::Init::Tinv::view::create(TinvData);

        VrtxCoords normal;
        VrtxCoords tangent1;
        VrtxCoords tangent2;
        MeshTools::normalAndTangents(element, side, vertices, normal, tangent1, tangent2);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);
        seissol::model::getFaceRotationMatrix<MaterialT>(normal, tangent1, tangent2, T, Tinv);

        // Evaluate easi boundary condition matrices if needed
        RealT* easiBoundaryMap = boundary[cell][side].easiBoundaryMap;
        RealT* easiBoundaryConstant = boundary[cell][side].easiBoundaryConstant;
        assert(easiBoundaryMap != nullptr);
        assert(easiBoundaryConstant != nullptr);
        if (cellInformation[cell].faceTypes[side] == FaceType::dirichlet) {
          easiBoundary->query(nodes, easiBoundaryMap, easiBoundaryConstant);
        } else {
          // Boundary should not be evaluated
          std::fill_n(easiBoundaryMap,
                      seissol::Yateto<Config>::Tensor::easiBoundaryMap::size(),
                      std::numeric_limits<RealT>::signaling_NaN());
          std::fill_n(easiBoundaryConstant,
                      seissol::Yateto<Config>::Tensor::easiBoundaryConstant::size(),
                      std::numeric_limits<RealT>::signaling_NaN());
        }

      }
    }
  }

  });
}

void seissol::initializers::initializeDynamicRuptureMatrices( seissol::geometry::MeshReader const&      i_meshReader,                                                    
                                            ClusterLTSForest& forest,
                                            const ClusterBackmap&                   i_ltsLut,
                                            DynRupLTSForest& dynrupforest,
                                            const DynrupBackmap&              ltsFaceToMeshFace,
                                            const GlobalDataStorage&      global)
{

  std::vector<Fault> const& fault = i_meshReader.getFault();
  std::vector<Element> const& elements = i_meshReader.getElements();

  forest.visitTwo([&](auto&& ltsview, auto&& dynrup) {
    using Config = typename std::decay_t<decltype(ltsview.lts)>::ConfigT;
    using RealT = typename Config::RealT;
    using MaterialT = typename Config::MaterialT;

    RealT TData[Yateto<Config>::Tensor::T::size()];
  RealT TinvData[Yateto<Config>::Tensor::Tinv::size()];
  RealT APlusData[Yateto<Config>::Tensor::star::size(0)];
  RealT AMinusData[Yateto<Config>::Tensor::star::size(0)];

  CellDRMapping (*drMapping)[4] = ltsview.tree.var(ltsview.lts.drMapping);
  CellMaterialData* material = ltsview.tree.var(ltsview.lts.material);
  RealT** derivatives = ltsview.tree.var(ltsview.lts.derivatives);
  RealT* (*faceNeighbors)[4] = ltsview.tree.var(ltsview.lts.faceNeighbors);
  CellLocalInformation* cellInformation = ltsview.tree.var(ltsview.lts.cellInformation);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (LTSTree::leaf_iterator it = dynrup.tree.beginLeaf(LayerMask(Ghost)); it != dynrup.tree.endLeaf(); ++it) {
    RealT**                                timeDerivativePlus                                        = it->var(dynrup.lts.timeDerivativePlus);
    RealT**                                timeDerivativeMinus                                       = it->var(dynrup.lts.timeDerivativeMinus);
    RealT                                (*imposedStatePlus)[Yateto<Config>::Tensor::QInterpolated::size()]          = it->var(dynrup.lts.imposedStatePlus);
    RealT                                (*imposedStateMinus)[Yateto<Config>::Tensor::QInterpolated::size()]         = it->var(dynrup.lts.imposedStateMinus);
    DRGodunovData<Config>*                        godunovData                                               = it->var(dynrup.lts.godunovData);
    RealT                                (*fluxSolverPlus)[Yateto<Config>::Tensor::fluxSolver::size()]               = it->var(dynrup.lts.fluxSolverPlus);
    RealT                                (*fluxSolverMinus)[Yateto<Config>::Tensor::fluxSolver::size()]              = it->var(dynrup.lts.fluxSolverMinus);
    DRFaceInformation*                    faceInformation                                           = it->var(dynrup.lts.faceInformation);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                            = it->var(dynrup.lts.waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                           = it->var(dynrup.lts.waveSpeedsMinus);
    seissol::dr::ImpedancesAndEta<Config>*        impAndEta                                                 = it->var(dynrup.lts.impAndEta);


#ifdef _OPENMP
  #pragma omp parallel for private(TData, TinvData, APlusData, AMinusData) schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      assert(fault[meshFace].element >= 0 || fault[meshFace].neighborElement >= 0);

      /// Face information
      faceInformation[ltsFace].meshFace = meshFace;
      faceInformation[ltsFace].plusSide = fault[meshFace].side;
      faceInformation[ltsFace].minusSide = fault[meshFace].neighborSide;
      if (fault[meshFace].element >= 0) {
        faceInformation[ltsFace].faceRelation = elements[ fault[meshFace].element ].sideOrientations[ fault[meshFace].side ] + 1;
        faceInformation[ltsFace].plusSideOnThisRank = true;
      } else {
        /// \todo check if this is correct
        faceInformation[ltsFace].faceRelation = elements[ fault[meshFace].neighborElement ].sideOrientations[ fault[meshFace].neighborSide ] + 1;
        faceInformation[ltsFace].plusSideOnThisRank = false;
      }

      /// Look for time derivative mapping in all duplicates
      int derivativesMeshId;
      int derivativesSide;
      if (fault[meshFace].element >= 0) {
        derivativesMeshId = fault[meshFace].element;
        derivativesSide = faceInformation[ltsFace].plusSide;
      } else {
        derivativesMeshId = fault[meshFace].neighborElement;
        derivativesSide = faceInformation[ltsFace].minusSide;
      }
      RealT* timeDerivative1 = NULL;
      RealT* timeDerivative2 = NULL;
      for (unsigned duplicate = 0; duplicate < Lut::MaxDuplicates; ++duplicate) {
        unsigned ltsId = i_ltsLut->ltsId(ltsview.lts.cellInformation.mask, derivativesMeshId, duplicate);
        if (timeDerivative1 == NULL && (cellInformation[ltsId].ltsSetup >> 9)%2 == 1) {
          timeDerivative1 = derivatives[ i_ltsLut->ltsId(ltsview.lts.derivatives.mask, derivativesMeshId, duplicate) ];
        }
        if (timeDerivative2 == NULL && (cellInformation[ltsId].ltsSetup >> derivativesSide)%2 == 1) {
          timeDerivative2 = faceNeighbors[ i_ltsLut->ltsId(ltsview.lts.faceNeighbors.mask, derivativesMeshId, duplicate) ][ derivativesSide ];
        }
      }

      assert(timeDerivative1 != NULL && timeDerivative2 != NULL);

      if (fault[meshFace].element >= 0) {
        timeDerivativePlus[ltsFace] = timeDerivative1;
        timeDerivativeMinus[ltsFace] = timeDerivative2;
      } else {
        timeDerivativePlus[ltsFace] = timeDerivative2;
        timeDerivativeMinus[ltsFace] = timeDerivative1;
      }

      assert(timeDerivativePlus[ltsFace] != NULL && timeDerivativeMinus[ltsFace] != NULL);

      /// DR mapping for elements
      for (unsigned duplicate = 0; duplicate < Lut::MaxDuplicates; ++duplicate) {
        unsigned plusLtsId = (fault[meshFace].element >= 0)          ? i_ltsLut->ltsId(ltsview.lts.drMapping.mask, fault[meshFace].element, duplicate) : std::numeric_limits<unsigned>::max();
        unsigned minusLtsId = (fault[meshFace].neighborElement >= 0) ? i_ltsLut->ltsId(ltsview.lts.drMapping.mask, fault[meshFace].neighborElement, duplicate) : std::numeric_limits<unsigned>::max();

        assert(duplicate != 0 || plusLtsId != std::numeric_limits<unsigned>::max() || minusLtsId != std::numeric_limits<unsigned>::max());

        if (plusLtsId != std::numeric_limits<unsigned>::max()) {
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
{
          CellDRMapping& mapping = drMapping[plusLtsId][ faceInformation[ltsFace].plusSide ];
          mapping.side = faceInformation[ltsFace].plusSide;
          mapping.faceRelation = 0;
          mapping.godunov = &imposedStatePlus[ltsFace][0];
          mapping.fluxSolver = &fluxSolverPlus[ltsFace][0];
}
        }
        if (minusLtsId != std::numeric_limits<unsigned>::max()) {
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
{
          CellDRMapping& mapping = drMapping[minusLtsId][ faceInformation[ltsFace].minusSide ];
          mapping.side = faceInformation[ltsFace].minusSide;
          mapping.faceRelation = faceInformation[ltsFace].faceRelation;
          mapping.godunov = &imposedStateMinus[ltsFace][0];
          mapping.fluxSolver = &fluxSolverMinus[ltsFace][0];
}
        }
      }

      /// Transformation matrix
      auto T = Yateto<Config>::Init::T::view::create(TData);
      auto Tinv = Yateto<Config>::Init::Tinv::view::create(TinvData);
      seissol::model::getFaceRotationMatrix<MaterialT>(fault[meshFace].normal, fault[meshFace].tangent1, fault[meshFace].tangent2, T, Tinv);

      /// Materials
      seissol::model::Material* plusMaterial;
      seissol::model::Material* minusMaterial;
      unsigned plusLtsId = (fault[meshFace].element >= 0)          ? i_ltsLut->ltsId(ltsview.lts.material.mask, fault[meshFace].element) : std::numeric_limits<unsigned>::max();
      unsigned minusLtsId = (fault[meshFace].neighborElement >= 0) ? i_ltsLut->ltsId(ltsview.lts.material.mask, fault[meshFace].neighborElement) : std::numeric_limits<unsigned>::max();

      assert(plusLtsId != std::numeric_limits<unsigned>::max() || minusLtsId != std::numeric_limits<unsigned>::max());

      if (plusLtsId != std::numeric_limits<unsigned>::max()) {
        plusMaterial = material[plusLtsId].local;
        minusMaterial = material[plusLtsId].neighbor[ faceInformation[ltsFace].plusSide ];
      } else {
        assert(minusLtsId != std::numeric_limits<unsigned>::max());
        plusMaterial = material[minusLtsId].neighbor[ faceInformation[ltsFace].minusSide ];
        minusMaterial = material[minusLtsId].local;
      }

      /// Wave speeds and Coefficient Matrices
      auto APlus = Yateto<Config>::Init::star::template view<0>::create(APlusData);
      auto AMinus = Yateto<Config>::Init::star::template view<0>::create(AMinusData);
      
      waveSpeedsPlus[ltsFace].density = plusMaterial->rho;
      waveSpeedsMinus[ltsFace].density = minusMaterial->rho;
      waveSpeedsPlus[ltsFace].pWaveVelocity = plusMaterial->getPWaveSpeed();
      waveSpeedsPlus[ltsFace].sWaveVelocity = plusMaterial->getSWaveSpeed();
      waveSpeedsMinus[ltsFace].pWaveVelocity = minusMaterial->getPWaveSpeed();
      waveSpeedsMinus[ltsFace].sWaveVelocity = minusMaterial->getSWaveSpeed();

      //calculate Impedances Z and eta
      impAndEta[ltsFace].zp = (waveSpeedsPlus[ltsFace].density * waveSpeedsPlus[ltsFace].pWaveVelocity);
      impAndEta[ltsFace].zpNeig = (waveSpeedsMinus[ltsFace].density * waveSpeedsMinus[ltsFace].pWaveVelocity);
      impAndEta[ltsFace].zs = (waveSpeedsPlus[ltsFace].density * waveSpeedsPlus[ltsFace].sWaveVelocity);
      impAndEta[ltsFace].zsNeig = (waveSpeedsMinus[ltsFace].density * waveSpeedsMinus[ltsFace].sWaveVelocity);

      impAndEta[ltsFace].invZp = 1/impAndEta[ltsFace].zp;
      impAndEta[ltsFace].invZpNeig = 1/impAndEta[ltsFace].zpNeig;
      impAndEta[ltsFace].invZs = 1/impAndEta[ltsFace].zs;
      impAndEta[ltsFace].invZsNeig = 1/impAndEta[ltsFace].zsNeig;

      impAndEta[ltsFace].etaP = 1.0 / (1.0 / impAndEta[ltsFace].zp + 1.0 / impAndEta[ltsFace].zpNeig);
      impAndEta[ltsFace].invEtaS = 1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig;
      impAndEta[ltsFace].etaS = 1.0 / (1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig);

      seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<MaterialT*>(plusMaterial), 0, APlus);
      seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<MaterialT*>(minusMaterial), 0, AMinus);

      /// Traction matrices for "average" traction
      auto tractionPlusMatrix = Yateto<Config>::Init::tractionPlusMatrix::view::create(godunovData[ltsFace].tractionPlusMatrix);
      auto tractionMinusMatrix = Yateto<Config>::Init::tractionMinusMatrix::view::create(godunovData[ltsFace].tractionMinusMatrix);
      double ZpP = plusMaterial->rho * waveSpeedsPlus[ltsFace].pWaveVelocity;
      double ZsP = plusMaterial->rho * waveSpeedsPlus[ltsFace].sWaveVelocity;
      double ZpM = minusMaterial->rho * waveSpeedsMinus[ltsFace].pWaveVelocity;
      double ZsM = minusMaterial->rho * waveSpeedsMinus[ltsFace].sWaveVelocity;
      double etaP = ZpP*ZpM / (ZpP + ZpM);
      double etaS = ZsP*ZsM / (ZsP + ZsM);

      tractionPlusMatrix.setZero();
      tractionPlusMatrix(0,0) = etaP / ZpP;
      tractionPlusMatrix(3,1) = etaS / ZsP;
      tractionPlusMatrix(5,2) = etaS / ZsP;

      tractionMinusMatrix.setZero();
      tractionMinusMatrix(0,0) = etaP / ZpM;
      tractionMinusMatrix(3,1) = etaS / ZsM;
      tractionMinusMatrix(5,2) = etaS / ZsM;

      /// Transpose Tinv
      typename Yateto<Config>::Kernel::dynamicRupture::transposeTinv ttKrnl;
      ttKrnl.Tinv = TinvData;
      ttKrnl.TinvT = godunovData[ltsFace].TinvT;
      ttKrnl.execute();

      double plusSurfaceArea, plusVolume, minusSurfaceArea, minusVolume, surfaceArea;
      if (fault[meshFace].element >= 0) {
        surfaceAreaAndVolume( i_meshReader, fault[meshFace].element, fault[meshFace].side, &plusSurfaceArea, &plusVolume );
        surfaceArea = plusSurfaceArea;
      } else {
        /// Blow up solution on purpose if used by mistake
        plusSurfaceArea = 1.e99; plusVolume = 1.0;
      }
      if (fault[meshFace].neighborElement >= 0) {
        surfaceAreaAndVolume( i_meshReader, fault[meshFace].neighborElement, fault[meshFace].neighborSide, &minusSurfaceArea, &minusVolume );
        surfaceArea = minusSurfaceArea;
      } else {
        /// Blow up solution on purpose if used by mistake
        minusSurfaceArea = 1.e99; minusVolume = 1.0;
      }
      godunovData[ltsFace].doubledSurfaceArea = 2.0 * surfaceArea;

      typename Yateto<Config>::Kernel::dynamicRupture::rotateFluxMatrix krnl;
      krnl.T = TData;

      krnl.fluxSolver = fluxSolverPlus[ltsFace];
      krnl.fluxScale = -2.0 * plusSurfaceArea / (6.0 * plusVolume);
      krnl.star(0) = APlusData;
      krnl.execute();

      krnl.fluxSolver = fluxSolverMinus[ltsFace];
      krnl.fluxScale = 2.0 * minusSurfaceArea / (6.0 * minusVolume);
      krnl.star(0) = AMinusData;
      krnl.execute();
    }

    layerLtsFaceToMeshFace += it->getNumberOfCells();
  }

  }, dynrupforest);
}

void seissol::initializers::copyCellMatricesToDevice(ClusterLTSForest& forest,
                                                     DynRupLTSForest& dynrupforest,
                                                     BoundaryLTSForest& boundaryforest) {
#ifdef ACL_DEVICE
  forest.visit([&](auto&& ltsview) {
    // Byte-copy of element compute-static data from the host to device
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    const std::vector<size_t > &variableSizes = ltsview.tree.getVariableSizes();

    device.api->copyTo(ltsview.tree.var(ltsview.lts.localIntegrationOnDevice),
                      ltsview.tree.var(ltsview.lts.localIntegration),
                      variableSizes[ltsview.lts.localIntegration.index]);

    device.api->copyTo(ltsview.tree.var(ltsview.lts.neighIntegrationOnDevice),
                      ltsview.tree.var(ltsview.lts.neighboringIntegration),
                      variableSizes[ltsview.lts.neighboringIntegration.index]);
  });
#endif // ACL_DEVICE
}
