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
#include <Numerical_aux/Transformation.h>
#include <Equations/Setup.h>
#include <Model/common.hpp>
#include <Geometry/MeshTools.h>
#include <generated_code/tensor.h>
#include <generated_code/kernel.h>
#include <utils/logger.h>

void setStarMatrix( real* i_AT,
                    real* i_BT,
                    real* i_CT,
                    real  i_grad[3],
                    real* o_starMatrix )
{
  for (unsigned idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
    o_starMatrix[idx] = i_grad[0] * i_AT[idx];
  }

  for (unsigned idx = 0; idx < seissol::tensor::star::size(1); ++idx) {
    o_starMatrix[idx] += i_grad[1] * i_BT[idx];
  }

  for (unsigned idx = 0; idx < seissol::tensor::star::size(2); ++idx) {
    o_starMatrix[idx] += i_grad[2] * i_CT[idx];
  }
}

void seissol::initializers::initializeCellLocalMatrices( MeshReader const&      i_meshReader,
                                                         LTSTree*               io_ltsTree,
                                                         LTS*                   i_lts,
                                                         Lut*                   i_ltsLut )
{
  std::vector<Element> const& elements = i_meshReader.getElements();
  std::vector<Vertex> const& vertices = i_meshReader.getVertices();

  assert(seissol::tensor::AplusT::Shape[0] == seissol::tensor::AminusT::Shape[0]);
  assert(seissol::tensor::AplusT::Shape[1] == seissol::tensor::AminusT::Shape[1]);

  unsigned* ltsToMesh = i_ltsLut->getLtsToMeshLut(i_lts->material.mask);

  assert(LayerMask(Ghost) == i_lts->material.mask);
  assert(LayerMask(Ghost) == i_lts->localIntegration.mask);
  assert(LayerMask(Ghost) == i_lts->neighboringIntegration.mask);

  assert(ltsToMesh      == i_ltsLut->getLtsToMeshLut(i_lts->localIntegration.mask));
  assert(ltsToMesh      == i_ltsLut->getLtsToMeshLut(i_lts->neighboringIntegration.mask));

  for (LTSTree::leaf_iterator it = io_ltsTree->beginLeaf(LayerMask(Ghost)); it != io_ltsTree->endLeaf(); ++it) {
    CellMaterialData*           material                = it->var(i_lts->material);
    LocalIntegrationData*       localIntegration        = it->var(i_lts->localIntegration);
    NeighboringIntegrationData* neighboringIntegration  = it->var(i_lts->neighboringIntegration);
    CellLocalInformation*       cellInformation         = it->var(i_lts->cellInformation);

#ifdef _OPENMP
  #pragma omp parallel
    {
#endif
    real ATData[tensor::star::size(0)];
    real ATtildeData[tensor::star::size(0)];
    real BTData[tensor::star::size(1)];
    real CTData[tensor::star::size(2)];
    auto AT = init::star::view<0>::create(ATData);
    // AT with elastic parameters in local coordinate system, used for flux kernel
    auto ATtilde = init::star::view<0>::create(ATtildeData);
    auto BT = init::star::view<0>::create(BTData);
    auto CT = init::star::view<0>::create(CTData);

    real TData[seissol::tensor::T::size()];
    real TinvData[seissol::tensor::Tinv::size()];
    auto T = init::T::view::create(TData);
    auto Tinv = init::Tinv::view::create(TinvData);

    real QgodLocalData[tensor::QgodLocal::size()];
    real QgodNeighborData[tensor::QgodNeighbor::size()];
    auto QgodLocal = init::QgodLocal::view::create(QgodLocalData);
    auto QgodNeighbor = init::QgodNeighbor::view::create(QgodNeighborData);
    
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
      unsigned meshId = ltsToMesh[cell];

      real x[4];
      real y[4];
      real z[4];
      real gradXi[3];
      real gradEta[3];
      real gradZeta[3];

      // Iterate over all 4 vertices of the tetrahedron
      for (unsigned vertex = 0; vertex < 4; ++vertex) {
        VrtxCoords const& coords = vertices[ elements[meshId].vertices[vertex] ].coords;
        x[vertex] = coords[0];
        y[vertex] = coords[1];
        z[vertex] = coords[2];
      }

      seissol::transformations::tetrahedronGlobalToReferenceJacobian( x, y, z, gradXi, gradEta, gradZeta );

      seissol::model::getTransposedCoefficientMatrix( material[cell].local, 0, AT );
      seissol::model::getTransposedCoefficientMatrix( material[cell].local, 1, BT );
      seissol::model::getTransposedCoefficientMatrix( material[cell].local, 2, CT );
      setStarMatrix(ATData, BTData, CTData, gradXi, localIntegration[cell].starMatrices[0]);
      setStarMatrix(ATData, BTData, CTData, gradEta, localIntegration[cell].starMatrices[1]);
      setStarMatrix(ATData, BTData, CTData, gradZeta, localIntegration[cell].starMatrices[2]);

      double volume = MeshTools::volume(elements[meshId], vertices);

      for (unsigned side = 0; side < 4; ++side) {
        seissol::model::getTransposedGodunovState(  material[cell].local,
                                                    material[cell].neighbor[side],
                                                    cellInformation[cell].faceTypes[side],
                                                    QgodLocal,
                                                    QgodNeighbor );

        VrtxCoords normal;
        VrtxCoords tangent1;
        VrtxCoords tangent2;
        MeshTools::normalAndTangents(elements[meshId], side, vertices, normal, tangent1, tangent2);
        double surface = MeshTools::surface(normal);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);

        real NLocalData[6*6];
        seissol::model::getBondMatrix(normal, tangent1, tangent2, NLocalData);
        if (material[cell].local.getMaterialType() == seissol::model::MaterialType::anisotropic) {
          seissol::model::getTransposedGodunovState(  seissol::model::getRotatedMaterialCoefficients(NLocalData, *dynamic_cast<seissol::model::AnisotropicMaterial*>(&material[cell].local)),
                                                      seissol::model::getRotatedMaterialCoefficients(NLocalData, *dynamic_cast<seissol::model::AnisotropicMaterial*>(&material[cell].neighbor[side])),
                                                      cellInformation[cell].faceTypes[side],
                                                      QgodLocal,
                                                      QgodNeighbor );
          seissol::model::getTransposedCoefficientMatrix( seissol::model::getRotatedMaterialCoefficients(NLocalData, *dynamic_cast<seissol::model::AnisotropicMaterial*>(&material[cell].local)), 0, ATtilde );
        } else {
          seissol::model::getTransposedGodunovState(  material[cell].local,
                                                      material[cell].neighbor[side],     
                                                      cellInformation[cell].faceTypes[side],
                                                      QgodLocal,
                                                      QgodNeighbor );
          seissol::model::getTransposedCoefficientMatrix( material[cell].local, 0, ATtilde );
        }

        // Calculate transposed T instead
        seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, T, Tinv);

        // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
        // must be subtracted.
        real fluxScale = -2.0 * surface / (6.0 * volume);

        kernel::computeFluxSolverLocal localKrnl;
        localKrnl.fluxScale = fluxScale;
        localKrnl.AplusT = localIntegration[cell].nApNm1[side];
        localKrnl.QgodLocal = QgodLocalData;
        localKrnl.T = TData;
        localKrnl.Tinv = TinvData;
        localKrnl.star(0) = ATtildeData;
        localKrnl.execute();
        
        kernel::computeFluxSolverNeighbor neighKrnl;
        neighKrnl.fluxScale = fluxScale;
        neighKrnl.AminusT = neighboringIntegration[cell].nAmNm1[side];
        neighKrnl.QgodNeighbor = QgodNeighborData;
        neighKrnl.T = TData;
        neighKrnl.Tinv = TinvData;
        neighKrnl.star(0) = ATtildeData;
        if (cellInformation[cell].faceTypes[side] == FaceType::dirichlet) {
          // Already rotated!
          neighKrnl.Tinv = init::identityT::Values;
        }
        neighKrnl.execute();
      }

      seissol::model::initializeSpecificLocalData(  material[cell].local,
                                                    &localIntegration[cell].specific );

      seissol::model::initializeSpecificNeighborData( material[cell].local,
                                                      &neighboringIntegration[cell].specific );
    }
#ifdef _OPENMP
    }
#endif
    ltsToMesh += it->getNumberOfCells();
  }
}

void surfaceAreaAndVolume(  MeshReader const&      i_meshReader,
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

void seissol::initializers::initializeBoundaryMappings(const MeshReader& i_meshReader,
                                                       const EasiBoundary* easiBoundary,
                                                       LTSTree* io_ltsTree,
                                                       LTS* i_lts,
                                                       Lut* i_ltsLut) {
  std::vector<Element> const& elements = i_meshReader.getElements();
  std::vector<Vertex> const& vertices = i_meshReader.getVertices();

  unsigned* ltsToMesh = i_ltsLut->getLtsToMeshLut(i_lts->material.mask);

  for (LTSTree::leaf_iterator it = io_ltsTree->beginLeaf(LayerMask(Ghost)); it != io_ltsTree->endLeaf(); ++it) {
    auto* cellInformation = it->var(i_lts->cellInformation);
    auto* boundary = it->var(i_lts->boundaryMapping);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
      const auto& element = elements[ltsToMesh[cell]];
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
        real nodesReferenceData[nodal::tensor::nodes2D::Size];
        std::copy_n(nodal::init::nodes2D::Values,
                    nodal::tensor::nodes2D::Size,
                    nodesReferenceData);
        auto nodesReference = nodal::init::nodes2D::view::create(nodesReferenceData);
        auto nodes = boundary[cell][side].nodes;
        assert(nodes != nullptr);
        auto offset = 0;
        for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
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
        real* TData = boundary[cell][side].TData;
        real* TinvData = boundary[cell][side].TinvData;
        assert(TData != nullptr);
        assert(TinvData != nullptr);
        auto T = init::T::view::create(TData);
        auto Tinv = init::Tinv::view::create(TinvData);

        VrtxCoords normal;
        VrtxCoords tangent1;
        VrtxCoords tangent2;
        MeshTools::normalAndTangents(element, side, vertices, normal, tangent1, tangent2);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);
        seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, T, Tinv);

        // Evaluate easi boundary condition matrices if needed
        real* easiBoundaryMap = boundary[cell][side].easiBoundaryMap;
        real* easiBoundaryConstant = boundary[cell][side].easiBoundaryConstant;
        assert(easiBoundaryMap != nullptr);
        assert(easiBoundaryConstant != nullptr);
        if (cellInformation[cell].faceTypes[side] == FaceType::dirichlet) {
          easiBoundary->query(nodes, easiBoundaryMap, easiBoundaryConstant);
        } else {
          // Boundary should not be evaluated
          std::fill_n(easiBoundaryMap,
                      seissol::tensor::easiBoundaryMap::size(),
                      std::numeric_limits<real>::signaling_NaN());
          std::fill_n(easiBoundaryConstant,
                      seissol::tensor::easiBoundaryConstant::size(),
                      std::numeric_limits<real>::signaling_NaN());
        }

      }
    }
    ltsToMesh += it->getNumberOfCells();
  }
}

void seissol::initializers::initializeDynamicRuptureMatrices( MeshReader const&      i_meshReader,
                                                              LTSTree*               io_ltsTree,
                                                              LTS*                   i_lts,
                                                              Lut*                   i_ltsLut,
                                                              LTSTree*               dynRupTree,
                                                              DynamicRupture*        dynRup,
                                                              unsigned*              ltsFaceToMeshFace,
                                                              GlobalData const&      global,
                                                              TimeStepping const&/*    timeStepping*/ )
{
  real TData[tensor::T::size()];
  real TinvData[tensor::Tinv::size()];
  real APlusData[tensor::star::size(0)];
  real AMinusData[tensor::star::size(0)];

  std::vector<Fault> const& fault = i_meshReader.getFault();
  std::vector<Element> const& elements = i_meshReader.getElements();
  CellDRMapping (*drMapping)[4] = io_ltsTree->var(i_lts->drMapping);
  CellMaterialData* material = io_ltsTree->var(i_lts->material);
  real** derivatives = io_ltsTree->var(i_lts->derivatives);
  real* (*faceNeighbors)[4] = io_ltsTree->var(i_lts->faceNeighbors);
  CellLocalInformation* cellInformation = io_ltsTree->var(i_lts->cellInformation);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (LTSTree::leaf_iterator it = dynRupTree->beginLeaf(LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
    real**                                timeDerivativePlus                                        = it->var(dynRup->timeDerivativePlus);
    real**                                timeDerivativeMinus                                       = it->var(dynRup->timeDerivativeMinus);
    real                                (*imposedStatePlus)[tensor::QInterpolated::size()]          = it->var(dynRup->imposedStatePlus);
    real                                (*imposedStateMinus)[tensor::QInterpolated::size()]         = it->var(dynRup->imposedStateMinus);
    DRGodunovData*                        godunovData                                               = it->var(dynRup->godunovData);
    real                                (*fluxSolverPlus)[tensor::fluxSolver::size()]               = it->var(dynRup->fluxSolverPlus);
    real                                (*fluxSolverMinus)[tensor::fluxSolver::size()]              = it->var(dynRup->fluxSolverMinus);
    DRFaceInformation*                    faceInformation                                           = it->var(dynRup->faceInformation);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                            = it->var(dynRup->waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                           = it->var(dynRup->waveSpeedsMinus);

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
      } else {
        /// \todo check if this is correct
        faceInformation[ltsFace].faceRelation = elements[ fault[meshFace].neighborElement ].sideOrientations[ fault[meshFace].neighborSide ] + 1;
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
      real* timeDerivative1 = NULL;
      real* timeDerivative2 = NULL;
      for (unsigned duplicate = 0; duplicate < Lut::MaxDuplicates; ++duplicate) {
        unsigned ltsId = i_ltsLut->ltsId(i_lts->cellInformation.mask, derivativesMeshId, duplicate);
        if (timeDerivative1 == NULL && (cellInformation[ltsId].ltsSetup >> 9)%2 == 1) {
          timeDerivative1 = derivatives[ i_ltsLut->ltsId(i_lts->derivatives.mask, derivativesMeshId, duplicate) ];
        }
        if (timeDerivative2 == NULL && (cellInformation[ltsId].ltsSetup >> derivativesSide)%2 == 1) {
          timeDerivative2 = faceNeighbors[ i_ltsLut->ltsId(i_lts->faceNeighbors.mask, derivativesMeshId, duplicate) ][ derivativesSide ];
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
        unsigned plusLtsId = (fault[meshFace].element >= 0)          ? i_ltsLut->ltsId(i_lts->drMapping.mask, fault[meshFace].element, duplicate) : std::numeric_limits<unsigned>::max();
        unsigned minusLtsId = (fault[meshFace].neighborElement >= 0) ? i_ltsLut->ltsId(i_lts->drMapping.mask, fault[meshFace].neighborElement, duplicate) : std::numeric_limits<unsigned>::max();

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
      auto T = init::T::view::create(TData);
      auto Tinv = init::Tinv::view::create(TinvData);
      seissol::model::getFaceRotationMatrix(fault[meshFace].normal, fault[meshFace].tangent1, fault[meshFace].tangent2, T, Tinv);

      /// Materials
      seissol::model::Material* plusMaterial;
      seissol::model::Material* minusMaterial;
      unsigned plusLtsId = (fault[meshFace].element >= 0)          ? i_ltsLut->ltsId(i_lts->material.mask, fault[meshFace].element) : std::numeric_limits<unsigned>::max();
      unsigned minusLtsId = (fault[meshFace].neighborElement >= 0) ? i_ltsLut->ltsId(i_lts->material.mask, fault[meshFace].neighborElement) : std::numeric_limits<unsigned>::max();

      assert(plusLtsId != std::numeric_limits<unsigned>::max() || minusLtsId != std::numeric_limits<unsigned>::max());

      if (plusLtsId != std::numeric_limits<unsigned>::max()) {
        plusMaterial = &material[plusLtsId].local;
        minusMaterial = &material[plusLtsId].neighbor[ faceInformation[ltsFace].plusSide ];
      } else {
        assert(minusLtsId != std::numeric_limits<unsigned>::max());
        plusMaterial = &material[minusLtsId].neighbor[ faceInformation[ltsFace].minusSide ];
        minusMaterial = &material[minusLtsId].local;
      }

      /// Wave speeds and Coefficient Matrices
      auto APlus = init::star::view<0>::create(APlusData);
      auto AMinus = init::star::view<0>::create(AMinusData);
      
      waveSpeedsPlus[ltsFace].density = plusMaterial->rho;
      waveSpeedsMinus[ltsFace].density = minusMaterial->rho;
      waveSpeedsPlus[ltsFace].pWaveVelocity = plusMaterial->getPWaveSpeed();
      waveSpeedsPlus[ltsFace].sWaveVelocity = plusMaterial->getSWaveSpeed();
      waveSpeedsMinus[ltsFace].pWaveVelocity = minusMaterial->getPWaveSpeed();
      waveSpeedsMinus[ltsFace].sWaveVelocity = minusMaterial->getSWaveSpeed();

      switch (plusMaterial->getMaterialType()) {
        case seissol::model::MaterialType::anisotropic: {
          logError() << "Dynamic Rupture does not work with anisotropy yet.";
          //TODO(SW): Make DR work with anisotropy 
          break;
        }
        case seissol::model::MaterialType::elastic: {
          seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<seissol::model::ElasticMaterial*>(plusMaterial), 0, APlus);
          seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<seissol::model::ElasticMaterial*>(minusMaterial), 0, AMinus);
          break;
        }
        case seissol::model::MaterialType::viscoelastic: {
          seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<seissol::model::ViscoElasticMaterial*>(plusMaterial), 0, APlus);
          seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<seissol::model::ViscoElasticMaterial*>(minusMaterial), 0, AMinus);
          break;
        }
      }

      /// Transpose Tinv
      dynamicRupture::kernel::transposeTinv ttKrnl;
      ttKrnl.Tinv = TinvData;
      ttKrnl.TinvT = godunovData[ltsFace].TinvT;
      ttKrnl.execute();

      double plusSurfaceArea, plusVolume, minusSurfaceArea, minusVolume;
      if (fault[meshFace].element >= 0) {
        surfaceAreaAndVolume( i_meshReader, fault[meshFace].element, fault[meshFace].side, &plusSurfaceArea, &plusVolume );
      } else {
        /// Blow up solution on purpose if used by mistake
        plusSurfaceArea = 1.e99; plusVolume = 1.0;
      }
      if (fault[meshFace].neighborElement >= 0) {
        surfaceAreaAndVolume( i_meshReader, fault[meshFace].neighborElement, fault[meshFace].neighborSide, &minusSurfaceArea, &minusVolume );
      } else {
        /// Blow up solution on purpose if used by mistake
        minusSurfaceArea = 1.e99; minusVolume = 1.0;
      }

      dynamicRupture::kernel::rotateFluxMatrix krnl;
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
}
