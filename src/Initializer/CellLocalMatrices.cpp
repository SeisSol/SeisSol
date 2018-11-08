/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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

#include <Numerical_aux/Transformation.h>
#include <Model/Setup.h>
#include <Model/common.hpp>
#include <Geometry/MeshTools.h>
#include <generated_code/init.h>

void setStarMatrix( real* i_AT,
                    real* i_BT,
                    real* i_CT,
                    real  i_grad[3],
                    real* o_starMatrix )
{
  for (unsigned idx = 0; idx < seissol::model::AstarT::reals; ++idx) {
    o_starMatrix[idx] = i_grad[0] * i_AT[idx];
  }

  for (unsigned idx = 0; idx < seissol::model::BstarT::reals; ++idx) {
    o_starMatrix[idx] += i_grad[1] * i_BT[idx];
  }

  for (unsigned idx = 0; idx < seissol::model::CstarT::reals; ++idx) {
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

  assert(seissol::model::AplusT::rows == seissol::model::AminusT::rows);
  assert(seissol::model::AplusT::cols == seissol::model::AminusT::cols);
  assert(seissol::model::AplusT::rows == seissol::model::AplusT::cols);

  real AT[seissol::model::AstarT::reals];
  real BT[seissol::model::BstarT::reals];
  real CT[seissol::model::CstarT::reals];
  real FlocalData[seissol::model::AplusT::rows * seissol::model::AplusT::cols];
  real FneighborData[seissol::model::AminusT::rows * seissol::model::AminusT::cols];
  real TData[seissol::model::AplusT::rows * seissol::model::AplusT::cols];
  real TinvData[seissol::model::AplusT::rows * seissol::model::AplusT::cols];

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
  #pragma omp parallel for private(AT, BT, CT, FlocalData, FneighborData, TData, TinvData) schedule(static)
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
      setStarMatrix(AT, BT, CT, gradXi, localIntegration[cell].starMatrices[0]);
      setStarMatrix(AT, BT, CT, gradEta, localIntegration[cell].starMatrices[1]);
      setStarMatrix(AT, BT, CT, gradZeta, localIntegration[cell].starMatrices[2]);

      double volume = MeshTools::volume(elements[meshId], vertices);

      for (unsigned side = 0; side < 4; ++side) {
        DenseMatrixView<seissol::model::AplusT::rows, seissol::model::AplusT::cols> Flocal(FlocalData);
        DenseMatrixView<seissol::model::AminusT::rows, seissol::model::AminusT::cols> Fneighbor(FneighborData);
        DenseMatrixView<seissol::model::AplusT::rows, seissol::model::AplusT::cols> T(TData);
        DenseMatrixView<seissol::model::AplusT::cols, seissol::model::AplusT::rows> Tinv(TinvData);

        seissol::model::getTransposedRiemannSolver( material[cell].local,
                                                    material[cell].neighbor[side],
                                                    cellInformation[cell].faceTypes[side],
                                                    //AT,
                                                    Flocal,
                                                    Fneighbor );

        VrtxCoords normal;
        VrtxCoords tangent1;
        VrtxCoords tangent2;
        MeshTools::normalAndTangents(elements[meshId], side, vertices, normal, tangent1, tangent2);
        double surface = MeshTools::surface(normal);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);

        // Calculate transposed T instead
        seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, T, Tinv);

        MatrixView nApNm1(localIntegration[cell].nApNm1[side], seissol::model::AplusT::reals, seissol::model::AplusT::index);
        MatrixView nAmNm1(neighboringIntegration[cell].nAmNm1[side], seissol::model::AminusT::reals, seissol::model::AminusT::index);

        nApNm1.setZero();
        nAmNm1.setZero();

        // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
        // must be subtracted.
        real fluxScale = -2.0 * surface / (6.0 * volume);

        // \todo Generate a kernel for this
        // Calculates  Tinv^T * F * T^T
        for (unsigned j = 0; j < seissol::model::AplusT::cols; ++j) {
          for (unsigned i = 0; i < 9; ++i) {
            for (unsigned k = 0; k < seissol::model::AplusT::rows; ++k) {
              for (unsigned l = 0; l < seissol::model::AplusT::cols; ++l) {
                nApNm1(i, j) += Tinv(k, i) * Flocal(k, l) * T(j, l);
                nAmNm1(i, j) += Tinv(k, i) * Fneighbor(k, l) * T(j, l);
              }
            }
            nApNm1(i, j) *= fluxScale;
            nAmNm1(i, j) *= fluxScale;
          }
        }
      }

      seissol::model::initializeSpecificLocalData(  material[cell].local,
                                                    &localIntegration[cell].specific );

      seissol::model::initializeSpecificNeighborData( material[cell].local,
                                                      material[cell].neighbor,
                                                      &neighboringIntegration[cell].specific );
    }

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
  real TData[seissol::model::AplusT::rows * seissol::model::AplusT::cols];
  real TinvData[seissol::model::AplusT::rows * seissol::model::AplusT::cols];
  real QgodLocalData[9*9];
  real QgodNeighborData[9*9];
  real APlusData[seissol::model::AstarT::reals];
  real AMinusData[seissol::model::AstarT::reals];

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
    real                                (*imposedStatePlus)[seissol::model::godunovState::reals]    = it->var(dynRup->imposedStatePlus);
    real                                (*imposedStateMinus)[seissol::model::godunovState::reals]   = it->var(dynRup->imposedStateMinus);
    DRGodunovData*                        godunovData                                               = it->var(dynRup->godunovData);
    real                                (*fluxSolverPlus)[seissol::model::fluxSolver::reals]        = it->var(dynRup->fluxSolverPlus);
    real                                (*fluxSolverMinus)[seissol::model::fluxSolver::reals]       = it->var(dynRup->fluxSolverMinus);
    DRFaceInformation*                    faceInformation                                           = it->var(dynRup->faceInformation);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                            = it->var(dynRup->waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                           = it->var(dynRup->waveSpeedsMinus);

#ifdef _OPENMP
  #pragma omp parallel for private(TData, TinvData, QgodLocalData, QgodNeighborData, APlusData, AMinusData) schedule(static)
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
          mapping.fluxKernel = 4*faceInformation[ltsFace].plusSide;
          mapping.godunov = &imposedStatePlus[ltsFace][0];
          mapping.fluxSolver = &fluxSolverPlus[ltsFace][0];
          mapping.fluxMatrix = global.nodalFluxMatrices[ faceInformation[ltsFace].plusSide ][0];
}
        }
        if (minusLtsId != std::numeric_limits<unsigned>::max()) {
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
{
          CellDRMapping& mapping = drMapping[minusLtsId][ faceInformation[ltsFace].minusSide ];
          mapping.fluxKernel = 4*faceInformation[ltsFace].minusSide + faceInformation[ltsFace].faceRelation;
          mapping.godunov = &imposedStateMinus[ltsFace][0];
          mapping.fluxSolver = &fluxSolverMinus[ltsFace][0];
          mapping.fluxMatrix = global.nodalFluxMatrices[ faceInformation[ltsFace].minusSide ][ faceInformation[ltsFace].faceRelation ];
}
        }
      }

      /// Transformation matrix
      DenseMatrixView<seissol::model::AplusT::rows, seissol::model::AplusT::cols> T(TData);
      DenseMatrixView<seissol::model::AplusT::cols, seissol::model::AplusT::rows> Tinv(TinvData);
      seissol::model::getFaceRotationMatrix(fault[meshFace].normal, fault[meshFace].tangent1, fault[meshFace].tangent2, T, Tinv);

      /// Materials
      seissol::model::Material plusMaterial;
      seissol::model::Material minusMaterial;
      unsigned plusLtsId = (fault[meshFace].element >= 0)          ? i_ltsLut->ltsId(i_lts->material.mask, fault[meshFace].element) : std::numeric_limits<unsigned>::max();
      unsigned minusLtsId = (fault[meshFace].neighborElement >= 0) ? i_ltsLut->ltsId(i_lts->material.mask, fault[meshFace].neighborElement) : std::numeric_limits<unsigned>::max();

      assert(plusLtsId != std::numeric_limits<unsigned>::max() || minusLtsId != std::numeric_limits<unsigned>::max());

      if (plusLtsId != std::numeric_limits<unsigned>::max()) {
        plusMaterial = material[plusLtsId].local;
        minusMaterial = material[plusLtsId].neighbor[ faceInformation[ltsFace].plusSide ];
      } else {
        assert(minusLtsId != std::numeric_limits<unsigned>::max());
        plusMaterial = material[minusLtsId].neighbor[ faceInformation[ltsFace].minusSide ];
        minusMaterial = material[minusLtsId].local;
      }

      /// Wave speeds
      waveSpeedsPlus[ltsFace].density = plusMaterial.rho;
      waveSpeedsPlus[ltsFace].pWaveVelocity = sqrt( (plusMaterial.lambda + 2.0*plusMaterial.mu) / plusMaterial.rho);
      waveSpeedsPlus[ltsFace].sWaveVelocity = sqrt( plusMaterial.mu / plusMaterial.rho);
      waveSpeedsMinus[ltsFace].density = minusMaterial.rho;
      waveSpeedsMinus[ltsFace].pWaveVelocity = sqrt( (minusMaterial.lambda + 2.0*minusMaterial.mu) / minusMaterial.rho);
      waveSpeedsMinus[ltsFace].sWaveVelocity = sqrt( minusMaterial.mu / minusMaterial.rho);

      /// Godunov state
      DenseMatrixView<9, 9> QgodLocal(QgodLocalData);
      DenseMatrixView<9, 9> QgodNeighbor(QgodNeighborData);
      seissol::model::getTransposedElasticGodunovState( plusMaterial, minusMaterial, QgodLocal, QgodNeighbor );

      // \todo Generate a kernel for this
      MatrixView godunovMatrixPlus(godunovData[ltsFace].godunovMatrixPlus, seissol::model::godunovMatrix::reals, seissol::model::godunovMatrix::index);
      MatrixView godunovMatrixMinus(godunovData[ltsFace].godunovMatrixMinus, seissol::model::godunovMatrix::reals, seissol::model::godunovMatrix::index);
      for (unsigned j = 0; j < 9; ++j) {
        for (unsigned i = 0; i < 9; ++i) {
          for (unsigned k = 0; k < 9; ++k) {
            godunovMatrixPlus(i, j) += Tinv(k, i) * QgodLocal(k, j);
            godunovMatrixMinus(i, j) += Tinv(k, i) * QgodNeighbor(k, j);
          }
        }
      }

      MatrixView APlus(APlusData, seissol::model::AstarT::reals, seissol::model::AstarT::index);
      MatrixView AMinus(AMinusData, seissol::model::AstarT::reals, seissol::model::AstarT::index);
      seissol::model::getTransposedCoefficientMatrix(plusMaterial, 0, APlusData);
      seissol::model::getTransposedCoefficientMatrix(minusMaterial, 0, AMinusData);

      MatrixView fluxSolverPlusView(fluxSolverPlus[ltsFace], seissol::model::fluxSolver::reals, seissol::model::fluxSolver::index);
      MatrixView fluxSolverMinusView(fluxSolverMinus[ltsFace], seissol::model::fluxSolver::reals, seissol::model::fluxSolver::index);

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

      // \todo Generate a kernel for this
      double fluxScalePlus = -2.0 * plusSurfaceArea / (6.0 * plusVolume);
      double fluxScaleMinus = 2.0 * minusSurfaceArea / (6.0 * minusVolume);
      for (unsigned j = 0; j < seissol::model::fluxSolver::cols; ++j) {
        for (unsigned i = 0; i < 9; ++i) {
          for (unsigned k = 0; k < seissol::model::fluxSolver::cols; ++k) {
            fluxSolverPlusView(i, j) += APlus.value(i, k) * T(j, k);
            fluxSolverMinusView(i, j) += AMinus.value(i, k) * T(j, k);
          }
          fluxSolverPlusView(i, j) *= fluxScalePlus;
          fluxSolverMinusView(i, j) *= fluxScaleMinus;
        }
      }

      MatrixView tractionMatrix(godunovData[ltsFace].tractionMatrix, seissol::model::tractionMatrix::reals, seissol::model::tractionMatrix::index);
      tractionMatrix.setZero();
      for (unsigned d = 0; d < 3; ++d) {
        tractionMatrix(d,d) = fault[meshFace].normal[d];
        tractionMatrix(d+3,d) = fault[meshFace].normal[(d+1)%3];
        tractionMatrix(d+3,(d+1)%3) = fault[meshFace].normal[d];
      }
    }

    layerLtsFaceToMeshFace += it->getNumberOfCells();
  }
}
