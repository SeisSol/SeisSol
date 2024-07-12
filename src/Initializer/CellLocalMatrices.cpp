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

#include "Initializer/MemoryManager.h"
#include "Initializer/ParameterDB.h"
#include "tree/Layer.hpp"
#include "Numerical_aux/Transformation.h"
#include "Equations/Setup.h"
#include "Model/common.hpp"
#include "Geometry/MeshTools.h"
#include "generated_code/tensor.h"
#include "generated_code/kernel.h"
#include <utils/logger.h>
#ifdef ACL_DEVICE
#include <device.h>
#endif

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

void seissol::initializer::initializeCellLocalMatrices( seissol::geometry::MeshReader const&      i_meshReader,
                                                         LTSTree*               io_ltsTree,
                                                         LTS*                   i_lts,
                                                         Lut*                   i_ltsLut,
                                                         TimeStepping const&    timeStepping )
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
      unsigned clusterId = cellInformation[cell].clusterId;
      auto timeStepWidth = timeStepping.globalCflTimeStepWidths[clusterId];
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
        if (cellInformation[cell].faceTypes[side] == FaceType::dirichlet ||
            cellInformation[cell].faceTypes[side] == FaceType::freeSurfaceGravity) {
          // Already rotated!
          neighKrnl.Tinv = init::identityT::Values;
        }
        neighKrnl.execute();
      }

      seissol::model::initializeSpecificLocalData(  material[cell].local,
                                                    timeStepWidth,
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

void seissol::initializer::initializeBoundaryMappings(const seissol::geometry::MeshReader& i_meshReader,
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

/**
 * Copies an eigen3 matrix to a 2D yateto tensor
 */
template <typename T, int dim1, int dim2>
void copyEigenToYateto(Eigen::Matrix<T, dim1, dim2> const& matrix,
                       yateto::DenseTensorView<2, T>& tensorView) {
  assert(tensorView.shape(0) == dim1);
  assert(tensorView.shape(1) == dim2);

  tensorView.setZero();
  for (size_t row = 0; row < dim1; ++row) {
    for (size_t col = 0; col < dim2; ++col) {
      tensorView(row, col) = matrix(row, col);
    }
  }
}

constexpr int N = tensor::Zminus::Shape[0];
Eigen::Matrix<real, N, N> extractMatrix(eigenvalues::Eigenpair<std::complex<double>, seissol::model::Material_t::NumberOfQuantities> eigenpair) {
#ifdef USE_POROELASTIC
  constexpr std::array<int, 4> tractionIndices = {0,3,5,9};
  constexpr std::array<int, 4> velocityIndices = {6,7,8,10};
  constexpr std::array<int, 4> columnIndices = {0,1,2,3};
#else
  constexpr std::array<int, 3> tractionIndices = {0,3,5};
  constexpr std::array<int, 3> velocityIndices = {6,7,8};
  constexpr std::array<int, 3> columnIndices = {0,1,2};
#endif
  auto matrix = eigenpair.getVectorsAsMatrix();
  Eigen::Matrix<double, N, N> RT = matrix(tractionIndices, columnIndices).real();
  Eigen::Matrix<double, N, N> RT_inv = RT.inverse();
  Eigen::Matrix<double, N, N> RU = matrix(velocityIndices, columnIndices).real();
  Eigen::Matrix<double, N, N> M = RU * RT_inv;
  return M.cast<real>();
};

void seissol::initializer::initializeDynamicRuptureMatrices( seissol::geometry::MeshReader const&      i_meshReader,
                                                              LTSTree*               io_ltsTree,
                                                              LTS*                   i_lts,
                                                              Lut*                   i_ltsLut,
                                                              LTSTree*               dynRupTree,
                                                              DynamicRupture*        dynRup,
                                                              unsigned*              ltsFaceToMeshFace,
                                                              GlobalData const&      global,
                                                              double etaHack )
{
  real TData[tensor::T::size()];
  real TinvData[tensor::Tinv::size()];
  real APlusData[tensor::star::size(0)];
  real AMinusData[tensor::star::size(0)];

  std::vector<Fault> const& fault = i_meshReader.getFault();
  std::vector<Element> const& elements = i_meshReader.getElements();
  CellDRMapping (*drMapping)[4] = io_ltsTree->var(i_lts->drMapping);
  CellDRMapping (*drMappingDevice)[4] = io_ltsTree->var(i_lts->drMappingDevice);
  CellMaterialData* material = io_ltsTree->var(i_lts->material);
  real** derivatives = io_ltsTree->var(i_lts->derivatives);
  real* (*faceNeighbors)[4] = io_ltsTree->var(i_lts->faceNeighbors);
  real** derivativesDevice = io_ltsTree->var(i_lts->derivativesDevice);
  real* (*faceNeighborsDevice)[4] = io_ltsTree->var(i_lts->faceNeighborsDevice);
  CellLocalInformation* cellInformation = io_ltsTree->var(i_lts->cellInformation);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (LTSTree::leaf_iterator it = dynRupTree->beginLeaf(LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
    real**                                timeDerivativePlus                                        = it->var(dynRup->timeDerivativePlus);
    real**                                timeDerivativeMinus                                       = it->var(dynRup->timeDerivativeMinus);
    real**                                timeDerivativePlusDevice                                        = it->var(dynRup->timeDerivativePlusDevice);
    real**                                timeDerivativeMinusDevice                                       = it->var(dynRup->timeDerivativeMinusDevice);
    DRGodunovData*                        godunovData                                               = it->var(dynRup->godunovData);
    real                                (*imposedStatePlus)[tensor::QInterpolated::size()]          = it->var(dynRup->imposedStatePlus, AllocationPlace::Host);
    real                                (*imposedStateMinus)[tensor::QInterpolated::size()]         = it->var(dynRup->imposedStateMinus, AllocationPlace::Host);
    real                                (*fluxSolverPlus)[tensor::fluxSolver::size()]               = it->var(dynRup->fluxSolverPlus, AllocationPlace::Host);
    real                                (*fluxSolverMinus)[tensor::fluxSolver::size()]              = it->var(dynRup->fluxSolverMinus, AllocationPlace::Host);
    real                                (*imposedStatePlusDevice)[tensor::QInterpolated::size()]          = it->var(dynRup->imposedStatePlus, AllocationPlace::Device);
    real                                (*imposedStateMinusDevice)[tensor::QInterpolated::size()]         = it->var(dynRup->imposedStateMinus, AllocationPlace::Device);
    real                                (*fluxSolverPlusDevice)[tensor::fluxSolver::size()]               = it->var(dynRup->fluxSolverPlus, AllocationPlace::Device);
    real                                (*fluxSolverMinusDevice)[tensor::fluxSolver::size()]              = it->var(dynRup->fluxSolverMinus, AllocationPlace::Device);
    DRFaceInformation*                    faceInformation                                           = it->var(dynRup->faceInformation);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                            = it->var(dynRup->waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                           = it->var(dynRup->waveSpeedsMinus);
    seissol::dr::ImpedancesAndEta*        impAndEta                                                 = it->var(dynRup->impAndEta);
    seissol::dr::ImpedanceMatrices*       impedanceMatrices                                         = it->var(dynRup->impedanceMatrices);


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
      real* timeDerivative1 = nullptr;
      real* timeDerivative2 = nullptr;
      real* timeDerivative1Device = nullptr;
      real* timeDerivative2Device = nullptr;
      for (unsigned duplicate = 0; duplicate < Lut::MaxDuplicates; ++duplicate) {
        unsigned ltsId = i_ltsLut->ltsId(i_lts->cellInformation.mask, derivativesMeshId, duplicate);
        if (timeDerivative1 == nullptr && (cellInformation[ltsId].ltsSetup >> 9)%2 == 1) {
          timeDerivative1 = derivatives[ i_ltsLut->ltsId(i_lts->derivatives.mask, derivativesMeshId, duplicate) ];
          timeDerivative1Device = derivativesDevice[ i_ltsLut->ltsId(i_lts->derivatives.mask, derivativesMeshId, duplicate) ];
        }
        if (timeDerivative2 == nullptr && (cellInformation[ltsId].ltsSetup >> derivativesSide)%2 == 1) {
          timeDerivative2 = faceNeighbors[ i_ltsLut->ltsId(i_lts->faceNeighbors.mask, derivativesMeshId, duplicate) ][ derivativesSide ];
          timeDerivative2Device = faceNeighborsDevice[ i_ltsLut->ltsId(i_lts->faceNeighbors.mask, derivativesMeshId, duplicate) ][ derivativesSide ];
        }
      }

      assert(timeDerivative1 != nullptr && timeDerivative2 != nullptr);

      if (fault[meshFace].element >= 0) {
        timeDerivativePlus[ltsFace] = timeDerivative1;
        timeDerivativeMinus[ltsFace] = timeDerivative2;
        timeDerivativePlusDevice[ltsFace] = timeDerivative1Device;
        timeDerivativeMinusDevice[ltsFace] = timeDerivative2Device;
      } else {
        timeDerivativePlus[ltsFace] = timeDerivative2;
        timeDerivativeMinus[ltsFace] = timeDerivative1;
        timeDerivativePlusDevice[ltsFace] = timeDerivative2Device;
        timeDerivativeMinusDevice[ltsFace] = timeDerivative1Device;
      }

      assert(timeDerivativePlus[ltsFace] != nullptr && timeDerivativeMinus[ltsFace] != nullptr);

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
          CellDRMapping& mappingDevice = drMappingDevice[plusLtsId][ faceInformation[ltsFace].plusSide ];
          mappingDevice.side = faceInformation[ltsFace].plusSide;
          mappingDevice.faceRelation = 0;
          mappingDevice.godunov = &imposedStatePlusDevice[ltsFace][0];
          mappingDevice.fluxSolver = &fluxSolverPlusDevice[ltsFace][0];
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
          CellDRMapping& mappingDevice = drMappingDevice[minusLtsId][ faceInformation[ltsFace].minusSide ];
          mappingDevice.side = faceInformation[ltsFace].minusSide;
          mappingDevice.faceRelation = faceInformation[ltsFace].faceRelation;
          mappingDevice.godunov = &imposedStateMinusDevice[ltsFace][0];
          mappingDevice.fluxSolver = &fluxSolverMinusDevice[ltsFace][0];
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

      //calculate Impedances Z and eta
      impAndEta[ltsFace].zp = (waveSpeedsPlus[ltsFace].density * waveSpeedsPlus[ltsFace].pWaveVelocity);
      impAndEta[ltsFace].zpNeig = (waveSpeedsMinus[ltsFace].density * waveSpeedsMinus[ltsFace].pWaveVelocity);
      impAndEta[ltsFace].zs = (waveSpeedsPlus[ltsFace].density * waveSpeedsPlus[ltsFace].sWaveVelocity);
      impAndEta[ltsFace].zsNeig = (waveSpeedsMinus[ltsFace].density * waveSpeedsMinus[ltsFace].sWaveVelocity);

      impAndEta[ltsFace].invZp = 1/impAndEta[ltsFace].zp;
      impAndEta[ltsFace].invZpNeig = 1/impAndEta[ltsFace].zpNeig;
      impAndEta[ltsFace].invZs = 1/impAndEta[ltsFace].zs;
      impAndEta[ltsFace].invZsNeig = 1/impAndEta[ltsFace].zsNeig;

      impAndEta[ltsFace].etaP = etaHack / (1.0 / impAndEta[ltsFace].zp + 1.0 / impAndEta[ltsFace].zpNeig);
      impAndEta[ltsFace].invEtaS = 1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig;
      impAndEta[ltsFace].etaS = 1.0 / (1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig);

      switch (plusMaterial->getMaterialType()) {
        case seissol::model::MaterialType::acoustic: {
          logError() << "Dynamic Rupture does not work with an acoustic material.";
          break;
        }
        case seissol::model::MaterialType::poroelastic: {
          // TODO (SW) Extract this into a function
          seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<seissol::model::PoroElasticMaterial*>(plusMaterial), 0, APlus);
          seissol::model::getTransposedCoefficientMatrix(*dynamic_cast<seissol::model::PoroElasticMaterial*>(minusMaterial), 0, AMinus);

          auto plusEigenpair = seissol::model::getEigenDecomposition(*dynamic_cast<seissol::model::PoroElasticMaterial*>(plusMaterial));
          auto minusEigenpair = seissol::model::getEigenDecomposition(*dynamic_cast<seissol::model::PoroElasticMaterial*>(minusMaterial));

          // The impedance matrices are diagonal in the (visco)elastic case, so we only store
          // the values Zp, Zs. In the poroelastic case, the fluid pressure and normal component
          // of the traction depend on each other, so we need a more complicated matrix structure.
          Eigen::Matrix<real, N, N> impedanceMatrix = extractMatrix(plusEigenpair);
          Eigen::Matrix<real, N, N> impedanceNeigMatrix = extractMatrix(minusEigenpair);
          Eigen::Matrix<real, N, N> etaMatrix = (impedanceMatrix + impedanceNeigMatrix).inverse();

          auto impedanceView = init::Zplus::view::create(impedanceMatrices[ltsFace].impedance);
          auto impedanceNeigView = init::Zminus::view::create(impedanceMatrices[ltsFace].impedanceNeig);
          auto etaView = init::eta::view::create(impedanceMatrices[ltsFace].eta);

          copyEigenToYateto(impedanceMatrix, impedanceView);
          copyEigenToYateto(impedanceNeigMatrix, impedanceNeigView);
          copyEigenToYateto(etaMatrix, etaView);

          break;
        }
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
      /// Traction matrices for "average" traction
      auto tractionPlusMatrix = init::tractionPlusMatrix::view::create(godunovData[ltsFace].tractionPlusMatrix);
      auto tractionMinusMatrix = init::tractionMinusMatrix::view::create(godunovData[ltsFace].tractionMinusMatrix);
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
      dynamicRupture::kernel::transposeTinv ttKrnl;
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

      dynamicRupture::kernel::rotateFluxMatrix krnl;
      krnl.T = TData;

      real                                (*fluxSolverPlusHost)[tensor::fluxSolver::size()]               = it->var(dynRup->fluxSolverPlus);
      real                                (*fluxSolverMinusHost)[tensor::fluxSolver::size()]              = it->var(dynRup->fluxSolverMinus);

      krnl.fluxSolver = fluxSolverPlusHost[ltsFace];
      krnl.fluxScaleDR = -2.0 * plusSurfaceArea / (6.0 * plusVolume);
      krnl.star(0) = APlusData;
      krnl.execute();

      krnl.fluxSolver = fluxSolverMinusHost[ltsFace];
      krnl.fluxScaleDR = 2.0 * minusSurfaceArea / (6.0 * minusVolume);
      krnl.star(0) = AMinusData;
      krnl.execute();
    }

    layerLtsFaceToMeshFace += it->getNumberOfCells();
  }
}

void seissol::initializer::copyCellMatricesToDevice(LTSTree*          ltsTree,
                                                     LTS*              lts,
                                                     LTSTree*          dynRupTree,
                                                     DynamicRupture*   dynRup,
                                                     LTSTree*          boundaryTree,
                                                     Boundary*         boundary) {
#ifdef ACL_DEVICE

#endif // ACL_DEVICE
}
