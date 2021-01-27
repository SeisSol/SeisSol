/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 */

#include "FreeSurfaceIntegrator.h"

#include <Initializer/MemoryAllocator.h>
#include <Initializer/MemoryManager.h>
#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>
#include <Numerical_aux/Functions.h>
#include <Numerical_aux/Quadrature.h>
#include <Numerical_aux/Transformation.h>
#include <Parallel/MPI.h>
#include <generated_code/kernel.h>
#include <utils/logger.h>

void seissol::solver::FreeSurfaceIntegrator::SurfaceLTS::addTo(seissol::initializers::LTSTree& surfaceLtsTree)
{
  seissol::initializers::LayerMask ghostMask(Ghost);
  surfaceLtsTree.addVar(             dofs, ghostMask,                 1,      seissol::memory::Standard );
  surfaceLtsTree.addVar( displacementDofs, ghostMask,                 1,      seissol::memory::Standard );
  surfaceLtsTree.addVar(             side, ghostMask,                 1,      seissol::memory::Standard );
  surfaceLtsTree.addVar(           meshId, ghostMask,                 1,      seissol::memory::Standard );
  surfaceLtsTree.addVar(  boundaryMapping, ghostMask,                 1,      seissol::memory::Standard );
}

seissol::solver::FreeSurfaceIntegrator::FreeSurfaceIntegrator()
  : projectionMatrixMemory(nullptr), numberOfSubTriangles(0), numberOfAlignedSubTriangles(0), m_enabled(false), totalNumberOfTriangles(0)
{
  for (auto& face : projectionMatrix) {
    face = nullptr;
  }

  for (unsigned dim = 0; dim < FREESURFACE_NUMBER_OF_COMPONENTS; ++dim) {
    velocities[dim] = nullptr;
    displacements[dim] = nullptr;
  }

  surfaceLts.addTo(surfaceLtsTree);
}

seissol::solver::FreeSurfaceIntegrator::~FreeSurfaceIntegrator()
{
  for (unsigned dim = 0; dim < FREESURFACE_NUMBER_OF_COMPONENTS; ++dim) {
    seissol::memory::free(velocities[dim]);
    seissol::memory::free(displacements[dim]);
  }
}


void seissol::solver::FreeSurfaceIntegrator::initialize(  unsigned maxRefinementDepth,
                                                          GlobalData* globalData,
                                                          seissol::initializers::LTS* lts,
                                                          seissol::initializers::LTSTree* ltsTree,
                                                          seissol::initializers::Lut* ltsLut  )
{
  if (maxRefinementDepth > FREESURFACE_MAX_REFINEMENT) {
    logError() << "Free surface integrator: Currently more than 3 levels of refinements are unsupported." << std::endl;
  }

  m_enabled = true;

	int const rank = seissol::MPI::mpi.rank();
	logInfo(rank) << "Initializing free surface integrator.";
  initializeProjectionMatrices(maxRefinementDepth);
  initializeSurfaceLTSTree(lts, ltsTree, ltsLut);
	logInfo(rank) << "Initializing free surface integrator. Done.";
}

void seissol::solver::FreeSurfaceIntegrator::calculateOutput()
{
  unsigned offset = 0;
  seissol::initializers::LayerMask ghostMask(Ghost);
  for (auto surfaceLayer = surfaceLtsTree.beginLeaf(ghostMask);
       surfaceLayer != surfaceLtsTree.endLeaf(); ++surfaceLayer) {
    real** dofs = surfaceLayer->var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer->var(surfaceLts.displacementDofs);
    unsigned* side = surfaceLayer->var(surfaceLts.side);
    auto boundaryMapping = surfaceLayer->var(surfaceLts.boundaryMapping);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif // _OPENMP
    for (unsigned face = 0; face < surfaceLayer->getNumberOfCells(); ++face) {
      real subTriangleDofs[tensor::subTriangleDofs::size(FREESURFACE_MAX_REFINEMENT)] __attribute__((aligned(ALIGNMENT)));

      kernel::subTriangleVelocity vkrnl;
      vkrnl.Q = dofs[face];
      vkrnl.selectVelocity = init::selectVelocity::Values;
      vkrnl.subTriangleProjection(triRefiner.maxDepth) = projectionMatrix[ side[face] ];
      vkrnl.subTriangleDofs(triRefiner.maxDepth) = subTriangleDofs;
      vkrnl.execute(triRefiner.maxDepth);

      auto addOutput = [&] (double* output[FREESURFACE_NUMBER_OF_COMPONENTS]) {
        for (unsigned component = 0; component < FREESURFACE_NUMBER_OF_COMPONENTS; ++component) {
          double* target = output[component] + offset + face * numberOfSubTriangles;
          /// @yateto_todo fix for multiple simulations
          real* source = subTriangleDofs + component * numberOfAlignedSubTriangles; 
          for (unsigned subtri = 0; subtri < numberOfSubTriangles; ++subtri) {
            target[subtri] = source[subtri];
          }
        }
      };

      addOutput(velocities);
      auto T = init::T::view::create(boundaryMapping[face]->TData);
      alignas(ALIGNMENT) real rotateVelocityToGlobalData[init::rotateVelocityToGlobal::Size];
      auto rotateVelocityToGlobal = init::rotateVelocityToGlobal::view::create(rotateVelocityToGlobalData);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          // Extract part that rotates velocity from T
          rotateVelocityToGlobal(i,j) = T(i+6, j+6);
        }
      }

      kernel::subTriangleDisplacement dkrnl = subTriangleDisplacementPrototype;
      dkrnl.faceDisplacement = displacementDofs[face];
      dkrnl.V2nTo2m = nodal::init::V2nTo2m::Values;
      dkrnl.subTriangleProjectionFromFace(triRefiner.maxDepth) = projectionMatrixFromFace[ side[face] ];
      dkrnl.rotateVelocityToGlobal = rotateVelocityToGlobalData;
      dkrnl.subTriangleDofs(triRefiner.maxDepth) = subTriangleDofs;
      dkrnl.execute(triRefiner.maxDepth);

      addOutput(displacements);
    }
    offset += surfaceLayer->getNumberOfCells() * numberOfSubTriangles;
  }
}


void seissol::solver::FreeSurfaceIntegrator::initializeProjectionMatrices(unsigned maxRefinementDepth)
{
  // Sub triangles
  triRefiner.refine(maxRefinementDepth);

  numberOfSubTriangles = triRefiner.subTris.size();
  numberOfAlignedSubTriangles = seissol::kernels::getNumberOfAlignedReals(numberOfSubTriangles);

  assert(numberOfSubTriangles == (1u << (2u*maxRefinementDepth)));

  const auto projectionMatrixMemorySize = 4 * tensor::subTriangleProjection::size(maxRefinementDepth) * sizeof(real);
  const auto projectionMatrixFromFaceMemorySize =
      4 * tensor::subTriangleProjectionFromFace::size(maxRefinementDepth) * sizeof(real);

  projectionMatrixMemory =
      std::unique_ptr<real>(static_cast<real*>(seissol::memory::allocate(projectionMatrixMemorySize, ALIGNMENT)));
  projectionMatrixFromFaceMemory =
      std::unique_ptr<real>(static_cast<real*>(seissol::memory::allocate(projectionMatrixFromFaceMemorySize, ALIGNMENT)));

  std::fill_n(projectionMatrixMemory.get(), 0, projectionMatrixMemorySize);
  std::fill_n(projectionMatrixFromFaceMemory.get(), 0, projectionMatrixFromFaceMemorySize);

  for (unsigned face = 0; face < 4; ++face) {
    projectionMatrix[face] = projectionMatrixMemory.get() + face * tensor::subTriangleProjection::size(maxRefinementDepth);
    projectionMatrixFromFace[face] =
        projectionMatrixFromFaceMemory.get() + face * tensor::subTriangleProjectionFromFace::size(maxRefinementDepth);
  }

  // Triangle quadrature points and weights
  auto points = new double[numQuadraturePoints][2];
  auto weights = new double[numQuadraturePoints];
  seissol::quadrature::TriangleQuadrature(points, weights, polyDegree);

  auto points3D = std::array<std::array<double, 3>, numQuadraturePoints>{}; // Points for eval of 2D basis
  auto points2D = std::array<std::array<double, 2>, numQuadraturePoints>{}; // Points for eval of 2D basis

  // Compute projection matrices
  for (unsigned face = 0; face < 4; ++face) {
    for (unsigned tri = 0; tri < numberOfSubTriangles; ++tri) {
      for (unsigned qp = 0; qp < numQuadraturePoints; ++qp) {
        seissol::refinement::Triangle const& subTri = triRefiner.subTris[tri];
        const auto chiTau = std::array<double, 2>{
            points[qp][0] * (subTri.x[1][0] - subTri.x[0][0]) + points[qp][1] * (subTri.x[2][0] - subTri.x[0][0]) +
            subTri.x[0][0],
            points[qp][0] * (subTri.x[1][1] - subTri.x[0][1]) + points[qp][1] * (subTri.x[2][1] - subTri.x[0][1]) +
            subTri.x[0][1]
        };
        seissol::transformations::chiTau2XiEtaZeta(face, chiTau.data(), points3D[qp].data());
      }
      computeSubTriangleAverages(projectionMatrix[face] + tri, points3D, weights);
      computeSubTriangleAveragesFromFaces(projectionMatrixFromFace[face] + tri, points2D, weights);
    }
  }




  delete[] points;
  delete[] weights;
}

void seissol::solver::FreeSurfaceIntegrator::computeSubTriangleAverages(real* projectionMatrixRow,
                                                                        const std::array<std::array<double, 3>,numQuadraturePoints>& bfPoints,
                                                                        double const* weights) const
{
  unsigned nbf = 0;
  for (unsigned d = 0; d < CONVERGENCE_ORDER; ++d) {
    for (unsigned k = 0; k <= d; ++k) {
      for (unsigned j = 0; j <= d-k; ++j) {
        unsigned i = d-k-j;

        // Compute subtriangle average via quadrature
        double average = 0.0;
        for (unsigned qp = 0; qp < numQuadraturePoints; ++qp) {
          average += weights[qp] * seissol::functions::TetraDubinerP({i, j, k}, {bfPoints[qp][0], bfPoints[qp][1], bfPoints[qp][2]});
        }
        // We have a factor J / area. As J = 2*area we have to multiply the average by 2.
        average *= 2.0;

        projectionMatrixRow[nbf * numberOfAlignedSubTriangles] = average;

        ++nbf;
      }
    }
  }
}

void seissol::solver::FreeSurfaceIntegrator::computeSubTriangleAveragesFromFaces(real* projectionMatrixFromFaceRow,
                                                                                 const std::array<std::array<double, 2>, numQuadraturePoints>& bfPoints,
                                                                                 double const* weights) const {
  unsigned nbf = 0;
  for (unsigned d = 0; d < CONVERGENCE_ORDER; ++d) {
    for (unsigned j = 0; j <= d; ++j) {

      // Compute subtriangle average via quadrature
      double average = 0.0;
      for (unsigned qp = 0; qp < numQuadraturePoints; ++qp) {
        average += weights[qp] * seissol::functions::TriDubinerP({d - j, j}, {bfPoints[qp][0], bfPoints[qp][1]});
      }
      // We have a factor J / area. As J = 2*area we have to multiply the average by 2.
      average *= 4.0; // TODO(Lukas) Why 4?

      projectionMatrixFromFaceRow[nbf * numberOfAlignedSubTriangles] = average;
      ++nbf;
    }
  }
}

void seissol::solver::FreeSurfaceIntegrator::initializeSurfaceLTSTree(  seissol::initializers::LTS* lts,
                                                                        seissol::initializers::LTSTree* ltsTree,
                                                                        seissol::initializers::Lut* ltsLut )
{
  seissol::initializers::LayerMask ghostMask(Ghost);

  surfaceLtsTree.setNumberOfTimeClusters(ltsTree->numChildren());
  surfaceLtsTree.fixate();

  totalNumberOfFreeSurfaces = 0;
  for ( seissol::initializers::LTSTree::leaf_iterator layer = ltsTree->beginLeaf(ghostMask), surfaceLayer = surfaceLtsTree.beginLeaf(ghostMask);
        layer != ltsTree->endLeaf() && surfaceLayer != surfaceLtsTree.endLeaf();
        ++layer, ++surfaceLayer) {
    CellLocalInformation* cellInformation = layer->var(lts->cellInformation);
    CellMaterialData* cellMaterialData = layer->var(lts->material);

    unsigned numberOfFreeSurfaces = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+ : numberOfFreeSurfaces)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
        for (unsigned face = 0; face < 4; ++face) {
        if (cellInformation[cell].faceTypes[face] == FaceType::freeSurface
        || cellInformation[cell].faceTypes[face] == FaceType::freeSurfaceGravity
        || initializers::isAtElasticAcousticInterface(cellMaterialData[cell], face)) {
          ++numberOfFreeSurfaces;
        }
      }
    }
    surfaceLayer->setNumberOfCells(numberOfFreeSurfaces);
    totalNumberOfFreeSurfaces += numberOfFreeSurfaces;
  }
  totalNumberOfTriangles = totalNumberOfFreeSurfaces * numberOfSubTriangles;

  surfaceLtsTree.allocateVariables();
  surfaceLtsTree.touchVariables();

  for (unsigned dim = 0; dim < FREESURFACE_NUMBER_OF_COMPONENTS; ++dim) {
    velocities[dim]     = (double*) seissol::memory::allocate(totalNumberOfTriangles * sizeof(double), ALIGNMENT);
    displacements[dim]  = (double*) seissol::memory::allocate(totalNumberOfTriangles * sizeof(double), ALIGNMENT);
  }

  /// @ yateto_todo
  unsigned* ltsToMesh = ltsLut->getLtsToMeshLut(ghostMask);
  for ( seissol::initializers::LTSTree::leaf_iterator layer = ltsTree->beginLeaf(ghostMask), surfaceLayer = surfaceLtsTree.beginLeaf(ghostMask);
        layer != ltsTree->endLeaf() && surfaceLayer != surfaceLtsTree.endLeaf();
        ++layer, ++surfaceLayer) {
    CellLocalInformation* cellInformation = layer->var(lts->cellInformation);
    real (*dofs)[tensor::Q::size()] = layer->var(lts->dofs);
    real* (*displacements)[4] = layer->var(lts->faceDisplacements);
    real** surfaceDofs = surfaceLayer->var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer->var(surfaceLts.displacementDofs);
    CellMaterialData* cellMaterialData = layer->var(lts->material);
    auto surfaceBoundaryMapping = surfaceLayer->var(surfaceLts.boundaryMapping);
    auto boundaryMapping = layer->var(lts->boundaryMapping);

    unsigned* side = surfaceLayer->var(surfaceLts.side);
    unsigned* meshId = surfaceLayer->var(surfaceLts.meshId);
    unsigned surfaceCell = 0;
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (initializers::requiresDisplacement(cellInformation[cell], cellMaterialData[cell], face)) {
          assert(displacements[cell][face] != nullptr);

          surfaceDofs[surfaceCell]      = dofs[cell];
          displacementDofs[surfaceCell] = displacements[cell][face];
          side[surfaceCell]             = face;
          meshId[surfaceCell]           = ltsToMesh[cell];
          surfaceBoundaryMapping[surfaceCell] = &boundaryMapping[cell][face];

          ++surfaceCell;
        }
      }
    }
    ltsToMesh += layer->getNumberOfCells();
  }
}
