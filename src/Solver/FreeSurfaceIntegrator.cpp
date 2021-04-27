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
}

seissol::solver::FreeSurfaceIntegrator::FreeSurfaceIntegrator()
  : projectionMatrixMemory(NULL), numberOfSubTriangles(0), numberOfAlignedSubTriangles(0), m_enabled(false), totalNumberOfTriangles(0)
{
  for (unsigned face = 0; face < 4; ++face) {
    projectionMatrix[face] = NULL;
  }

  for (unsigned dim = 0; dim < FREESURFACE_NUMBER_OF_COMPONENTS; ++dim) {
    velocities[dim] = NULL;
    displacements[dim] = NULL;
  }

  surfaceLts.addTo(surfaceLtsTree);
}

seissol::solver::FreeSurfaceIntegrator::~FreeSurfaceIntegrator()
{
  seissol::memory::free(projectionMatrixMemory);

  for (unsigned dim = 0; dim < FREESURFACE_NUMBER_OF_COMPONENTS; ++dim) {
    seissol::memory::free(velocities[dim]);
    seissol::memory::free(displacements[dim]);
  }
}


void seissol::solver::FreeSurfaceIntegrator::initialize(  unsigned maxRefinementDepth,
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
  for ( seissol::initializers::LTSTree::leaf_iterator surfaceLayer = surfaceLtsTree.beginLeaf(ghostMask);
        surfaceLayer != surfaceLtsTree.endLeaf();
        ++surfaceLayer)
  {
    real** dofs             = surfaceLayer->var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer->var(surfaceLts.displacementDofs);
    unsigned* side = surfaceLayer->var(surfaceLts.side);

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

      auto addOutput = [&] (real* output[FREESURFACE_NUMBER_OF_COMPONENTS]) {
        for (unsigned component = 0; component < FREESURFACE_NUMBER_OF_COMPONENTS; ++component) {
          real* target = output[component] + offset + face * numberOfSubTriangles;
          /// @yateto_todo fix for multiple simulations
          real* source = subTriangleDofs + component * numberOfAlignedSubTriangles; 
          for (unsigned subtri = 0; subtri < numberOfSubTriangles; ++subtri) {
            target[subtri] = source[subtri];
            if (!std::isfinite(source[subtri])) {
              logError() << "Detected Inf/NaN in free surface output. Aborting.";
            }
          }

        }
      };

      addOutput(velocities);

      kernel::subTriangleDisplacement dkrnl;
      dkrnl.displacement = displacementDofs[face];
      dkrnl.subTriangleProjection(triRefiner.maxDepth) = projectionMatrix[ side[face] ];
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

  size_t projectionMatrixMemorySize = 4 * tensor::subTriangleProjection::size(maxRefinementDepth) * sizeof(real);
  projectionMatrixMemory = (real*) seissol::memory::allocate(projectionMatrixMemorySize, ALIGNMENT);
  memset(projectionMatrixMemory, 0, projectionMatrixMemorySize);

  for (unsigned face = 0; face < 4; ++face) {
    projectionMatrix[face] = projectionMatrixMemory + face * tensor::subTriangleProjection::size(maxRefinementDepth);
  }

  // Triangle quadrature points and weights
  unsigned polyDegree = CONVERGENCE_ORDER-1;
  unsigned numQuadraturePoints = polyDegree*polyDegree;
  double (*points)[2] = new double[numQuadraturePoints][2];
  double*  weights = new double[numQuadraturePoints];
  seissol::quadrature::TriangleQuadrature(points, weights, polyDegree);

  double (*bfPoints)[3] = new double[numQuadraturePoints][3];
  // Compute projection matrices
  for (unsigned face = 0; face < 4; ++face) {
    for (unsigned tri = 0; tri < numberOfSubTriangles; ++tri) {
      for (unsigned qp = 0; qp < numQuadraturePoints; ++qp) {
        double chiTau[2];
        seissol::refinement::Triangle const& subTri = triRefiner.subTris[tri];
        chiTau[0] = points[qp][0] * (subTri.x[1][0] - subTri.x[0][0]) + points[qp][1] * (subTri.x[2][0] - subTri.x[0][0]) + subTri.x[0][0];
        chiTau[1] = points[qp][0] * (subTri.x[1][1] - subTri.x[0][1]) + points[qp][1] * (subTri.x[2][1] - subTri.x[0][1]) + subTri.x[0][1];
        seissol::transformations::chiTau2XiEtaZeta(face, chiTau, bfPoints[qp]);
      }
      computeSubTriangleAverages(projectionMatrix[face] + tri, bfPoints, weights, numQuadraturePoints);
    }
  }

  delete[] points;
  delete[] weights;
  delete[] bfPoints;
}

void seissol::solver::FreeSurfaceIntegrator::computeSubTriangleAverages(real* projectionMatrixRow, double const (*bfPoints)[3], double const* weights, unsigned numQuadraturePoints)
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
    velocities[dim]     = (real*) seissol::memory::allocate(totalNumberOfTriangles * sizeof(real), ALIGNMENT);
    displacements[dim]  = (real*) seissol::memory::allocate(totalNumberOfTriangles * sizeof(real), ALIGNMENT);
  }

  /// @ yateto_todo
  unsigned* ltsToMesh = ltsLut->getLtsToMeshLut(ghostMask);
  for ( seissol::initializers::LTSTree::leaf_iterator layer = ltsTree->beginLeaf(ghostMask), surfaceLayer = surfaceLtsTree.beginLeaf(ghostMask);
        layer != ltsTree->endLeaf() && surfaceLayer != surfaceLtsTree.endLeaf();
        ++layer, ++surfaceLayer) {
    CellLocalInformation* cellInformation = layer->var(lts->cellInformation);
    real (*dofs)[tensor::Q::size()] = layer->var(lts->dofs);
    real** displacements = layer->var(lts->displacements);
    real** surfaceDofs = surfaceLayer->var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer->var(surfaceLts.displacementDofs);
    CellMaterialData* cellMaterialData = layer->var(lts->material);

    unsigned* side = surfaceLayer->var(surfaceLts.side);
    unsigned* meshId = surfaceLayer->var(surfaceLts.meshId);
    unsigned surfaceCell = 0;
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (cellInformation[cell].faceTypes[face] == FaceType::freeSurface
        || cellInformation[cell].faceTypes[face] == FaceType::freeSurfaceGravity
        || initializers::isAtElasticAcousticInterface(cellMaterialData[cell], face)) {
          assert(displacements[cell] != nullptr);

          surfaceDofs[surfaceCell]      = dofs[cell];
          displacementDofs[surfaceCell] = displacements[cell];
          side[surfaceCell]             = face;
          meshId[surfaceCell]           = ltsToMesh[cell];
          ++surfaceCell;
        }
      }
    }
    ltsToMesh += layer->getNumberOfCells();
  }
}
