// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "FreeSurfaceIntegrator.h"

#include "Memory/MemoryAllocator.h"
#include "Initializer/MemoryManager.h"
#include "Kernels/Common.h"
#include "Kernels/DenseMatrixOps.h"
#include "Numerical/Functions.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "Parallel/MPI.h"
#include "generated_code/kernel.h"
#include <Common/Iterator.h>
#include <utils/logger.h>

void seissol::solver::FreeSurfaceIntegrator::SurfaceLTS::addTo(seissol::initializer::LTSTree& surfaceLtsTree)
{
  seissol::initializer::LayerMask ghostMask(Ghost);
  surfaceLtsTree.addVar(             dofs, ghostMask,                 1,      initializer::AllocationMode::HostOnly );
  surfaceLtsTree.addVar( displacementDofs, ghostMask,                 1,      initializer::AllocationMode::HostOnly );
  surfaceLtsTree.addVar(             side, ghostMask,                 1,      initializer::AllocationMode::HostOnly );
  surfaceLtsTree.addVar(           meshId, ghostMask,                 1,      initializer::AllocationMode::HostOnly );
  surfaceLtsTree.addVar(  boundaryMapping, ghostMask,                 1,      initializer::AllocationMode::HostOnly );
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
                                                          seissol::initializer::LTS* lts,
                                                          seissol::initializer::LTSTree* ltsTree,
                                                          seissol::initializer::Lut* ltsLut  )
{
  if (maxRefinementDepth > FREESURFACE_MAX_REFINEMENT) {
    logError() << "Free surface integrator: Currently more than 3 levels of refinements are unsupported." << std::endl;
  }

  m_enabled = true;

	int const rank = seissol::MPI::mpi.rank();
	logInfo() << "Initializing free surface integrator.";
  initializeProjectionMatrices(maxRefinementDepth);
  initializeSurfaceLTSTree(lts, ltsTree, ltsLut);
	logInfo() << "Initializing free surface integrator. Done.";
}

void seissol::solver::FreeSurfaceIntegrator::calculateOutput()
{
  unsigned offset = 0;
  seissol::initializer::LayerMask ghostMask(Ghost);
  for (auto& surfaceLayer : surfaceLtsTree.leaves(ghostMask)) {
    real** dofs = surfaceLayer.var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer.var(surfaceLts.displacementDofs);
    unsigned* side = surfaceLayer.var(surfaceLts.side);

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
    #pragma omp parallel for schedule(static) default(none) shared(offset, surfaceLayer, dofs, displacementDofs, side)
#endif // _OPENMP
    for (unsigned face = 0; face < surfaceLayer.getNumberOfCells(); ++face) {
      alignas(Alignment) real subTriangleDofs[tensor::subTriangleDofs::size(FREESURFACE_MAX_REFINEMENT)];

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
      dkrnl.faceDisplacement = displacementDofs[face];
      dkrnl.MV2nTo2m = nodal::init::MV2nTo2m::Values;
      dkrnl.subTriangleProjectionFromFace(triRefiner.maxDepth) = projectionMatrixFromFace.get();
      dkrnl.subTriangleDofs(triRefiner.maxDepth) = subTriangleDofs;
      dkrnl.execute(triRefiner.maxDepth);

      addOutput(displacements);
    }
    offset += surfaceLayer.getNumberOfCells() * numberOfSubTriangles;
  }
}


void seissol::solver::FreeSurfaceIntegrator::initializeProjectionMatrices(unsigned maxRefinementDepth)
{
  // Sub triangles
  triRefiner.refine(maxRefinementDepth);

  const auto projectionMatrixNCols = tensor::subTriangleProjection::Shape[tensor::subTriangleProjection::index(maxRefinementDepth)][1];

  numberOfSubTriangles = triRefiner.subTris.size();
  numberOfAlignedSubTriangles = tensor::subTriangleProjection::size(maxRefinementDepth) / projectionMatrixNCols;

  assert(numberOfAlignedSubTriangles * projectionMatrixNCols == tensor::subTriangleProjection::size(maxRefinementDepth));
  assert(numberOfSubTriangles == (1u << (2u*maxRefinementDepth)));

  const auto projectionMatrixNumberOfReals = 4 * tensor::subTriangleProjection::size(maxRefinementDepth);
  const auto projectionMatrixMemorySize = projectionMatrixNumberOfReals * sizeof(real);
  const auto projectionMatrixFromFaceMemoryNumberOfReals = tensor::subTriangleProjectionFromFace::size(maxRefinementDepth);
  const auto projectionMatrixFromFaceMemorySize = projectionMatrixFromFaceMemoryNumberOfReals * sizeof(real);

  projectionMatrixMemory =
      std::unique_ptr<real>(static_cast<real*>(seissol::memory::allocate(projectionMatrixMemorySize, Alignment)));
  projectionMatrixFromFace =
      std::unique_ptr<real>(static_cast<real*>(seissol::memory::allocate(projectionMatrixFromFaceMemorySize, Alignment)));

  std::fill_n(projectionMatrixMemory.get(), 0, projectionMatrixNumberOfReals);
  std::fill_n(projectionMatrixFromFace.get(), 0, projectionMatrixFromFaceMemoryNumberOfReals);

  for (unsigned face = 0; face < 4; ++face) {
    projectionMatrix[face] = projectionMatrixMemory.get() + face * tensor::subTriangleProjection::size(maxRefinementDepth);
  }

  // Triangle quadrature points and weights
  auto points = new double[numQuadraturePoints][2];
  auto weights = new double[numQuadraturePoints];
  // TODO(SW): Use the same quadrature rule, which is used for Dynamic Rupture
  seissol::quadrature::TriangleQuadrature(points, weights, polyDegree);

  auto points3D = std::array<std::array<double, 3>, numQuadraturePoints>{}; // Points for eval of 3D basis
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
        points2D[qp] = chiTau;
      }
      computeSubTriangleAverages(projectionMatrix[face] + tri, points3D, weights);
      if (face == 0) {
        computeSubTriangleAveragesFromFaces(projectionMatrixFromFace.get() + tri, points2D, weights);
      }
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
  for (unsigned d = 0; d < ConvergenceOrder; ++d) {
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
  for (unsigned d = 0; d < ConvergenceOrder; ++d) {
    for (unsigned j = 0; j <= d; ++j) {
      // Compute subtriangle average via quadrature
      double average = 0.0;
      for (unsigned qp = 0; qp < numQuadraturePoints; ++qp) {
        average += weights[qp] * seissol::functions::TriDubinerP({d - j, j}, {bfPoints[qp][0], bfPoints[qp][1]});
      }
      // We have a factor J / area. As J = 2*area we have to multiply the average by 2.
      average *= 2.0;

      projectionMatrixFromFaceRow[nbf * numberOfAlignedSubTriangles] = average;
      ++nbf;
    }
  }
}


seissol::solver::FreeSurfaceIntegrator::LocationFlag seissol::solver::FreeSurfaceIntegrator::getLocationFlag(
    CellMaterialData materialData, FaceType faceType, unsigned int face) {
  if (initializer::isAcousticSideOfElasticAcousticInterface(materialData, face)) {
    return LocationFlag::Acoustic;
  } else if (initializer::isElasticSideOfElasticAcousticInterface(materialData, face)) {
    return LocationFlag::Elastic;
  } else if (faceType == FaceType::FreeSurface) {
    return LocationFlag::FreeSurface;
  } else if (faceType == FaceType::FreeSurfaceGravity) {
    return LocationFlag::FreeSurfaceWithGravity;
  } else {
    logError() << "Internal error in free surface integrator. Called for invalid cell.";
    std::abort(); // logError aborts, but this hides the compiler warning about missing return statement
  }
}

void seissol::solver::FreeSurfaceIntegrator::initializeSurfaceLTSTree(  seissol::initializer::LTS* lts,
                                                                        seissol::initializer::LTSTree* ltsTree,
                                                                        seissol::initializer::Lut* ltsLut )
{
  seissol::initializer::LayerMask ghostMask(Ghost);

  surfaceLtsTree.setNumberOfTimeClusters(ltsTree->numChildren());
  surfaceLtsTree.fixate();

  totalNumberOfFreeSurfaces = 0;
  for (auto [layer, surfaceLayer] : seissol::common::zip(ltsTree->leaves(ghostMask), surfaceLtsTree.leaves(ghostMask))) {
    auto* cellInformation = layer.var(lts->cellInformation);
    auto* secondaryInformation = layer.var(lts->secondaryInformation);
    auto* cellMaterialData = layer.var(lts->material);

    unsigned numberOfFreeSurfaces = 0;
    const auto layerSize = layer.getNumberOfCells();
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+ : numberOfFreeSurfaces)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layerSize; ++cell) {
      if (secondaryInformation[cell].duplicate == 0) {
        for (unsigned face = 0; face < 4; ++face) {
          if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurface
          || cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity
          || initializer::isAtElasticAcousticInterface(cellMaterialData[cell], face)) {
            ++numberOfFreeSurfaces;
          }
        }
      }
    }
    surfaceLayer.setNumberOfCells(numberOfFreeSurfaces);
    totalNumberOfFreeSurfaces += numberOfFreeSurfaces;
  }
  totalNumberOfTriangles = totalNumberOfFreeSurfaces * numberOfSubTriangles;

  surfaceLtsTree.allocateVariables();
  surfaceLtsTree.touchVariables();

  for (unsigned dim = 0; dim < FREESURFACE_NUMBER_OF_COMPONENTS; ++dim) {
    velocities[dim]     = (real*) seissol::memory::allocate(totalNumberOfTriangles * sizeof(real), Alignment);
    displacements[dim]  = (real*) seissol::memory::allocate(totalNumberOfTriangles * sizeof(real), Alignment);
  }
  locationFlags = std::vector<unsigned int>(totalNumberOfTriangles, 0);

  /// @ yateto_todo
  unsigned surfaceCellOffset = 0; // Counts all surface cells of all layers
  for (auto [layer, surfaceLayer] : seissol::common::zip(ltsTree->leaves(ghostMask), surfaceLtsTree.leaves(ghostMask))) {
    auto* cellInformation = layer.var(lts->cellInformation);
    real (*dofs)[tensor::Q::size()] = layer.var(lts->dofs);
    real* (*faceDisplacements)[4] = layer.var(lts->faceDisplacements);
    real** surfaceDofs = surfaceLayer.var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer.var(surfaceLts.displacementDofs);
    auto* cellMaterialData = layer.var(lts->material);
    auto* surfaceBoundaryMapping = surfaceLayer.var(surfaceLts.boundaryMapping);
    auto* boundaryMapping = layer.var(lts->boundaryMapping);
    auto* secondaryInformation = layer.var(lts->secondaryInformation);

    unsigned* side = surfaceLayer.var(surfaceLts.side);
    unsigned* meshId = surfaceLayer.var(surfaceLts.meshId);
    unsigned surfaceCell = 0;
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      if (secondaryInformation[cell].duplicate == 0) {
        for (unsigned face = 0; face < 4; ++face) {
          if (initializer::requiresDisplacement(cellInformation[cell], cellMaterialData[cell], face)) {
            assert(faceDisplacements[cell][face] != nullptr);

            surfaceDofs[surfaceCell]      = dofs[cell];
            displacementDofs[surfaceCell] = faceDisplacements[cell][face];

            side[surfaceCell]             = face;
            meshId[surfaceCell]           = secondaryInformation[cell].meshId;
            surfaceBoundaryMapping[surfaceCell] = &boundaryMapping[cell][face];

            for (unsigned i = 0; i < numberOfSubTriangles; ++i) {
              locationFlags[surfaceCellOffset * numberOfSubTriangles + i] = static_cast<unsigned int>(getLocationFlag(
                  cellMaterialData[cell], cellInformation[cell].faceTypes[face], face));
            }
            ++surfaceCell;
            ++surfaceCellOffset;
          }
        }
      }
    }
  }
}

