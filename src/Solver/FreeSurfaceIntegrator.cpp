// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "FreeSurfaceIntegrator.h"

#include "Initializer/MemoryManager.h"
#include "Memory/MemoryAllocator.h"
#include "Numerical/Functions.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "generated_code/kernel.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <Common/Iterator.h>
#include <Geometry/Refinement/TriangleRefiner.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/PreProcessorMacros.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <init.h>
#include <ostream>
#include <tensor.h>
#include <utils/logger.h>
#include <vector>

namespace seissol::solver {

void FreeSurfaceIntegrator::SurfaceLTS::addTo(seissol::initializer::LTSTree& surfaceLtsTree) {
  const seissol::initializer::LayerMask ghostMask(Ghost);
  surfaceLtsTree.add(dofs, ghostMask, 1, initializer::AllocationMode::HostOnly);
  surfaceLtsTree.add(displacementDofs, ghostMask, 1, initializer::AllocationMode::HostOnly);
  surfaceLtsTree.add(side, ghostMask, 1, initializer::AllocationMode::HostOnly);
  surfaceLtsTree.add(meshId, ghostMask, 1, initializer::AllocationMode::HostOnly);
  surfaceLtsTree.add(boundaryMapping, ghostMask, 1, initializer::AllocationMode::HostOnly);
}

FreeSurfaceIntegrator::FreeSurfaceIntegrator() {
  for (auto& face : projectionMatrix) {
    face = nullptr;
  }

  for (unsigned dim = 0; dim < NumComponents; ++dim) {
    velocities[dim] = nullptr;
    displacements[dim] = nullptr;
  }

  surfaceLts.addTo(surfaceLtsTree);
}

FreeSurfaceIntegrator::~FreeSurfaceIntegrator() {
  for (unsigned dim = 0; dim < NumComponents; ++dim) {
    seissol::memory::free(velocities[dim]);
    seissol::memory::free(displacements[dim]);
  }

  seissol::memory::free(projectionMatrixMemory);
  seissol::memory::free(projectionMatrixFromFace);
}

void FreeSurfaceIntegrator::initialize(unsigned maxRefinementDepth,
                                       GlobalData* globalData,
                                       seissol::initializer::LTS* lts,
                                       seissol::initializer::LTSTree* ltsTree) {
  if (maxRefinementDepth > MaxRefinement) {
    logError()
        << "Free surface integrator: Currently more than 3 levels of refinements are unsupported."
        << std::endl;
  }

  m_enabled = true;

  logInfo() << "Initializing free surface integrator.";
  initializeProjectionMatrices(maxRefinementDepth);
  initializeSurfaceLTSTree(lts, ltsTree);
  logInfo() << "Initializing free surface integrator. Done.";
}

void FreeSurfaceIntegrator::calculateOutput() {
  unsigned offset = 0;
  const seissol::initializer::LayerMask ghostMask(Ghost);
  for (auto& surfaceLayer : surfaceLtsTree.leaves(ghostMask)) {
    real** dofs = surfaceLayer.var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer.var(surfaceLts.displacementDofs);
    unsigned* side = surfaceLayer.var(surfaceLts.side);

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for schedule(static) default(none)                                            \
    shared(offset, surfaceLayer, dofs, displacementDofs, side)
#endif // _OPENMP
    for (unsigned face = 0; face < surfaceLayer.size(); ++face) {
      alignas(Alignment) real subTriangleDofs[tensor::subTriangleDofs::size(MaxRefinement)];

      kernel::subTriangleVelocity vkrnl;
      vkrnl.Q = dofs[face];
      vkrnl.selectVelocity = init::selectVelocity::Values;
      vkrnl.subTriangleProjection(triRefiner.maxDepth) = projectionMatrix[side[face]];
      vkrnl.subTriangleDofs(triRefiner.maxDepth) = subTriangleDofs;
      vkrnl.execute(triRefiner.maxDepth);

      auto addOutput = [&](const std::array<real*, NumComponents>& output) {
        for (unsigned component = 0; component < NumComponents; ++component) {
          real* target =
              output[component] + offset + static_cast<size_t>(face * numberOfSubTriangles);
          /// @yateto_todo fix for multiple simulations
          real* source =
              subTriangleDofs + static_cast<size_t>(component * numberOfAlignedSubTriangles);
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
      dkrnl.subTriangleProjectionFromFace(triRefiner.maxDepth) = projectionMatrixFromFace;
      dkrnl.subTriangleDofs(triRefiner.maxDepth) = subTriangleDofs;
      dkrnl.execute(triRefiner.maxDepth);

      addOutput(displacements);
    }
    offset += surfaceLayer.size() * numberOfSubTriangles;
  }
}

void FreeSurfaceIntegrator::initializeProjectionMatrices(unsigned maxRefinementDepth) {
  // Sub triangles
  triRefiner.refine(maxRefinementDepth);

  const auto projectionMatrixNCols =
      tensor::subTriangleProjection::Shape[tensor::subTriangleProjection::index(maxRefinementDepth)]
                                          [1];

  numberOfSubTriangles = triRefiner.subTris.size();
  numberOfAlignedSubTriangles =
      tensor::subTriangleProjection::size(maxRefinementDepth) / projectionMatrixNCols;

  assert(numberOfAlignedSubTriangles * projectionMatrixNCols ==
         tensor::subTriangleProjection::size(maxRefinementDepth));
  assert(numberOfSubTriangles == (1U << (2U * maxRefinementDepth)));

  const auto projectionMatrixNumberOfReals =
      4 * tensor::subTriangleProjection::size(maxRefinementDepth);
  const auto projectionMatrixFromFaceMemoryNumberOfReals =
      tensor::subTriangleProjectionFromFace::size(maxRefinementDepth);

  projectionMatrixMemory =
      seissol::memory::allocTyped<real>(projectionMatrixNumberOfReals, Alignment);
  projectionMatrixFromFace =
      seissol::memory::allocTyped<real>(projectionMatrixFromFaceMemoryNumberOfReals, Alignment);

  std::fill_n(projectionMatrixMemory, 0, projectionMatrixNumberOfReals);
  std::fill_n(projectionMatrixFromFace, 0, projectionMatrixFromFaceMemoryNumberOfReals);

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    projectionMatrix[face] =
        projectionMatrixMemory +
        static_cast<size_t>(face * tensor::subTriangleProjection::size(maxRefinementDepth));
  }

  // Triangle quadrature points and weights
  // TODO: Use the same quadrature rule, which is used for Dynamic Rupture
  const auto [points, weights] = seissol::quadrature::quadrature<2>(PolyDegree);

  auto points3D =
      std::array<std::array<double, 3>, NumQuadraturePoints>{}; // Points for eval of 3D basis
  auto points2D =
      std::array<std::array<double, 2>, NumQuadraturePoints>{}; // Points for eval of 2D basis

  // Compute projection matrices
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (unsigned tri = 0; tri < numberOfSubTriangles; ++tri) {
      for (unsigned qp = 0; qp < NumQuadraturePoints; ++qp) {
        const seissol::refinement::Triangle& subTri = triRefiner.subTris[tri];
        const auto chiTau = std::array<double, 2>{
            points[qp][0] * (subTri.x[1][0] - subTri.x[0][0]) +
                points[qp][1] * (subTri.x[2][0] - subTri.x[0][0]) + subTri.x[0][0],
            points[qp][0] * (subTri.x[1][1] - subTri.x[0][1]) +
                points[qp][1] * (subTri.x[2][1] - subTri.x[0][1]) + subTri.x[0][1]};
        seissol::transformations::chiTau2XiEtaZeta(face, chiTau, points3D[qp]);
        points2D[qp] = chiTau;
      }
      computeSubTriangleAverages(projectionMatrix[face] + tri, points3D, weights.data());
      if (face == 0) {
        computeSubTriangleAveragesFromFaces(
            projectionMatrixFromFace + tri, points2D, weights.data());
      }
    }
  }
}

void FreeSurfaceIntegrator::computeSubTriangleAverages(
    real* projectionMatrixRow,
    const std::array<std::array<double, 3>, NumQuadraturePoints>& bfPoints,
    const double* weights) const {
  unsigned nbf = 0;
  for (unsigned d = 0; d < ConvergenceOrder; ++d) {
    for (unsigned k = 0; k <= d; ++k) {
      for (unsigned j = 0; j <= d - k; ++j) {
        const unsigned i = d - k - j;

        // Compute subtriangle average via quadrature
        double average = 0.0;
        for (unsigned qp = 0; qp < NumQuadraturePoints; ++qp) {
          average +=
              weights[qp] * seissol::functions::TetraDubinerP(
                                {i, j, k}, {bfPoints[qp][0], bfPoints[qp][1], bfPoints[qp][2]});
        }
        // We have a factor J / area. As J = 2*area we have to multiply the average by 2.
        average *= 2.0;

        projectionMatrixRow[static_cast<size_t>(nbf * numberOfAlignedSubTriangles)] = average;

        ++nbf;
      }
    }
  }
}

void FreeSurfaceIntegrator::computeSubTriangleAveragesFromFaces(
    real* projectionMatrixFromFaceRow,
    const std::array<std::array<double, 2>, NumQuadraturePoints>& bfPoints,
    const double* weights) const {
  unsigned nbf = 0;
  for (unsigned d = 0; d < ConvergenceOrder; ++d) {
    for (unsigned j = 0; j <= d; ++j) {
      // Compute subtriangle average via quadrature
      double average = 0.0;
      for (unsigned qp = 0; qp < NumQuadraturePoints; ++qp) {
        average += weights[qp] *
                   seissol::functions::TriDubinerP({d - j, j}, {bfPoints[qp][0], bfPoints[qp][1]});
      }
      // We have a factor J / area. As J = 2*area we have to multiply the average by 2.
      average *= 2.0;

      projectionMatrixFromFaceRow[static_cast<size_t>(nbf * numberOfAlignedSubTriangles)] = average;
      ++nbf;
    }
  }
}

FreeSurfaceIntegrator::LocationFlag FreeSurfaceIntegrator::getLocationFlag(
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
    std::abort(); // logError aborts, but this hides the compiler warning about missing return
                  // statement
  }
}

void FreeSurfaceIntegrator::initializeSurfaceLTSTree(seissol::initializer::LTS* lts,
                                                     seissol::initializer::LTSTree* ltsTree) {
  const seissol::initializer::LayerMask ghostMask(Ghost);

  surfaceLtsTree.setNumberOfTimeClusters(ltsTree->numChildren());
  surfaceLtsTree.fixate();

  totalNumberOfFreeSurfaces = 0;
  for (auto [layer, surfaceLayer] :
       seissol::common::zip(ltsTree->leaves(ghostMask), surfaceLtsTree.leaves(ghostMask))) {
    auto* cellInformation = layer.var(lts->cellInformation);
    auto* secondaryInformation = layer.var(lts->secondaryInformation);
    auto* cellMaterialData = layer.var(lts->material);

    unsigned numberOfFreeSurfaces = 0;
    const auto layerSize = layer.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : numberOfFreeSurfaces)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layerSize; ++cell) {
      if (secondaryInformation[cell].duplicate == 0) {
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurface ||
              cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity ||
              initializer::isAtElasticAcousticInterface(cellMaterialData[cell], face)) {
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

  for (unsigned dim = 0; dim < NumComponents; ++dim) {
    velocities[dim] = seissol::memory::allocTyped<real>(totalNumberOfTriangles, Alignment);
    displacements[dim] = seissol::memory::allocTyped<real>(totalNumberOfTriangles, Alignment);
  }
  locationFlags = std::vector<unsigned int>(totalNumberOfTriangles, 0);

  /// @ yateto_todo
  unsigned surfaceCellOffset = 0; // Counts all surface cells of all layers
  for (auto [layer, surfaceLayer] :
       seissol::common::zip(ltsTree->leaves(ghostMask), surfaceLtsTree.leaves(ghostMask))) {
    auto* cellInformation = layer.var(lts->cellInformation);
    real(*dofs)[tensor::Q::size()] = layer.var(lts->dofs);
    real*(*faceDisplacements)[4] = layer.var(lts->faceDisplacements);
    real** surfaceDofs = surfaceLayer.var(surfaceLts.dofs);
    real** displacementDofs = surfaceLayer.var(surfaceLts.displacementDofs);
    auto* cellMaterialData = layer.var(lts->material);
    auto* surfaceBoundaryMapping = surfaceLayer.var(surfaceLts.boundaryMapping);
    auto* boundaryMapping = layer.var(lts->boundaryMapping);
    auto* secondaryInformation = layer.var(lts->secondaryInformation);

    unsigned* side = surfaceLayer.var(surfaceLts.side);
    unsigned* meshId = surfaceLayer.var(surfaceLts.meshId);
    unsigned surfaceCell = 0;
    for (unsigned cell = 0; cell < layer.size(); ++cell) {
      if (secondaryInformation[cell].duplicate == 0) {
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          if (initializer::requiresDisplacement(
                  cellInformation[cell], cellMaterialData[cell], face)) {
            assert(faceDisplacements[cell][face] != nullptr);

            surfaceDofs[surfaceCell] = dofs[cell];
            displacementDofs[surfaceCell] = faceDisplacements[cell][face];

            side[surfaceCell] = face;
            meshId[surfaceCell] = secondaryInformation[cell].meshId;
            surfaceBoundaryMapping[surfaceCell] = &boundaryMapping[cell][face];

            for (unsigned i = 0; i < numberOfSubTriangles; ++i) {
              locationFlags[surfaceCellOffset * numberOfSubTriangles + i] =
                  static_cast<unsigned int>(getLocationFlag(
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

} // namespace seissol::solver
