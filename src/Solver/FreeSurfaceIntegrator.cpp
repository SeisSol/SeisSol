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
#include <Memory/Descriptor/Surface.h>
#include <Memory/Tree/Layer.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <init.h>
#include <limits>
#include <tensor.h>
#include <utils/logger.h>
#include <vector>

namespace seissol::solver {

FreeSurfaceIntegrator::FreeSurfaceIntegrator() {
  for (auto& face : projectionMatrix) {
    face = nullptr;
  }

  for (unsigned dim = 0; dim < NumComponents; ++dim) {
    velocities[dim] = nullptr;
    displacements[dim] = nullptr;
  }
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
                                       LTS::Tree* ltsTree,
                                       SurfaceLTS::Tree* surfaceltsTree) {
  this->surfaceLtsTree = surfaceltsTree;
  if (maxRefinementDepth > MaxRefinement) {
    logError()
        << "Free surface integrator: Currently more than 3 levels of refinements are unsupported.";
  }

  m_enabled = true;

  logInfo() << "Initializing free surface integrator.";
  initializeProjectionMatrices(maxRefinementDepth);
  initializeSurfaceLTSTree(ltsTree);
  logInfo() << "Initializing free surface integrator. Done.";
}

void FreeSurfaceIntegrator::calculateOutput() const {
  const seissol::initializer::LayerMask ghostMask(Ghost);
  for (auto& surfaceLayer : surfaceLtsTree->leaves(ghostMask)) {
    real** dofs = surfaceLayer.var<SurfaceLTS::Dofs>();
    auto* displacementDofs = surfaceLayer.var<SurfaceLTS::DisplacementDofs>();
    auto* side = surfaceLayer.var<SurfaceLTS::Side>();
    auto* outputPosition = surfaceLayer.var<SurfaceLTS::OutputPosition>();

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for schedule(static) default(none)                                            \
    shared(surfaceLayer, dofs, displacementDofs, side, outputPosition)
#endif // _OPENMP
    for (std::size_t face = 0; face < surfaceLayer.size(); ++face) {
      if (outputPosition[face] != std::numeric_limits<std::size_t>::max()) {
        alignas(Alignment) real subTriangleDofs[tensor::subTriangleDofs::size(MaxRefinement)];

        kernel::subTriangleVelocity vkrnl;
        vkrnl.Q = dofs[face];
        vkrnl.selectVelocity = init::selectVelocity::Values;
        vkrnl.subTriangleProjection(triRefiner.maxDepth) = projectionMatrix[side[face]];
        vkrnl.subTriangleDofs(triRefiner.maxDepth) = subTriangleDofs;
        vkrnl.execute(triRefiner.maxDepth);

        auto addOutput = [&](const std::array<real*, NumComponents>& output) {
          for (std::size_t component = 0; component < NumComponents; ++component) {
            real* target = output[component] + outputPosition[face] * numberOfSubTriangles;
            /// @yateto_todo fix for multiple simulations
            real* source =
                subTriangleDofs + static_cast<size_t>(component * numberOfAlignedSubTriangles);
            for (std::size_t subtri = 0; subtri < numberOfSubTriangles; ++subtri) {
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
    }
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
  auto* points = new double[NumQuadraturePoints][2];
  auto* weights = new double[NumQuadraturePoints];
  // TODO(SW): Use the same quadrature rule, which is used for Dynamic Rupture
  seissol::quadrature::TriangleQuadrature(points, weights, PolyDegree);

  auto points3D =
      std::array<std::array<double, 3>, NumQuadraturePoints>{}; // Points for eval of 3D basis
  auto points2D =
      std::array<std::array<double, 2>, NumQuadraturePoints>{}; // Points for eval of 2D basis

  // Compute projection matrices
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (std::size_t tri = 0; tri < numberOfSubTriangles; ++tri) {
      for (std::size_t qp = 0; qp < NumQuadraturePoints; ++qp) {
        const seissol::refinement::Triangle& subTri = triRefiner.subTris[tri];
        const auto chiTau = std::array<double, 2>{
            points[qp][0] * (subTri.x[1][0] - subTri.x[0][0]) +
                points[qp][1] * (subTri.x[2][0] - subTri.x[0][0]) + subTri.x[0][0],
            points[qp][0] * (subTri.x[1][1] - subTri.x[0][1]) +
                points[qp][1] * (subTri.x[2][1] - subTri.x[0][1]) + subTri.x[0][1]};
        seissol::transformations::chiTau2XiEtaZeta(face, chiTau.data(), points3D[qp].data());
        points2D[qp] = chiTau;
      }
      computeSubTriangleAverages(projectionMatrix[face] + tri, points3D, weights);
      if (face == 0) {
        computeSubTriangleAveragesFromFaces(projectionMatrixFromFace + tri, points2D, weights);
      }
    }
  }

  delete[] points;
  delete[] weights;
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

void FreeSurfaceIntegrator::initializeSurfaceLTSTree(LTS::Tree* ltsTree) {
  const seissol::initializer::LayerMask ghostMask(Ghost);

  surfaceLtsTree->setLayerCount(ltsTree->getColorMap());
  surfaceLtsTree->fixate();

  totalNumberOfFreeSurfaces = 0;
  for (auto [layer, surfaceLayer] :
       seissol::common::zip(ltsTree->leaves(ghostMask), surfaceLtsTree->leaves(ghostMask))) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* cellMaterialData = layer.var<LTS::Material>();

    std::size_t numberOfFreeSurfaces = 0;
    std::size_t numberOfOutputFreeSurfaces = 0;
    const auto layerSize = layer.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)                                                          \
    reduction(+ : numberOfFreeSurfaces, numberOfOutputFreeSurfaces)
#endif // _OPENMP
    for (std::size_t cell = 0; cell < layerSize; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (initializer::requiresDisplacement(
                cellInformation[cell], cellMaterialData[cell], face)) {
          ++numberOfFreeSurfaces;

          if (secondaryInformation[cell].duplicate == 0) {
            ++numberOfOutputFreeSurfaces;
          }
        }
      }
    }
    surfaceLayer.setNumberOfCells(numberOfFreeSurfaces);
    totalNumberOfFreeSurfaces += numberOfOutputFreeSurfaces;
  }
  totalNumberOfTriangles = totalNumberOfFreeSurfaces * numberOfSubTriangles;
  backmap.resize(totalNumberOfFreeSurfaces);

  surfaceLtsTree->allocateVariables();
  surfaceLtsTree->touchVariables();

  for (std::size_t dim = 0; dim < NumComponents; ++dim) {
    velocities[dim] = seissol::memory::allocTyped<real>(totalNumberOfTriangles, Alignment);
    displacements[dim] = seissol::memory::allocTyped<real>(totalNumberOfTriangles, Alignment);
  }
  locationFlags = std::vector<std::uint8_t>(totalNumberOfTriangles, 0);
  globalIds.resize(totalNumberOfTriangles);

  // NOTE: we store also for space tree duplicates here
  // thus, we need a non-duplicate lookup table (backmap)

  std::size_t surfaceCellOffset = 0; // Counts all surface cells of all layers
  std::size_t surfaceCellGlobal = 0;
  for (auto [layer, surfaceLayer] :
       seissol::common::zip(ltsTree->leaves(ghostMask), surfaceLtsTree->leaves(ghostMask))) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    real(*dofs)[tensor::Q::size()] = layer.var<LTS::Dofs>();
    real*(*faceDisplacements)[4] = layer.var<LTS::FaceDisplacements>();
    real*(*faceDisplacementsDevice)[4] = layer.var<LTS::FaceDisplacementsDevice>();
    real** surfaceDofs = surfaceLayer.var<SurfaceLTS::Dofs>();
    auto* displacementDofs = surfaceLayer.var<SurfaceLTS::DisplacementDofs>();
    auto* displacementDofsDevice =
        surfaceLayer.var<SurfaceLTS::DisplacementDofs>(initializer::AllocationPlace::Device);
    auto* cellMaterialData = layer.var<LTS::Material>();
    auto* surfaceBoundaryMapping = surfaceLayer.var<SurfaceLTS::BoundaryMapping>();
    auto* boundaryMapping = layer.var<LTS::BoundaryMapping>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* locationFlagLayer = surfaceLayer.var<SurfaceLTS::LocationFlag>();

    auto* side = surfaceLayer.var<SurfaceLTS::Side>();
    auto* meshId = surfaceLayer.var<SurfaceLTS::MeshId>();
    auto* outputPosition = surfaceLayer.var<SurfaceLTS::OutputPosition>();
    std::size_t surfaceCell = 0;
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (initializer::requiresDisplacement(
                cellInformation[cell], cellMaterialData[cell], face)) {
          surfaceDofs[surfaceCell] = dofs[cell];

          // NOTE: assign LTSTree data here
          faceDisplacements[cell][face] = displacementDofs[surfaceCell];
          faceDisplacementsDevice[cell][face] = displacementDofsDevice[surfaceCell];

          side[surfaceCell] = face;
          meshId[surfaceCell] = secondaryInformation[cell].meshId;
          surfaceBoundaryMapping[surfaceCell] = &boundaryMapping[cell][face];
          locationFlagLayer[surfaceCell] = static_cast<std::uint8_t>(
              getLocationFlag(cellMaterialData[cell], cellInformation[cell].faceTypes[face], face));

          const auto globalId = secondaryInformation[cell].globalId * 4 + face;

          if (secondaryInformation[cell].duplicate == 0) {
            for (std::size_t i = 0; i < numberOfSubTriangles; ++i) {
              locationFlags[surfaceCellOffset * numberOfSubTriangles + i] =
                  locationFlagLayer[surfaceCell];
              globalIds[surfaceCellOffset * numberOfSubTriangles + i] = globalId;
            }
            outputPosition[surfaceCell] = surfaceCellOffset;
            backmap[surfaceCellOffset] = surfaceCellGlobal;
            ++surfaceCellOffset;
          } else {
            outputPosition[surfaceCell] = std::numeric_limits<std::size_t>::max();
          }
          ++surfaceCell;
          ++surfaceCellGlobal;
        }
      }
    }
  }
}

} // namespace seissol::solver
