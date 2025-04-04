// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_SOLVER_FREESURFACEINTEGRATOR_H_
#define SEISSOL_SRC_SOLVER_FREESURFACEINTEGRATOR_H_

#include <memory>

#include "Geometry/MeshReader.h"
#include "Geometry/Refinement/TriangleRefiner.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Lut.h"

#define FREESURFACE_MAX_REFINEMENT 3
#define FREESURFACE_NUMBER_OF_COMPONENTS 3

namespace seissol::solver {
class FreeSurfaceIntegrator;
} // namespace seissol::solver

class seissol::solver::FreeSurfaceIntegrator {
  private:
  enum class LocationFlag {
    Elastic = 0,
    Acoustic = 1,
    FreeSurface = 2,
    FreeSurfaceWithGravity = 3
  };
  struct SurfaceLTS {
    seissol::initializer::Variable<real*> dofs;
    seissol::initializer::Variable<real*> displacementDofs;
    seissol::initializer::Variable<unsigned> side;
    seissol::initializer::Variable<unsigned> meshId;
    seissol::initializer::Variable<CellBoundaryMapping*> boundaryMapping;

    void addTo(seissol::initializer::LTSTree& surfaceLtsTree);
  };

  real* projectionMatrixMemory{nullptr};
  real* projectionMatrix[4]{};
  real* projectionMatrixFromFace{nullptr};
  unsigned numberOfSubTriangles{0};
  unsigned numberOfAlignedSubTriangles{0};

  static constexpr auto polyDegree = ConvergenceOrder - 1;
  static constexpr auto numQuadraturePoints = polyDegree * polyDegree;
  bool mEnabled{false};

  void initializeProjectionMatrices(unsigned maxRefinementDepth);
  void computeSubTriangleAverages(
      real* projectionMatrixRow,
      const std::array<std::array<double, 3>, numQuadraturePoints>& bfPoints,
      const double* weights) const;
  void computeSubTriangleAveragesFromFaces(
      real* projectionMatrixFromFaceRow,
      const std::array<std::array<double, 2>, numQuadraturePoints>& bfPoints,
      const double* weights) const;
  void initializeSurfaceLTSTree(seissol::initializer::LTS* lts,
                                seissol::initializer::LTSTree* ltsTree,
                                seissol::initializer::Lut* ltsLut);

  static LocationFlag
      getLocationFlag(CellMaterialData materialData, FaceType faceType, unsigned face);

  public:
  real* velocities[FREESURFACE_NUMBER_OF_COMPONENTS];
  real* displacements[FREESURFACE_NUMBER_OF_COMPONENTS];

  std::vector<unsigned int> locationFlags;
  unsigned totalNumberOfFreeSurfaces;
  unsigned totalNumberOfTriangles{0};

  SurfaceLTS surfaceLts;
  seissol::initializer::LTSTree surfaceLtsTree;
  seissol::refinement::TriangleRefiner triRefiner;

  explicit FreeSurfaceIntegrator();
  ~FreeSurfaceIntegrator();

  void initialize(unsigned maxRefinementDepth,
                  GlobalData* globalData,
                  seissol::initializer::LTS* lts,
                  seissol::initializer::LTSTree* ltsTree,
                  seissol::initializer::Lut* ltsLut);

  void calculateOutput();

  [[nodiscard]] bool enabled() const { return mEnabled; }
};

#endif // SEISSOL_SRC_SOLVER_FREESURFACEINTEGRATOR_H_
