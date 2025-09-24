// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_SOLVER_FREESURFACEINTEGRATOR_H_
#define SEISSOL_SRC_SOLVER_FREESURFACEINTEGRATOR_H_

#include <Memory/Descriptor/Surface.h>
#include <Memory/Tree/Layer.h>
#include <memory>

#include "Geometry/MeshReader.h"
#include "Geometry/Refinement/TriangleRefiner.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::solver {
class FreeSurfaceIntegrator {
  public:
  static constexpr std::size_t MaxRefinement = 3;
  static constexpr std::size_t NumComponents = 3;

  private:
  enum class LocationFlag : std::uint8_t {
    Elastic = 0,
    Acoustic = 1,
    FreeSurface = 2,
    FreeSurfaceWithGravity = 3
  };

  real* projectionMatrixMemory{nullptr};
  real* projectionMatrix[4]{};
  real* projectionMatrixFromFace{nullptr};
  std::size_t numberOfSubTriangles{0};
  std::size_t numberOfAlignedSubTriangles{0};

  static constexpr auto PolyDegree = ConvergenceOrder - 1;
  static constexpr auto NumQuadraturePoints = PolyDegree * PolyDegree;
  bool m_enabled{false};

  void initializeProjectionMatrices(unsigned maxRefinementDepth);
  void computeSubTriangleAverages(
      real* projectionMatrixRow,
      const std::array<std::array<double, 3>, NumQuadraturePoints>& bfPoints,
      const double* weights) const;
  void computeSubTriangleAveragesFromFaces(
      real* projectionMatrixFromFaceRow,
      const std::array<std::array<double, 2>, NumQuadraturePoints>& bfPoints,
      const double* weights) const;
  void initializeSurfaceStorage(LTS::Storage& ltsStorage);

  static LocationFlag
      getLocationFlag(CellMaterialData materialData, FaceType faceType, unsigned face);

  public:
  std::array<real*, NumComponents> velocities;
  std::array<real*, NumComponents> displacements;

  std::vector<std::uint8_t> locationFlags;
  std::size_t totalNumberOfFreeSurfaces{0};
  std::size_t totalNumberOfTriangles{0};
  std::vector<std::size_t> backmap;
  std::vector<std::size_t> globalIds;

  SurfaceLTS::Storage* surfaceStorage{nullptr};
  seissol::refinement::TriangleRefiner triRefiner;

  explicit FreeSurfaceIntegrator();
  ~FreeSurfaceIntegrator();

  FreeSurfaceIntegrator(const FreeSurfaceIntegrator&) = delete;
  auto operator=(const FreeSurfaceIntegrator&) -> FreeSurfaceIntegrator& = delete;

  FreeSurfaceIntegrator(FreeSurfaceIntegrator&&) = delete;
  auto operator=(FreeSurfaceIntegrator&&) -> FreeSurfaceIntegrator& = delete;

  void initialize(unsigned maxRefinementDepth,
                  GlobalData* globalData,
                  LTS::Storage& ltsStorage,
                  SurfaceLTS::Storage& surfaceStorage);

  void calculateOutput() const;

  [[nodiscard]] bool enabled() const { return m_enabled; }
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_FREESURFACEINTEGRATOR_H_
