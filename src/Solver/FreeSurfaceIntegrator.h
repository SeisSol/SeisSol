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
#include <Memory/GlobalData.h>
#include <Memory/Tree/Layer.h>
#include <memory>

#include "Geometry/MeshReader.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::solver {
class FreeSurfaceIntegrator {
  public:
  static constexpr std::size_t NumComponents = 3;

  private:
  enum class LocationFlag : std::uint8_t {
    Elastic = 0,
    Acoustic = 1,
    FreeSurface = 2,
    FreeSurfaceWithGravity = 3
  };

  bool m_enabled{false};

  void initializeSurfaceStorage(LTS::Storage& ltsStorage);

  static LocationFlag
      getLocationFlag(CellMaterialData materialData, FaceType faceType, unsigned face);

  public:
  std::size_t totalNumberOfFreeSurfaces{0};
  std::vector<std::size_t> backmap;

  SurfaceLTS::Storage* surfaceStorage{nullptr};

  explicit FreeSurfaceIntegrator();
  ~FreeSurfaceIntegrator();

  FreeSurfaceIntegrator(const FreeSurfaceIntegrator&) = delete;
  auto operator=(const FreeSurfaceIntegrator&) -> FreeSurfaceIntegrator& = delete;

  FreeSurfaceIntegrator(FreeSurfaceIntegrator&&) = delete;
  auto operator=(FreeSurfaceIntegrator&&) -> FreeSurfaceIntegrator& = delete;

  void initialize(unsigned maxRefinementDepth,
                  LTS::Storage& ltsStorage,
                  SurfaceLTS::Storage& surfaceStorage);

  void calculateOutput() const;

  [[nodiscard]] bool enabled() const { return m_enabled; }
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_FREESURFACEINTEGRATOR_H_
