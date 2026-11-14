// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "BoundaryMappings.h"

#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BoundarySetup.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

template <typename Cfg>
std::array<double, 2> referenceFaceNode(std::size_t node) {
  // The generated nodes2D run over the nodes along the dimension of the basis functions, which
  // comes second where simulations are fused, and leave out the entries outside their bounds,
  // which are zero.
  const auto nodes = nodal::init::nodes2D<Cfg>::view::create(nodal::init::nodes2D<Cfg>::Values);
  std::array<double, 2> point{};
  for (std::size_t d = 0; d < point.size(); ++d) {
    point[d] = multisim::multisimTranspose<Cfg>(
        [&](std::size_t i, std::size_t j) { return nodes.isInRange(i, j) ? nodes(i, j) : 0.0; },
        node,
        d);
  }
  return point;
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg)                                                            \
  template std::array<double, 2> referenceFaceNode<Cfg>(std::size_t);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

namespace {

/// A boundary condition which fused simulations of a copy of the cells take.
struct FusedCondition {
  /// the condition, among the conditions of the run
  std::size_t condition{};
  /// the first simulation of the run that takes it in the copy
  std::size_t simulation{};
  /// the fused simulations of the copy that take it
  std::vector<std::size_t> fused;
};

/// The boundary conditions which the fused simulations of every copy of the cells take, of which
/// `conditionOfSimulation` names the one of every simulation of the run. The fused simulations
/// which only fill up the last copy repeat the first simulations.
std::vector<std::vector<FusedCondition>>
    fusedConditions(const std::vector<std::size_t>& conditionOfSimulation, std::size_t fusedCount) {
  const auto simulationCount = conditionOfSimulation.size();
  std::vector<std::vector<FusedCondition>> conditions((simulationCount + fusedCount - 1) /
                                                      fusedCount);
  for (std::size_t copy = 0; copy < conditions.size(); ++copy) {
    for (std::size_t fused = 0; fused < fusedCount; ++fused) {
      const auto simulation = (copy * fusedCount + fused) % simulationCount;
      const auto condition = conditionOfSimulation[simulation];
      auto known = std::find_if(conditions[copy].begin(),
                                conditions[copy].end(),
                                [&](const auto& taken) { return taken.condition == condition; });
      if (known == conditions[copy].end()) {
        conditions[copy].push_back(FusedCondition{condition, simulation, {fused}});
      } else {
        known->fused.push_back(fused);
      }
    }
  }
  return conditions;
}

/// The face data of the boundary faces of a layer of the configuration `Cfg`.
template <typename Cfg>
void initializeBoundaryMappingsOfLayer(
    LTS::Layer& layer,
    const seissol::geometry::MeshReader& meshReader,
    const std::vector<std::optional<DirichletCondition>>& dirichletConditions,
    const std::vector<std::size_t>& conditionOfSimulation) {
  using real = Real<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();

  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* boundary = layer.var<LTS::BoundaryMapping>(Cfg());
  auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();

  // A Dirichlet face samples each condition of its copy once, for all fused simulations that take
  // it.
  const auto conditionsOfCopy = fusedConditions(conditionOfSimulation, Cfg::NumSimulations);

#pragma omp for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    const auto meshId = secondaryInformation[cell].meshId;
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      if (!boundaryProperties(cellInformation[cell].faceTypes[side]).requiresFaceData) {
        continue;
      }

      const auto face =
          seissol::geometry::AffineFaceTransform::fromMeshCell(meshId, side, meshReader);

      // Compute nodal points in global coordinates for each side.
      auto* nodes = boundary[cell][side].nodes;
      assert(nodes != nullptr);
      auto offset = 0;
      for (std::size_t i = 0; i < nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisDim<Cfg>];
           ++i) {
        // Compute the global coordinates for the nodal points.
        const auto nodeReference = referenceFaceNode<Cfg>(i);
        const auto xyz = face.refToSpace(
            seissol::geometry::FaceTransform::FaceVectorT(nodeReference[0], nodeReference[1]));
        for (std::size_t d = 0; d < Cell::Dim; ++d) {
          nodes[offset++] = xyz(d);
        }
      }

      // Compute map that rotates to normal aligned coordinate system.
      real* matTData = boundary[cell][side].dataT;
      real* matTinvData = boundary[cell][side].dataTinv;
      assert(matTData != nullptr);
      assert(matTinvData != nullptr);
      auto matT = init::T<Cfg>::view::create(matTData);
      auto matTinv = init::Tinv<Cfg>::view::create(matTinvData);

      const auto basis = face.faceAlignedBasis();
      seissol::model::getFaceRotationMatrix<Cfg>(
          basis[0].normalized(), basis[1].normalized(), basis[2].normalized(), matT, matTinv);

      // Evaluate easi boundary condition matrices if needed
      real* dirichletMap = boundary[cell][side].dirichletMap;
      real* dirichletOffset = boundary[cell][side].dirichletOffset;
      assert(dirichletMap != nullptr);
      assert(dirichletOffset != nullptr);
      if (cellInformation[cell].faceTypes[side] == FaceType::Dirichlet) {
        const auto faceBarycenter = face.center();

        // Every fused simulation takes the constant terms of its own condition. The map of the
        // face is one for all of them, and so is the frame it is stated in, which their conditions
        // therefore have to agree on.
        const auto& conditions = conditionsOfCopy[secondaryInformation[cell].copyId];
        real globalMapData[tensor::dirichletMapGlobal<Cfg>::size()]{};
        real globalConstantData[tensor::dirichletOffsetGlobal<Cfg>::size()]{};
        BoundaryFrame frame{};
        for (std::size_t i = 0; i < conditions.size(); ++i) {
          const auto& dirichletCondition = dirichletConditions[conditions[i].condition];
          if (!dirichletCondition.has_value()) {
            logError() << "Dirichlet face found, but no boundary condition definition given for "
                          "the simulation"
                       << conditions[i].simulation;
          }
          real mapData[tensor::dirichletMapGlobal<Cfg>::size()]{};
          const auto conditionFrame = dirichletCondition->template query<Cfg>(
              faceBarycenter.data(), mapData, globalConstantData, conditions[i].fused);
          if (i == 0) {
            frame = conditionFrame;
            std::copy_n(mapData, tensor::dirichletMapGlobal<Cfg>::size(), globalMapData);
          } else if (conditionFrame != frame ||
                     !std::equal(mapData,
                                 mapData + tensor::dirichletMapGlobal<Cfg>::size(),
                                 globalMapData)) {
            logError() << "The simulations" << conditions[0].simulation << "and"
                       << conditions[i].simulation
                       << "are fused into the same cell, whose Dirichlet faces have one map and "
                          "frame for all of them: their boundary conditions can differ in their "
                          "constant terms only, but differ otherwise at the face around"
                       << faceBarycenter[0] << faceBarycenter[1] << faceBarycenter[2];
          }
        }

        if (frame == BoundaryFrame::FaceAligned) {
          std::copy_n(globalMapData, tensor::dirichletMap<Cfg>::size(), dirichletMap);
          std::copy_n(globalConstantData, tensor::dirichletOffset<Cfg>::size(), dirichletOffset);
        } else {
          runtime::kernel::rotateBoundaryCondition rotateKrnl;
          rotateKrnl.dirichletMapGlobal =
              runtime::init::dirichletMapGlobal::view(Variant, globalMapData);
          rotateKrnl.dirichletOffsetGlobal =
              runtime::init::dirichletOffsetGlobal::view(Variant, globalConstantData);
          rotateKrnl.dirichletMap = runtime::init::dirichletMap::view(Variant, dirichletMap);
          rotateKrnl.dirichletOffset =
              runtime::init::dirichletOffset::view(Variant, dirichletOffset);
          rotateKrnl.T = runtime::init::T::view(Variant, matTData);
          rotateKrnl.Tinv = runtime::init::Tinv::view(Variant, matTinvData);
          rotateKrnl.execute(Variant);
        }
      } else {
        // Boundary should not be evaluated
        std::fill_n(dirichletMap,
                    seissol::tensor::dirichletMap<Cfg>::size(),
                    std::numeric_limits<real>::signaling_NaN());
        std::fill_n(dirichletOffset,
                    seissol::tensor::dirichletOffset<Cfg>::size(),
                    std::numeric_limits<real>::signaling_NaN());
      }
    }
  }
}

} // namespace

void initializeBoundaryMappings(
    const seissol::geometry::MeshReader& meshReader,
    const std::vector<std::optional<DirichletCondition>>& dirichletConditions,
    const std::vector<std::size_t>& conditionOfSimulation,
    LTS::Storage& ltsStorage) {
  for (auto& layer : ltsStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      initializeBoundaryMappingsOfLayer<decltype(cfg)>(
          layer, meshReader, dirichletConditions, conditionOfSimulation);
    });
  }
}

} // namespace seissol::initializer
