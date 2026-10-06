// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Scratchpads.h"

#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Common.h"
#include "Kernels/ConfigBoundary.h"
#include "Kernels/SolverSelector.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"

#include <algorithm>
#include <cstddef>
#include <unordered_set>

namespace seissol::initializer::internal {

namespace {

/// Sets the sizes of the scratchpads of a layer of cells of the configuration `Cfg`.
template <typename Cfg>
void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Layer& layer) {
  using real = Real<Cfg>;
  constexpr size_t TotalDerivativesSize = kernels::SolverOf<Cfg>::DerivativesSize;
  constexpr size_t NodalDisplacementsSize = tensor::averageNormalDisplacement<Cfg>::size();

  const auto* cellInformation = layer.var<LTS::CellInformation>();

  // look at const pointers (instead of non-const) to make clang-tidy happy
  std::unordered_set<const real*> registry{};
  auto* faceNeighbors = layer.var<LTS::FaceNeighborsDevice>();

  std::size_t derivativesCounter{0};
  std::size_t integratedDofsCounterLocal{0};
  std::size_t integratedDofsCounterNeighbor{0};
  std::size_t nodalDisplacementsCounter{0};
  std::size_t analyticCounter = 0;
  std::size_t numPlasticCells = 0;
  // at most, as each face counts here, but a neighbor (and side) only once in the recorder
  std::size_t configBoundaryBytes = 0;

  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    const bool needsScratchMemForDerivatives =
        !cellInformation[cell].ltsSetup.hasBuffer(BufferType::Derivatives);
    const bool needsScratchMemForStepIntegral =
        !cellInformation[cell].ltsSetup.hasBuffer(BufferType::StepIntegrals);
    if (needsScratchMemForDerivatives) {
      ++derivativesCounter;
    }
    if (needsScratchMemForStepIntegral) {
      ++integratedDofsCounterLocal;
    }

    // include data provided by ghost layers
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {

      const auto* neighborBuffer = static_cast<const real*>(faceNeighbors[cell][face]);
      const auto neighborConfig = cellInformation[cell].neighborConfigIds[face];
      const bool configBoundary = cellInformation[cell].faceTypes[face] == FaceType::Regular &&
                                  neighborBuffer != nullptr && neighborConfig != configIdOf<Cfg>();

      if (configBoundary) {
        // its time integral goes into the scratch of the faces between configurations
        configBoundaryBytes += kernels::ConfigBoundary<Cfg>::scratchBytes(neighborConfig);
      } else if (registry.find(neighborBuffer) == registry.end()) {
        // the time integral of each neighbor counts once

        // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
        if (neighborBuffer != nullptr) {
          if (cellInformation[cell].faceTypes[face] == FaceType::Regular) {

            const bool isNeighbProvidesDerivatives =
                cellInformation[cell].ltsSetup.neighborBuffer(face) == BufferType::Derivatives;
            if (isNeighbProvidesDerivatives) {
              ++integratedDofsCounterNeighbor;
            }
            registry.insert(neighborBuffer);
          }
        }
      }

      if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity) {
        ++nodalDisplacementsCounter;
      }

      if (cellInformation[cell].faceTypes[face] == FaceType::Analytical) {
        ++analyticCounter;
      }

      if (cellInformation[cell].plasticityEnabled) {
        ++numPlasticCells;
      }
    }
  }

  const auto integratedDofsCounter =
      std::max(integratedDofsCounterLocal, integratedDofsCounterNeighbor);

  layer.setEntrySize<LTS::IntegratedDofsScratch>(
      integratedDofsCounter * kernels::SolverOf<Cfg>::IntegralsSize * sizeof(real));
  layer.setEntrySize<LTS::DerivativesScratch>(derivativesCounter * TotalDerivativesSize *
                                              sizeof(real));
  layer.setEntrySize<LTS::NodalAvgDisplacements>(nodalDisplacementsCounter *
                                                 NodalDisplacementsSize * sizeof(real));

  if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
    layer.setEntrySize<LTS::IDofsAneScratch>(layer.size() * kernels::size<tensor::Iane<Cfg>>() *
                                             sizeof(real));
    layer.setEntrySize<LTS::DerivativesExtScratch>(
        layer.size() *
        (kernels::size<tensor::dQext<Cfg>>(1) + kernels::size<tensor::dQext<Cfg>>(2)) *
        sizeof(real));
    layer.setEntrySize<LTS::DerivativesAneScratch>(
        layer.size() *
        (kernels::size<tensor::dQane<Cfg>>(1) + kernels::size<tensor::dQane<Cfg>>(2)) *
        sizeof(real));
    layer.setEntrySize<LTS::DofsExtScratch>(layer.size() * kernels::size<tensor::Qext<Cfg>>() *
                                            sizeof(real));
  }

  layer.setEntrySize<LTS::AnalyticScratch>(analyticCounter * tensor::INodal<Cfg>::size() *
                                           sizeof(real));
  layer.setEntrySize<LTS::ConfigBoundaryScratch>(configBoundaryBytes);

  if (plasticity) {
    layer.setEntrySize<LTS::FlagScratch>(numPlasticCells * sizeof(unsigned));
    layer.setEntrySize<LTS::QStressNodalScratch>(numPlasticCells * tensor::QStressNodal<Cfg>::Size *
                                                 sizeof(real));
  }

  if constexpr (Cfg::MaterialType == model::MaterialType::Poroelastic) {
    layer.setEntrySize<LTS::ZinvExtra>(layer.size() * kernels::familySize<tensor::Zinv<Cfg>>() *
                                       sizeof(real));
  }
}

/// Sets the sizes of the scratchpads of a layer of fault faces of the configuration `Cfg`.
template <typename Cfg>
void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Layer& layer) {
  constexpr size_t IdofsSize = tensor::Q<Cfg>::size() * sizeof(Real<Cfg>);
  const auto layerSize = layer.size();
  layer.setEntrySize<DynamicRupture::IdofsPlusOnDevice>(IdofsSize * layerSize);
  layer.setEntrySize<DynamicRupture::IdofsMinusOnDevice>(IdofsSize * layerSize);
}

} // namespace

void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage) {
  for (auto& layer : ltsStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      deriveRequiredScratchpadMemoryForWp<decltype(cfg)>(plasticity, layer);
    });
  }
}

void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage) {
  for (auto& layer : drStorage.leaves()) {
    dispatchConfig(layer.getIdentifier().config,
                   [&](auto cfg) { deriveRequiredScratchpadMemoryForDr<decltype(cfg)>(layer); });
  }
}

} // namespace seissol::initializer::internal
