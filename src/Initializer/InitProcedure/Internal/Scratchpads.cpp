// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Scratchpads.h"

#include "Common/Constants.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <unordered_set>

namespace seissol::tensor {
struct Iane;
struct Qext;
struct dQext;
struct dQane;
} // namespace seissol::tensor

namespace seissol::initializer::internal {

void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage) {
  constexpr size_t TotalDerivativesSize = kernels::Solver::DerivativesSize;
  constexpr size_t NodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto& layer : ltsStorage.leaves(Ghost)) {

    const auto* cellInformation = layer.var<LTS::CellInformation>();

    // look at const pointers (instead of non-const) to make clang-tidy happy
    std::unordered_set<const real*> registry{};
    auto* faceNeighbors = layer.var<LTS::FaceNeighborsDevice>();

    std::size_t derivativesCounter{0};
    std::size_t integratedDofsCounter{0};
    std::size_t nodalDisplacementsCounter{0};
    std::size_t analyticCounter = 0;
    std::size_t numPlasticCells = 0;

    std::array<std::size_t, 4> freeSurfacePerFace{};
    std::array<std::size_t, 4> dirichletPerFace{};

    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const bool needsScratchMemForDerivatives = !cellInformation[cell].ltsSetup.hasDerivatives();
      if (needsScratchMemForDerivatives) {
        ++derivativesCounter;
      }
      ++integratedDofsCounter;

      // include data provided by ghost layers
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        const real* neighborBuffer = faceNeighbors[cell][face];

        // check whether a neighbor element idofs has not been counted twice
        if ((registry.find(neighborBuffer) == registry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighborBuffer != nullptr) {
            if (cellInformation[cell].faceTypes[face] == FaceType::Regular) {

              const bool isNeighbProvidesDerivatives =
                  cellInformation[cell].ltsSetup.neighborHasDerivatives(face);
              if (isNeighbProvidesDerivatives) {
                ++integratedDofsCounter;
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

        if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity) {
          ++freeSurfacePerFace[face];
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::Dirichlet) {
          ++dirichletPerFace[face];
        }

        if (cellInformation[cell].plasticityEnabled) {
          ++numPlasticCells;
        }
      }
    }
    const auto freeSurfaceCount =
        *std::max_element(freeSurfacePerFace.begin(), freeSurfacePerFace.end());
    const auto dirichletCountPre =
        *std::max_element(dirichletPerFace.begin(), dirichletPerFace.end());

    // FSG also counts as Dirichlet
    const auto dirichletCount = std::max(dirichletCountPre, freeSurfaceCount);

    layer.setEntrySize<LTS::IntegratedDofsScratch>(integratedDofsCounter *
                                                   kernels::Solver::BuffersSize * sizeof(real));
    layer.setEntrySize<LTS::DerivativesScratch>(derivativesCounter * TotalDerivativesSize *
                                                sizeof(real));
    layer.setEntrySize<LTS::NodalAvgDisplacements>(nodalDisplacementsCounter *
                                                   NodalDisplacementsSize * sizeof(real));

    if constexpr (Config::ViscoMode == ViscoImplementation::AnelasticTensor) {
      layer.setEntrySize<LTS::IDofsAneScratch>(layer.size() * kernels::size<tensor::Iane>() *
                                               sizeof(real));
      layer.setEntrySize<LTS::DerivativesExtScratch>(
          layer.size() * (kernels::size<tensor::dQext>(1) + kernels::size<tensor::dQext>(2)) *
          sizeof(real));
      layer.setEntrySize<LTS::DerivativesAneScratch>(
          layer.size() * (kernels::size<tensor::dQane>(1) + kernels::size<tensor::dQane>(2)) *
          sizeof(real));
      layer.setEntrySize<LTS::DofsExtScratch>(layer.size() * kernels::size<tensor::Qext>() *
                                              sizeof(real));
    }

    layer.setEntrySize<LTS::AnalyticScratch>(analyticCounter * tensor::INodal::size() *
                                             sizeof(real));
    if (plasticity) {
      layer.setEntrySize<LTS::FlagScratch>(numPlasticCells * sizeof(unsigned));
      layer.setEntrySize<LTS::QStressNodalScratch>(numPlasticCells * tensor::QStressNodal::Size *
                                                   sizeof(real));
    }

    layer.setEntrySize<LTS::DofsFaceBoundaryNodalScratch>(sizeof(real) * dirichletCount *
                                                          tensor::INodal::size());

    layer.setEntrySize<LTS::RotateDisplacementToFaceNormalScratch>(
        sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize<LTS::RotateDisplacementToGlobalScratch>(
        sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize<LTS::RotatedFaceDisplacementScratch>(sizeof(real) * freeSurfaceCount *
                                                            init::rotatedFaceDisplacement::Size);
    layer.setEntrySize<LTS::DofsFaceNodalScratch>(sizeof(real) * freeSurfaceCount *
                                                  tensor::INodal::size());
    layer.setEntrySize<LTS::PrevCoefficientsScratch>(sizeof(real) * freeSurfaceCount *
                                                     NodalDisplacementsSize);
  }
}

void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage) {
  constexpr size_t IdofsSize = tensor::Q::size() * sizeof(real);
  for (auto& layer : drStorage.leaves()) {
    const auto layerSize = layer.size();
    layer.setEntrySize<DynamicRupture::IdofsPlusOnDevice>(IdofsSize * layerSize);
    layer.setEntrySize<DynamicRupture::IdofsMinusOnDevice>(IdofsSize * layerSize);
  }
}

} // namespace seissol::initializer::internal
