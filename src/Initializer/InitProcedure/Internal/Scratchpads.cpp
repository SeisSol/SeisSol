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
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"

#include <algorithm>
#include <cstddef>

namespace seissol::tensor {
struct Iane;
struct Qext;
struct dQext;
struct dQane;
struct Zinv;
} // namespace seissol::tensor

namespace seissol::initializer::internal {

void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage) {
  constexpr size_t TotalDerivativesSize = kernels::Solver::DerivativesSize;
  constexpr size_t NodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto& layer : ltsStorage.leaves(Ghost)) {

    const auto* cellInformation = layer.var<LTS::CellInformation>();

    // the traces of the neighbors that write the traces of their time derivatives, evaluated over
    // the time step of the cell (see NeighIntegrationRecorder)
    std::size_t neighborTraceCounter{0};
    std::size_t nodalDisplacementsCounter{0};
    std::size_t analyticCounter = 0;
    std::size_t numPlasticCells = 0;

    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (cellInformation[cell].faceTypes[face] == FaceType::Regular &&
            cellInformation[cell].ltsSetup.neighborBuffer(face) == BufferType::Derivatives) {
          ++neighborTraceCounter;
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

    // The time integral of every cell in the local integral, which writes the traces from it, and
    // the traces of the neighbors in the neighbor integral; the two do not overlap in time. The
    // time derivatives of every cell, which it writes the traces of where a face takes them.
    const auto integratedDofsSize = std::max(layer.size() * kernels::Solver::IntegralsSize,
                                             neighborTraceCounter * kernels::Solver::FaceTraceSize);
    layer.setEntrySize<LTS::IntegratedDofsScratch>(integratedDofsSize * sizeof(real));
    layer.setEntrySize<LTS::DerivativesScratch>(layer.size() * TotalDerivativesSize * sizeof(real));
    layer.setEntrySize<LTS::NodalAvgDisplacements>(nodalDisplacementsCounter *
                                                   NodalDisplacementsSize * sizeof(real));

    if constexpr (Config::Solver == SolverType::LinearCKAnelastic) {
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

    if constexpr (Config::MaterialType == model::MaterialType::Poroelastic) {
      layer.setEntrySize<LTS::ZinvExtra>(layer.size() * kernels::familySize<tensor::Zinv>() *
                                         sizeof(real));
    }
  }
}

void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage) {
  constexpr size_t TraceSize = kernels::Solver::FaceTraceSize * sizeof(real);
  for (auto& layer : drStorage.leaves()) {
    // both sides of every face
    layer.setEntrySize<DynamicRupture::TracesOnDevice>(2 * TraceSize * layer.size());
  }
}

} // namespace seissol::initializer::internal
