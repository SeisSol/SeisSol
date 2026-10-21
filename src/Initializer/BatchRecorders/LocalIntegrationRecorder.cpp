// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "Common/Typedefs.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/Condition.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Model/OperatorLayout.h"
#include "Recorders.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace seissol::initializer;
using namespace seissol::recording;

namespace seissol::tensor {
struct Qane;
struct Qext;
struct dQane;
struct dQext;
struct Zinv;
} // namespace seissol::tensor

// NOLINTBEGIN (-misc-const-correctness)

namespace seissol::recording {

real* liftedFieldOf(LTS::Layer& layer, std::size_t cell) {
  if constexpr (Config::Solver == SolverType::LinearCKAnelastic) {
    auto* dofsExt = layer.var<LTS::DofsExtScratch>(AllocationPlace::Device);
    return static_cast<real*>(dofsExt) + kernels::size<tensor::Qext>() * cell;
  } else {
    return static_cast<real*>(layer.cellRef(cell, AllocationPlace::Device).get<LTS::Dofs>());
  }
}

real* fluxTargetOf(LTS::Layer& layer, std::size_t cell) {
  if (WeightAdjustable &&
      layer.cellRef(cell, AllocationPlace::Host).get<LTS::CellInformation>().weightAdjusted) {
    return static_cast<real*>(layer.cellRef(cell, AllocationPlace::Device).get<LTS::FluxBuffer>());
  }
  return liftedFieldOf(layer, cell);
}

} // namespace seissol::recording

void LocalIntegrationRecorder::record(LTS::Layer& layer) {
  setUpContext(layer);

  recordTimeAndVolumeIntegrals();
  recordFreeSurfaceGravityBc();
  recordDirichletBc();
  recordAnalyticalBc(layer);
  recordFaceTraces();
  recordDisplacements();
}

void LocalIntegrationRecorder::recordTimeAndVolumeIntegrals() {
  real* integratedDofsScratch =
      static_cast<real*>(currentLayer_->var<LTS::IntegratedDofsScratch>(AllocationPlace::Device));
  real* derivativesScratch =
      static_cast<real*>(currentLayer_->var<LTS::DerivativesScratch>(AllocationPlace::Device));

  const auto size = currentLayer_->size();
  if (size > 0) {
    std::vector<real*> dofsPtrs(size, nullptr);
    std::vector<real*> integralsPtrs(size, nullptr);
    std::vector<real*> dofsAnePtrs(size, nullptr);
    std::vector<real*> dofsExtPtrs(size, nullptr);
    std::vector<real*> localPtrs(size, nullptr);
    std::vector<real*> idofsAnePtrs(size, nullptr);
    std::vector<real*> derivativesAnePtrs(size, nullptr);
    std::vector<real*> derivativesExtPtrs(size, nullptr);
    std::vector<real*> zinvExtraPtrs(size, nullptr);
    // the weight-adjusted cells: their flux buffers, the field the weight comes back to, and the
    // material at their samples
    std::vector<real*> fluxBufferPtrs;
    std::vector<real*> weightTargetPtrs;
    std::vector<real*> weightLocalPtrs;

    idofsPtrs_.resize(size);
    dQPtrs_.resize(size);
    fluxTargets_.resize(size);

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      // dofs
      dofsPtrs[cell] = static_cast<real*>(data.get<LTS::Dofs>());
      integralsPtrs[cell] = static_cast<real*>(data.get<LTS::Integrals>());

      // the time integral and the time derivatives of a cell stay in the scratch memory of the
      // layer: what the cell leaves on its faces, it writes from there (recordFaceTraces)
      idofsPtrs_[cell] = integratedDofsScratch + cell * kernels::Solver::IntegralsSize;
      dQPtrs_[cell] = derivativesScratch + cell * kernels::Solver::DerivativesSize;

      // stars
      localPtrs[cell] = reinterpret_cast<real*>(&data.get<LTS::LocalIntegration>());
      if constexpr (Config::Solver == SolverType::LinearCKAnelastic) {
        auto* dofsAne = currentLayer_->var<LTS::DofsAne>(AllocationPlace::Device);
        dofsAnePtrs[cell] = dofsAne[cell];

        auto* idofsAne = currentLayer_->var<LTS::IDofsAneScratch>(AllocationPlace::Device);
        idofsAnePtrs[cell] = static_cast<real*>(idofsAne) + kernels::size<tensor::Iane>() * cell;

        auto* derivativesExt =
            currentLayer_->var<LTS::DerivativesExtScratch>(AllocationPlace::Device);
        derivativesExtPtrs[cell] =
            static_cast<real*>(derivativesExt) +
            (kernels::size<tensor::dQext>(1) + kernels::size<tensor::dQext>(2)) * cell;

        auto* derivativesAne =
            currentLayer_->var<LTS::DerivativesAneScratch>(AllocationPlace::Device);
        derivativesAnePtrs[cell] =
            static_cast<real*>(derivativesAne) +
            (kernels::size<tensor::dQane>(1) + kernels::size<tensor::dQane>(2)) * cell;

        auto* dofsExt = currentLayer_->var<LTS::DofsExtScratch>(AllocationPlace::Device);
        dofsExtPtrs[cell] = static_cast<real*>(dofsExt) + kernels::size<tensor::Qext>() * cell;
      }
      if constexpr (Config::MaterialType == model::MaterialType::Poroelastic) {
        auto* zinvExtraPtr = currentLayer_->var<LTS::ZinvExtra>(AllocationPlace::Device);
        zinvExtraPtrs[cell] = zinvExtraPtr + kernels::familySize<tensor::Zinv>() * cell;
      }

      fluxTargets_[cell] = fluxTargetOf(*currentLayer_, cell);
      if (WeightAdjustable && dataHost.get<LTS::CellInformation>().weightAdjusted) {
        fluxBufferPtrs.push_back(fluxTargets_[cell]);
        weightTargetPtrs.push_back(liftedFieldOf(*currentLayer_, cell));
        weightLocalPtrs.push_back(localPtrs[cell]);
      }
    }

    const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    checkKey(key);

    (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::Integrals, integralsPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::LocalIntegrationData, localPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs_);
    (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, dQPtrs_);
    (*currentTable_)[key].set(inner_keys::Wp::Id::FluxTarget, fluxTargets_);

    if (!fluxBufferPtrs.empty()) {
      // the local integral clears the buffers, the neighbour integral applies the weight to them
      const ConditionalKey weightKey(*KernelNames::WeightAdjust);
      checkKey(weightKey);
      (*currentTable_)[weightKey].set(inner_keys::Wp::Id::FluxBuffer, fluxBufferPtrs);
      (*currentTable_)[weightKey].set(Config::Solver == SolverType::LinearCKAnelastic
                                          ? inner_keys::Wp::Id::DofsExt
                                          : inner_keys::Wp::Id::Dofs,
                                      weightTargetPtrs);
      (*currentTable_)[weightKey].set(inner_keys::Wp::Id::LocalIntegrationData, weightLocalPtrs);
    }

    if constexpr (Config::Solver == SolverType::LinearCKAnelastic) {
      (*currentTable_)[key].set(inner_keys::Wp::Id::DofsAne, dofsAnePtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::DofsExt, dofsExtPtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::IdofsAne, idofsAnePtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::DerivativesAne, derivativesAnePtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::DerivativesExt, derivativesExtPtrs);
    }
    if constexpr (Config::MaterialType == model::MaterialType::Poroelastic) {
      (*currentTable_)[key].set(inner_keys::Wp::Id::ZinvExtra, zinvExtraPtrs);
    }
  }
}

void LocalIntegrationRecorder::recordFaceTraces() {
  auto* faceTraces = currentLayer_->var<LTS::FaceTracesDevice>();
  auto* faceBuffers = currentLayer_->var<LTS::FaceBuffersDevice>();
  const auto* cellInformation = currentLayer_->var<LTS::CellInformation>();

  // by the face and by the side of the face the cell is on (CellLocalInformation::traceSides)
  using PerFaceSide = std::array<std::array<std::vector<real*>, 2>, Cell::NumFaces>;
  PerFaceSide integralPtrs{};
  PerFaceSide tracePtrs{};
  PerFaceSide derivativesPtrs{};
  PerFaceSide faceDerivativesPtrs{};
  // a trace accumulated over the steps of a neighbor that takes fewer
  std::vector<real*> accumulatedTracePtrs;
  std::vector<real*> accumulatedBufferPtrs;

  const auto size = currentLayer_->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto side = cellInformation[cell].traceSides[face];
      integralPtrs[face][side].push_back(idofsPtrs_[cell]);
      tracePtrs[face][side].push_back(faceTraces[cell][face]);

      const auto type = cellInformation[cell].ltsSetup.faceBuffer(face);
      if (type == BufferType::Derivatives) {
        derivativesPtrs[face][side].push_back(dQPtrs_[cell]);
        faceDerivativesPtrs[face][side].push_back(faceBuffers[cell][face]);
      } else if (type == BufferType::AccumulatedIntegrals) {
        accumulatedTracePtrs.push_back(faceTraces[cell][face]);
        accumulatedBufferPtrs.push_back(faceBuffers[cell][face]);
      }
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (std::size_t side = 0; side < 2; ++side) {
      if (!tracePtrs[face][side].empty()) {
        const ConditionalKey key(*KernelNames::FaceTrace, *ComputationKind::None, face, side);
        checkKey(key);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, integralPtrs[face][side]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FaceTrace, tracePtrs[face][side]);
      }
      if (!faceDerivativesPtrs[face][side].empty()) {
        const ConditionalKey key(
            *KernelNames::FaceTrace, *ComputationKind::WithDerivatives, face, side);
        checkKey(key);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, derivativesPtrs[face][side]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FaceBuffer, faceDerivativesPtrs[face][side]);
      }
    }
  }

  if (!accumulatedBufferPtrs.empty()) {
    const ConditionalKey key(*KernelNames::FaceTrace, *ComputationKind::WithLtsBuffers);
    checkKey(key);
    (*currentTable_)[key].set(inner_keys::Wp::Id::FaceTrace, accumulatedTracePtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::FaceBuffer, accumulatedBufferPtrs);
  }
}

void LocalIntegrationRecorder::recordDisplacements() {
  auto* faceDisplacements = currentLayer_->var<LTS::FaceDisplacementsDevice>();
  std::array<std::vector<real*>, Cell::NumFaces> iVelocitiesPtrs{{}};
  std::array<std::vector<real*>, Cell::NumFaces> displacementsPtrs{};

  const auto size = currentLayer_->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      auto isRequired = faceDisplacements[cell][face] != nullptr;
      auto notFreeSurfaceGravity =
          dataHost.get<LTS::CellInformation>().faceTypes[face] != FaceType::FreeSurfaceGravity;

      if (isRequired && notFreeSurfaceGravity) {
        auto iview = init::I::view::create(idofsPtrs_[cell]);
        // gpu_addVelocity reads the three integrated velocities as consecutive columns of I,
        // starting at this pointer. Where they start depends on the equation (column 6 for the
        // elastic ones, column 1 for the acoustic ones); a wrong column reads other quantities or,
        // past the last one, the integrals of other cells, without any error at runtime.
        static_assert(model::MaterialT::VelocityOffset + Cell::Dim <=
                      tensor::I::Shape[multisim::BasisFunctionDimension + 1]);
        iVelocitiesPtrs[face].push_back(
            &multisim::multisimWrap(iview, 0, 0, model::MaterialT::VelocityOffset));
        displacementsPtrs[face].push_back(faceDisplacements[cell][face]);
      }
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (!displacementsPtrs[face].empty()) {
      const ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
      checkKey(key);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Ivelocities, iVelocitiesPtrs[face]);
      (*currentTable_)[key].set(inner_keys::Wp::Id::FaceDisplacement, displacementsPtrs[face]);
    }
  }
}

void LocalIntegrationRecorder::recordFreeSurfaceGravityBc() {
  const auto size = currentLayer_->size();
  constexpr size_t NodalAvgDisplacementsSize = tensor::averageNormalDisplacement::size();

  real* nodalAvgDisplacements =
      static_cast<real*>(currentLayer_->var<LTS::NodalAvgDisplacements>(AllocationPlace::Device));

  if (size > 0) {
    std::array<std::vector<unsigned>, Cell::NumFaces> cellIndices{};
    std::array<std::vector<real*>, Cell::NumFaces> nodalAvgDisplacementsPtrs{};
    std::array<std::vector<real*>, Cell::NumFaces> displacementsPtrs{};

    std::array<std::vector<real*>, 4> derivatives{};
    std::array<std::vector<real*>, 4> dofsPtrs{};
    std::array<std::vector<real*>, 4> fluxTargetPtrs{};
    std::array<std::vector<real*>, 4> neighPtrs{};
    std::array<std::vector<real*>, 4> t{};
    std::array<std::vector<real*>, 4> tInv{};
    std::array<std::vector<real*>, 4> fsgdata{};

    std::array<std::vector<double>, 4> invImpedances{};

    std::array<std::size_t, Cell::NumFaces> counter{};

    size_t nodalAvgDisplacementsCounter{0};

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < 4; ++face) {
        if (dataHost.get<LTS::CellInformation>().faceTypes[face] == FaceType::FreeSurfaceGravity) {
          assert(dataHost.get<LTS::FaceDisplacementsDevice>()[face] != nullptr);
          cellIndices[face].push_back(cell);

          derivatives[face].push_back(dQPtrs_[cell]);
          dofsPtrs[face].push_back(static_cast<real*>(data.get<LTS::Dofs>()));
          fluxTargetPtrs[face].push_back(fluxTargets_[cell]);
          neighPtrs[face].push_back(
              reinterpret_cast<real*>(&data.get<LTS::NeighboringIntegration>()));
          displacementsPtrs[face].push_back(dataHost.get<LTS::FaceDisplacementsDevice>()[face]);
          t[face].push_back(dataHost.get<LTS::BoundaryMappingDevice>()[face].dataT);
          tInv[face].push_back(dataHost.get<LTS::BoundaryMappingDevice>()[face].dataTinv);

          real* displ{&nodalAvgDisplacements[nodalAvgDisplacementsCounter]};
          nodalAvgDisplacementsPtrs[face].push_back(displ);
          nodalAvgDisplacementsCounter += NodalAvgDisplacementsSize;

          fsgdata[face].push_back(dataHost.get<LTS::BoundaryMappingDevice>()[face].fsgData);

          const auto rho = dataHost.get<LTS::Material>().local->getDensity();
          const auto lambda = dataHost.get<LTS::Material>().local->getLambdaBar();

          auto* fsgdataPtrHost = dataHost.get<LTS::BoundaryMapping>()[face].fsgData;

          fsgdataPtrHost[0] = 1.0 / std::sqrt(rho * lambda);
          fsgdataPtrHost[1] = rho * g_;
          fsgdataPtrHost[2] = rho;

          ++counter[face];
        }
      }
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!cellIndices[face].empty()) {
        const ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, face);
        checkKey(key);
        (*currentIndicesTable_)[key].set(inner_keys::Indices::Id::Cells, cellIndices[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, derivatives[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FluxTarget, fluxTargetPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::T, t[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Tinv, tInv[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::FaceDisplacement, displacementsPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::NodalAvgDisplacements,
                                  nodalAvgDisplacementsPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::FSGData, fsgdata[face]);
      }
    }
  }
}

void LocalIntegrationRecorder::recordDirichletBc() {
  const auto size = currentLayer_->size();
  if (size > 0) {
    std::array<std::vector<real*>, Cell::NumFaces> dofsPtrs{};
    std::array<std::vector<real*>, Cell::NumFaces> fluxTargetPtrs{};
    std::array<std::vector<real*>, Cell::NumFaces> neighPtrs{};

    std::array<std::vector<real*>, Cell::NumFaces> dirichletOffsetPtrs{};

    std::array<std::size_t, 4> counter{};

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < 4; ++face) {
        if (dataHost.get<LTS::CellInformation>().faceTypes[face] == FaceType::Dirichlet) {

          dofsPtrs[face].push_back(static_cast<real*>(data.get<LTS::Dofs>()));
          fluxTargetPtrs[face].push_back(fluxTargets_[cell]);
          neighPtrs[face].push_back(
              reinterpret_cast<real*>(&data.get<LTS::NeighboringIntegration>()));

          dirichletOffsetPtrs[face].push_back(
              dataHost.get<LTS::BoundaryMappingDevice>()[face].dirichletOffset);

          ++counter[face];
        }
      }
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!dofsPtrs[face].empty()) {
        const ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
        checkKey(key);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FluxTarget, fluxTargetPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::DirichletOffset, dirichletOffsetPtrs[face]);
      }
    }
  }
}

void LocalIntegrationRecorder::recordAnalyticalBc(LTS::Layer& layer) {
  const auto size = currentLayer_->size();
  if (size > 0) {
    std::array<std::vector<real*>, Cell::NumFaces> dofsPtrs{};
    std::array<std::vector<real*>, Cell::NumFaces> fluxTargetPtrs{};
    std::array<std::vector<real*>, Cell::NumFaces> neighPtrs{};
    std::array<std::vector<unsigned>, Cell::NumFaces> cellIndices{};
    std::array<std::vector<real*>, Cell::NumFaces> analytical{};

    real* analyticScratch =
        reinterpret_cast<real*>(layer.var<LTS::AnalyticScratch>(AllocationPlace::Device));

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < 4; ++face) {
        if (dataHost.get<LTS::CellInformation>().faceTypes[face] == FaceType::Analytical) {
          cellIndices[face].push_back(cell);
          dofsPtrs[face].push_back(data.get<LTS::Dofs>());
          fluxTargetPtrs[face].push_back(fluxTargets_[cell]);
          neighPtrs[face].push_back(
              reinterpret_cast<real*>(&data.get<LTS::NeighboringIntegration>()));
          analytical[face].push_back(analyticScratch + cell * tensor::INodal::size());
        }
      }
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!cellIndices[face].empty()) {
        const ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
        checkKey(key);
        (*currentIndicesTable_)[key].set(inner_keys::Indices::Id::Cells, cellIndices[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FluxTarget, fluxTargetPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Analytical, analytical[face]);
      }
    }
  }
}

// NOLINTEND (-misc-const-correctness)
