// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/OperatorLayout.h"
#include "Recorders.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>
#include <yateto.h>

// NOLINTBEGIN (-misc-const-correctness)

using namespace seissol::initializer;
using namespace seissol::recording;

void NeighIntegrationRecorder::record(LTS::Layer& layer) {
  setUpContext(layer);

  recordNeighborTraces();
  recordFaceCorrectors();
  recordDynamicRuptureFluxes();
}

void NeighIntegrationRecorder::recordNeighborTraces() {
  const auto* faceNeighborsDevice = currentLayer_->var<LTS::FaceNeighborsDevice>();
  const auto* cellInformation = currentLayer_->var<LTS::CellInformation>();
  real* tracesScratch =
      static_cast<real*>(currentLayer_->var<LTS::IntegratedDofsScratch>(AllocationPlace::Device));
  real* zeroTrace = currentLayer_->var<LTS::ZeroFaceTrace>(AllocationPlace::Device);

  const auto size = currentLayer_->size();
  neighborTraces_.assign(size, {});

  // A neighbor that writes the traces of its time derivatives on a face takes a larger time step
  // (see LtsSetup); they are evaluated over the part of it the cell takes, into the scratch memory
  // of the layer.
  std::vector<real*> derivativesPtrs;
  std::vector<real*> tracePtrs;
  std::size_t evaluated = 0;

  for (std::size_t cell = 0; cell < size; ++cell) {
    const auto& ltsSetup = cellInformation[cell].ltsSetup;
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (cellInformation[cell].faceTypes[face] != FaceType::Regular) {
        neighborTraces_[cell][face] = zeroTrace;
      } else if (ltsSetup.neighborBuffer(face) != BufferType::Derivatives) {
        // the neighbor wrote the trace over the time step of the cell
        neighborTraces_[cell][face] = faceNeighborsDevice[cell][face];
      } else {
        assert(!ltsSetup.neighborGTSRelation(face));
        real* trace = tracesScratch + evaluated * kernels::Solver::FaceTraceSize;
        ++evaluated;
        neighborTraces_[cell][face] = trace;
        derivativesPtrs.push_back(faceNeighborsDevice[cell][face]);
        tracePtrs.push_back(trace);
      }
    }
  }

  if (!tracePtrs.empty()) {
    const ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
    checkKey(key);
    (*currentTable_)[key].set(inner_keys::Wp::Id::FaceBuffer, derivativesPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::FaceTraceNeighbor, tracePtrs);
  }
}

void NeighIntegrationRecorder::recordFaceCorrectors() {
  auto* faceTraces = currentLayer_->var<LTS::FaceTracesDevice>();
  const auto* cellInformation = currentLayer_->var<LTS::CellInformation>();

  // by the face and by the side of the face the cell is on (CellLocalInformation::traceSides)
  using PerFaceSide = std::array<std::array<std::vector<real*>, 2>, Cell::NumFaces>;
  PerFaceSide tracePtrs{};
  PerFaceSide neighborTracePtrs{};
  // what the faces lift into, see fluxTargetOf
  PerFaceSide fluxTargetPtrs{};
  PerFaceSide localPtrs{};

  const auto size = currentLayer_->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
    real* fluxTarget = fluxTargetOf(*currentLayer_, cell);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!appliesFaceCorrector(cellInformation[cell].faceTypes[face])) {
        continue;
      }
      const auto side = cellInformation[cell].traceSides[face];
      tracePtrs[face][side].push_back(faceTraces[cell][face]);
      neighborTracePtrs[face][side].push_back(neighborTraces_[cell][face]);
      fluxTargetPtrs[face][side].push_back(fluxTarget);
      localPtrs[face][side].push_back(reinterpret_cast<real*>(&data.get<LTS::LocalIntegration>()));
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (std::size_t side = 0; side < 2; ++side) {
      if (!tracePtrs[face][side].empty()) {
        const ConditionalKey key(*KernelNames::FaceCorrector, *ComputationKind::None, face, side);
        checkKey(key);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FaceTrace, tracePtrs[face][side]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FaceTraceNeighbor,
                                  neighborTracePtrs[face][side]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FluxTarget, fluxTargetPtrs[face][side]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::LocalIntegrationData, localPtrs[face][side]);
      }
    }
  }
}

void NeighIntegrationRecorder::recordDynamicRuptureFluxes() {
  // what the faces lift into, see fluxTargetOf
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drFluxTargets {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drDofs {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drGodunov {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drFluxSolver {};

  const auto* drMappingDevice = currentLayer_->var<LTS::DRMappingDevice>();
  const auto* cellInformation = currentLayer_->var<LTS::CellInformation>();

  const auto size = currentLayer_->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);

    for (std::size_t face = 0; face < Cell::NumFaces; face++) {
      if (cellInformation[cell].faceTypes[face] == FaceType::DynamicRupture) {
        const std::size_t faceRelation =
            drMappingDevice[cell][face].side + 4 * drMappingDevice[cell][face].faceRelation;
        assert((*DrFaceRelations::Count) > faceRelation &&
               "incorrect face relation count in dyn. rupture has been detected");
        drDofs[face][faceRelation].push_back(static_cast<real*>(data.get<LTS::Dofs>()));
        drFluxTargets[face][faceRelation].push_back(fluxTargetOf(*currentLayer_, cell));
        drGodunov[face][faceRelation].push_back(drMappingDevice[cell][face].godunov);
        drFluxSolver[face][faceRelation].push_back(drMappingDevice[cell][face].fluxSolver);
      }
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (std::size_t faceRelation = 0; faceRelation < (*DrFaceRelations::Count); ++faceRelation) {
      if (!drDofs[face][faceRelation].empty()) {
        const ConditionalKey key(
            *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);
        checkKey(key);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, drDofs[face][faceRelation]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FluxTarget,
                                  drFluxTargets[face][faceRelation]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Godunov, drGodunov[face][faceRelation]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FluxSolver, drFluxSolver[face][faceRelation]);
      }
    }
  }
}

// NOLINTEND (-misc-const-correctness)
