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
#include "Kernels/FaceTraces.h"
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
  real* zeroTrace = currentLayer_->var<LTS::ZeroFaceTrace>(AllocationPlace::Device);
  real* sideOperators = currentLayer_->var<LTS::TraceSideOperators>(AllocationPlace::Device);

  const auto size = currentLayer_->size();
  if (size == 0) {
    return;
  }

  // the four faces of every cell, face after face, each with the operator of the side the cell is
  // on there (CellLocalInformation::traceSides): the kernel that applies the correctors of the four
  // at once reads them so (kernels::applyFaceCorrectorsBatched). A face that takes no corrector
  // (appliesFaceCorrector) reads a trace of zeros on both sides, which its corrector turns into
  // nothing.
  std::vector<real*> tracePtrs(Cell::NumFaces * size);
  std::vector<real*> neighborTracePtrs(Cell::NumFaces * size);
  std::vector<real*> sidePtrs(Cell::NumFaces * size);
  // what the faces lift into, see fluxTargetOf, and where a cell keeps the operators of its faces
  std::vector<real*> fluxTargetPtrs(size);
  std::vector<real*> localPtrs(size);

  for (std::size_t cell = 0; cell < size; ++cell) {
    auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
    fluxTargetPtrs[cell] = fluxTargetOf(*currentLayer_, cell);
    localPtrs[cell] = reinterpret_cast<real*>(&data.get<LTS::LocalIntegration>());

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto index = face * size + cell;
      sidePtrs[index] =
          sideOperators + cellInformation[cell].traceSides[face] * kernels::TraceSideOperatorSize;
      if (appliesFaceCorrector(cellInformation[cell].faceTypes[face])) {
        tracePtrs[index] = faceTraces[cell][face];
        neighborTracePtrs[index] = neighborTraces_[cell][face];
      } else {
        tracePtrs[index] = zeroTrace;
        neighborTracePtrs[index] = zeroTrace;
      }
    }
  }

  const ConditionalKey key(*KernelNames::FaceCorrector, *ComputationKind::None);
  checkKey(key);
  (*currentTable_)[key].set(inner_keys::Wp::Id::FaceTrace, tracePtrs);
  (*currentTable_)[key].set(inner_keys::Wp::Id::FaceTraceNeighbor, neighborTracePtrs);
  (*currentTable_)[key].set(inner_keys::Wp::Id::TraceSide, sidePtrs);
  (*currentTable_)[key].set(inner_keys::Wp::Id::FluxTarget, fluxTargetPtrs);
  (*currentTable_)[key].set(inner_keys::Wp::Id::LocalIntegrationData, localPtrs);
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
