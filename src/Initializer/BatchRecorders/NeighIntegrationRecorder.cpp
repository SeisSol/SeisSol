// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GeneratedCode/tensor.h"
#include "Kernels/Interface.h"
#include "Recorders.h"
#include "utils/logger.h"
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void NeighIntegrationRecorder::record(LTS::Layer& layer) {
  setUpContext(layer);
  idofsAddressRegistry.clear();

  recordDofsTimeEvaluation();
  recordNeighborFluxIntegrals();
}

void NeighIntegrationRecorder::recordDofsTimeEvaluation() {
  real*(*faceNeighborsDevice)[4] = currentLayer->var<LTS::FaceNeighborsDevice>();
  real* integratedDofsScratch =
      static_cast<real*>(currentLayer->var<LTS::IntegratedDofsScratch>(AllocationPlace::Device));

  const auto size = currentLayer->size();
  if (size > 0) {
    std::vector<real*> ltsIDofsPtrs{};
    std::vector<real*> ltsDerivativesPtrs{};
    std::vector<real*> gtsDerivativesPtrs{};
    std::vector<real*> gtsIDofsPtrs{};

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto dataHost = currentLayer->cellRef(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        real* neighborBuffer = faceNeighborsDevice[cell][face];

        // check whether a neighbor element idofs has not been counted twice
        if ((idofsAddressRegistry.find(neighborBuffer) == idofsAddressRegistry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighborBuffer != nullptr) {
            if (dataHost.get<LTS::CellInformation>().faceTypes[face] != FaceType::Outflow &&
                dataHost.get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {

              const bool isNeighbProvidesDerivatives =
                  ((dataHost.get<LTS::CellInformation>().ltsSetup >> face) % 2) == 1;

              if (isNeighbProvidesDerivatives) {
                real* nextTempIDofsPtr = &integratedDofsScratch[integratedDofsAddressCounter];

                const bool isGtsNeighbor =
                    ((dataHost.get<LTS::CellInformation>().ltsSetup >> (face + 4)) % 2) == 1;
                if (isGtsNeighbor) {

                  idofsAddressRegistry[neighborBuffer] = nextTempIDofsPtr;
                  gtsIDofsPtrs.push_back(nextTempIDofsPtr);
                  gtsDerivativesPtrs.push_back(neighborBuffer);

                } else {
                  idofsAddressRegistry[neighborBuffer] = nextTempIDofsPtr;
                  ltsIDofsPtrs.push_back(nextTempIDofsPtr);
                  ltsDerivativesPtrs.push_back(neighborBuffer);
                }
                integratedDofsAddressCounter += tensor::I::size();
              } else {
                idofsAddressRegistry[neighborBuffer] = neighborBuffer;
              }
            }
          }
        }
      }
    }

    if (!gtsIDofsPtrs.empty()) {
      const ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Derivatives, gtsDerivativesPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, gtsIDofsPtrs);
    }

    if (!ltsIDofsPtrs.empty()) {
      const ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Derivatives, ltsDerivativesPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, ltsIDofsPtrs);
    }
  }
}

void NeighIntegrationRecorder::recordNeighborFluxIntegrals() {
  real*(*faceNeighborsDevice)[4] = currentLayer->var<LTS::FaceNeighborsDevice>();

  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicDofs {};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicIDofs {};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicAminusT {};

  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drDofs {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drGodunov {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drFluxSolver {};

  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularDofsExt {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drDofsExt {};

  CellDRMapping(*drMappingDevice)[4] = currentLayer->var<LTS::DRMappingDevice>();

#ifdef USE_VISCOELASTIC2
  auto* dofsExt = currentLayer->var<LTS::DofsExtScratch>(AllocationPlace::Device);
#endif

  const auto size = currentLayer->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto data = currentLayer->cellRef(cell, AllocationPlace::Device);
    auto dataHost = currentLayer->cellRef(cell, AllocationPlace::Host);

    for (std::size_t face = 0; face < Cell::NumFaces; face++) {
      switch (dataHost.get<LTS::CellInformation>().faceTypes[face]) {
      case FaceType::Regular:
        [[fallthrough]];
      case FaceType::Periodic: {
        // compute face type relation

        real* neighborBufferPtr = faceNeighborsDevice[cell][face];
        // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
        if (neighborBufferPtr != nullptr) {
          const auto faceRelation =
              dataHost.get<LTS::CellInformation>().faceRelations[face][1] +
              3 * dataHost.get<LTS::CellInformation>().faceRelations[face][0] + 12 * face;

          assert((*FaceRelations::Count) > faceRelation &&
                 "incorrect face relation count has been detected");

          regularPeriodicDofs[face][faceRelation].push_back(
              static_cast<real*>(data.get<LTS::Dofs>()));
          regularPeriodicIDofs[face][faceRelation].push_back(
              idofsAddressRegistry[neighborBufferPtr]);
          regularPeriodicAminusT[face][faceRelation].push_back(
              reinterpret_cast<real*>(&data.get<LTS::NeighboringIntegration>()));
#ifdef USE_VISCOELASTIC2
          regularDofsExt[face][faceRelation].push_back(static_cast<real*>(dofsExt) +
                                                       tensor::Qext::size() * cell);
#endif
        }
        break;
      }
      case FaceType::FreeSurface: {
        break;
      }
      case FaceType::DynamicRupture: {
        const auto faceRelation =
            drMappingDevice[cell][face].side + 4 * drMappingDevice[cell][face].faceRelation;
        assert((*DrFaceRelations::Count) > faceRelation &&
               "incorrect face relation count in dyn. rupture has been detected");
        drDofs[face][faceRelation].push_back(static_cast<real*>(data.get<LTS::Dofs>()));
        drGodunov[face][faceRelation].push_back(drMappingDevice[cell][face].godunov);
        drFluxSolver[face][faceRelation].push_back(drMappingDevice[cell][face].fluxSolver);
#ifdef USE_VISCOELASTIC2
        drDofsExt[face][faceRelation].push_back(static_cast<real*>(dofsExt) +
                                                tensor::Qext::size() * cell);
#endif
        break;
      }
      case FaceType::Outflow:
        [[fallthrough]];
      case FaceType::Analytical:
        [[fallthrough]];
      case FaceType::FreeSurfaceGravity:
        [[fallthrough]];
      case FaceType::Dirichlet: {
        // Do not need to compute anything in the neighboring macro-kernel
        // for outflow, analytical, freeSurfaceGravity and dirichlet
        // boundary conditions
        break;
      }
      default: {
        logError() << "unknown boundary condition type: "
                   << static_cast<int>(dataHost.get<LTS::CellInformation>().faceTypes[face]);
      }
      }
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // regular and periodic
    for (size_t faceRelation = 0; faceRelation < (*FaceRelations::Count); ++faceRelation) {
      if (!regularPeriodicDofs[face][faceRelation].empty()) {
        const ConditionalKey key(*KernelNames::NeighborFlux,
                                 (FaceKinds::Regular || FaceKinds::Periodic),
                                 face,
                                 faceRelation);
        checkKey(key);

        (*currentTable)[key].set(inner_keys::Wp::Id::Idofs,
                                 regularPeriodicIDofs[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, regularPeriodicDofs[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::NeighborIntegrationData,
                                 regularPeriodicAminusT[face][faceRelation]);
#ifdef USE_VISCOELASTIC2
        (*currentTable)[key].set(inner_keys::Wp::Id::DofsExt, regularDofsExt[face][faceRelation]);
#endif
      }
    }

    // dynamic rupture
    for (std::size_t faceRelation = 0; faceRelation < (*DrFaceRelations::Count); ++faceRelation) {
      if (!drDofs[face][faceRelation].empty()) {
        const ConditionalKey key(
            *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);
        checkKey(key);

        (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, drDofs[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Godunov, drGodunov[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::FluxSolver, drFluxSolver[face][faceRelation]);
#ifdef USE_VISCOELASTIC2
        (*currentTable)[key].set(inner_keys::Wp::Id::DofsExt, drDofsExt[face][faceRelation]);
#endif
      }
    }
  }
}
