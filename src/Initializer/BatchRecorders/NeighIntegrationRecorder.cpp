// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

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
#include <tensor.h>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void NeighIntegrationRecorder::record(LTS& handler, Layer& layer) {
  kernels::NeighborData::Loader loader, loaderHost;
  loader.load(handler, layer, AllocationPlace::Device);
  loaderHost.load(handler, layer, AllocationPlace::Host);
  setUpContext(handler, layer, loader, loaderHost);
  idofsAddressRegistry.clear();

  recordDofsTimeEvaluation();
  recordNeighbourFluxIntegrals();
}

void NeighIntegrationRecorder::recordDofsTimeEvaluation() {
  real*(*faceNeighborsDevice)[4] = currentLayer->var(currentHandler->faceNeighborsDevice);
  real* integratedDofsScratch = static_cast<real*>(currentLayer->getScratchpadMemory(
      currentHandler->integratedDofsScratch, AllocationPlace::Device));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real*> ltsIDofsPtrs{};
    std::vector<real*> ltsDerivativesPtrs{};
    std::vector<real*> gtsDerivativesPtrs{};
    std::vector<real*> gtsIDofsPtrs{};

    for (unsigned cell = 0; cell < size; ++cell) {
      auto dataHost = currentLoaderHost->entry(cell);

      for (unsigned face = 0; face < 4; ++face) {
        real* neighbourBuffer = faceNeighborsDevice[cell][face];

        // check whether a neighbour element idofs has not been counted twice
        if ((idofsAddressRegistry.find(neighbourBuffer) == idofsAddressRegistry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighbourBuffer != nullptr) {
            if (dataHost.cellInformation().faceTypes[face] != FaceType::Outflow &&
                dataHost.cellInformation().faceTypes[face] != FaceType::DynamicRupture) {

              const bool isNeighbProvidesDerivatives =
                  ((dataHost.cellInformation().ltsSetup >> face) % 2) == 1;

              if (isNeighbProvidesDerivatives) {
                real* nextTempIDofsPtr = &integratedDofsScratch[integratedDofsAddressCounter];

                const bool isGtsNeigbour =
                    ((dataHost.cellInformation().ltsSetup >> (face + 4)) % 2) == 1;
                if (isGtsNeigbour) {

                  idofsAddressRegistry[neighbourBuffer] = nextTempIDofsPtr;
                  gtsIDofsPtrs.push_back(nextTempIDofsPtr);
                  gtsDerivativesPtrs.push_back(neighbourBuffer);

                } else {
                  idofsAddressRegistry[neighbourBuffer] = nextTempIDofsPtr;
                  ltsIDofsPtrs.push_back(nextTempIDofsPtr);
                  ltsDerivativesPtrs.push_back(neighbourBuffer);
                }
                integratedDofsAddressCounter += tensor::I::size();
              } else {
                idofsAddressRegistry[neighbourBuffer] = neighbourBuffer;
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

void NeighIntegrationRecorder::recordNeighbourFluxIntegrals() {
  real*(*faceNeighborsDevice)[4] = currentLayer->var(currentHandler->faceNeighborsDevice);

  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicDofs {};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicIDofs {};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicAminusT {};

  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drDofs {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drGodunov {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drFluxSolver {};

  CellDRMapping(*drMappingDevice)[4] = currentLayer->var(currentHandler->drMappingDevice);

  const auto size = currentLayer->getNumberOfCells();
  for (unsigned cell = 0; cell < size; ++cell) {
    auto data = currentLoader->entry(cell);
    auto dataHost = currentLoaderHost->entry(cell);

    for (unsigned int face = 0; face < 4; face++) {
      switch (dataHost.cellInformation().faceTypes[face]) {
      case FaceType::Regular:
        [[fallthrough]];
      case FaceType::Periodic: {
        // compute face type relation

        real* neighbourBufferPtr = faceNeighborsDevice[cell][face];
        // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
        if (neighbourBufferPtr != nullptr) {
          const unsigned faceRelation = dataHost.cellInformation().faceRelations[face][1] +
                                        3 * dataHost.cellInformation().faceRelations[face][0] +
                                        12 * face;

          assert((*FaceRelations::Count) > faceRelation &&
                 "incorrect face relation count has been detected");

          regularPeriodicDofs[face][faceRelation].push_back(static_cast<real*>(data.dofs()));
          regularPeriodicIDofs[face][faceRelation].push_back(
              idofsAddressRegistry[neighbourBufferPtr]);
          regularPeriodicAminusT[face][faceRelation].push_back(
              static_cast<real*>(data.neighboringIntegration().nAmNm1[face]));
        }
        break;
      }
      case FaceType::FreeSurface: {
        break;
      }
      case FaceType::DynamicRupture: {
        const unsigned faceRelation =
            drMappingDevice[cell][face].side + 4 * drMappingDevice[cell][face].faceRelation;
        assert((*DrFaceRelations::Count) > faceRelation &&
               "incorrect face relation count in dyn. rupture has been detected");
        drDofs[face][faceRelation].push_back(static_cast<real*>(data.dofs()));
        drGodunov[face][faceRelation].push_back(drMappingDevice[cell][face].godunov);
        drFluxSolver[face][faceRelation].push_back(drMappingDevice[cell][face].fluxSolver);

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
                   << static_cast<int>(dataHost.cellInformation().faceTypes[face]);
      }
      }
    }
  }

  for (unsigned int face = 0; face < 4; ++face) {
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
        (*currentTable)[key].set(inner_keys::Wp::Id::AminusT,
                                 regularPeriodicAminusT[face][faceRelation]);
      }
    }

    // dynamic rupture
    for (unsigned faceRelation = 0; faceRelation < (*DrFaceRelations::Count); ++faceRelation) {
      if (!drDofs[face][faceRelation].empty()) {
        const ConditionalKey key(
            *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);
        checkKey(key);

        (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, drDofs[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Godunov, drGodunov[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::FluxSolver, drFluxSolver[face][faceRelation]);
      }
    }
  }
}
