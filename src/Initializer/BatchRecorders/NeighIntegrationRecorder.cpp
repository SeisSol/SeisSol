// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Constants.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/ConfigBoundary.h"
#include "Kernels/SolverSelector.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yateto.h>

// NOLINTBEGIN (-misc-const-correctness)

using namespace seissol::initializer;
using namespace seissol::recording;

template <typename Cfg>
void NeighIntegrationRecorder<Cfg>::record(LTS::Layer& layer) {
  setUpContext(layer);
  idofsAddressRegistry_.clear();

  recordDofsTimeEvaluation();
  recordConfigBoundaryBatches();
  recordNeighborFluxIntegrals();
}

template <typename Cfg>
void NeighIntegrationRecorder<Cfg>::recordDofsTimeEvaluation() {
  auto* faceNeighborsDevice = currentLayer_->var<LTS::FaceNeighborsDevice>();
  real* integratedDofsScratch = static_cast<real*>(
      currentLayer_->var<LTS::IntegratedDofsScratch>(Cfg(), AllocationPlace::Device));

  const auto size = currentLayer_->size();
  if (size > 0) {
    std::vector<real*> ltsIDofsPtrs{};
    std::vector<real*> ltsDerivativesPtrs{};
    std::vector<real*> gtsDerivativesPtrs{};
    std::vector<real*> gtsIDofsPtrs{};

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto dataHost = currentLayer_->cellRef<Cfg>(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        auto* neighborBuffer = static_cast<real*>(faceNeighborsDevice[cell][face]);

        const auto& cellInformation = dataHost.template get<LTS::CellInformation>();
        if (cellInformation.faceTypes[face] == FaceType::Regular && neighborBuffer != nullptr &&
            cellInformation.neighborConfigIds[face] != configIdOf<Cfg>()) {
          recordConfigBoundaryFace(cell, face);
          continue;
        }

        // check whether a neighbor element idofs has not been counted twice
        if ((idofsAddressRegistry_.find(neighborBuffer) == idofsAddressRegistry_.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighborBuffer != nullptr) {
            if (dataHost.template get<LTS::CellInformation>().faceTypes[face] ==
                FaceType::Regular) {

              const bool isNeighbProvidesDerivatives =
                  dataHost.template get<LTS::CellInformation>().ltsSetup.neighborBuffer(face) ==
                  BufferType::Derivatives;

              if (isNeighbProvidesDerivatives) {
                real* nextTempIDofsPtr = &integratedDofsScratch[integratedDofsAddressCounter_];

                const bool isGtsNeighbor =
                    dataHost.template get<LTS::CellInformation>().ltsSetup.neighborGTSRelation(
                        face);
                if (isGtsNeighbor) {

                  // might effectively not occur anymore; but we keep it anyways (for now)
                  // (formerly, this was due to AccumulatedIntegrals and StepIntegrals not being
                  // allowed to coexist; so we used Derivatives instead)

                  idofsAddressRegistry_[neighborBuffer] = nextTempIDofsPtr;
                  gtsIDofsPtrs.push_back(nextTempIDofsPtr);
                  gtsDerivativesPtrs.push_back(neighborBuffer);

                } else {
                  idofsAddressRegistry_[neighborBuffer] = nextTempIDofsPtr;
                  ltsIDofsPtrs.push_back(nextTempIDofsPtr);
                  ltsDerivativesPtrs.push_back(neighborBuffer);
                }
                integratedDofsAddressCounter_ += kernels::SolverOf<Cfg>::IntegralsSize;
              } else {
                idofsAddressRegistry_[neighborBuffer] = neighborBuffer;
              }
            }
          }
        }
      }
    }

    if (!gtsIDofsPtrs.empty()) {
      const ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
      checkKey(key);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, gtsDerivativesPtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, gtsIDofsPtrs);
    }

    if (!ltsIDofsPtrs.empty()) {
      const ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
      checkKey(key);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, ltsDerivativesPtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, ltsIDofsPtrs);
    }
  }
}

template <typename Cfg>
void* NeighIntegrationRecorder<Cfg>::allocateConfigBoundaryScratch(std::size_t bytes) {
  auto* scratch = static_cast<std::uint8_t*>(
      currentLayer_->var<LTS::ConfigBoundaryScratch>(AllocationPlace::Device));
  void* region = scratch + configBoundaryScratchCounter_;
  configBoundaryScratchCounter_ += kernels::configboundary::alignScratch(bytes);
  return region;
}

template <typename Cfg>
void NeighIntegrationRecorder<Cfg>::recordConfigBoundaryFace(std::size_t cell, std::size_t face) {
  const auto& cellInformation = currentLayer_->var<LTS::CellInformation>()[cell];
  const void* neighborBuffer = currentLayer_->var<LTS::FaceNeighborsDevice>()[cell][face];
  const auto neighborConfig = cellInformation.neighborConfigIds[face];
  const auto side = static_cast<std::size_t>(cellInformation.faceRelations[face][0]);

  dispatchConfig(neighborConfig, [&](auto neighborCfg) {
    using NeighborCfg = decltype(neighborCfg);
    if constexpr (!std::is_same_v<NeighborCfg, Cfg> &&
                  kernels::DeviceConvertible<Cfg, NeighborCfg>) {
      // the neighbor in the canonical form, once per neighbor
      auto found = canonicalRegistry_.find(neighborBuffer);
      if (found == canonicalRegistry_.end()) {
        auto& batches = configBoundaryBatches_.at(neighborConfig);
        void* integral = const_cast<void*>(neighborBuffer);
        if (cellInformation.ltsSetup.neighborBuffer(face) == BufferType::Derivatives) {
          // the time integral of the neighbor in its configuration first
          const auto relation = cellInformation.ltsSetup.neighborGTSRelation(face) ? 0 : 1;
          integral = allocateConfigBoundaryScratch(kernels::SolverOf<NeighborCfg>::IntegralsSize *
                                                   sizeof(Real<NeighborCfg>));
          batches.derivatives[relation].push_back(const_cast<void*>(neighborBuffer));
          batches.integrals[relation].push_back(integral);
        }
        auto* canonical = static_cast<double*>(allocateConfigBoundaryScratch(
            tensor::canonicalI<NeighborCfg>::size() * sizeof(double)));
        batches.toCanonical.push_back(integral);
        batches.canonical.push_back(canonical);
        found = canonicalRegistry_.emplace(neighborBuffer, canonical).first;
      }

      // converted into the configuration of the layer, once per neighbor and side
      const auto key = std::make_pair(neighborBuffer, side);
      if (convertedRegistry_.find(key) == convertedRegistry_.end()) {
        auto* converted = static_cast<real*>(
            allocateConfigBoundaryScratch(kernels::SolverOf<Cfg>::IntegralsSize * sizeof(real)));
        fromCanonical_.at(side).push_back(found->second);
        converted_.at(side).push_back(converted);
        convertedRegistry_.emplace(key, converted);
      }
    } else {
      logError() << "The faces between the configurations"
                 << configName(configValue(neighborConfig)) << "and"
                 << configName(configValue(configIdOf<Cfg>())) << "cannot be computed on GPUs.";
    }
  });
}

template <typename Cfg>
void NeighIntegrationRecorder<Cfg>::recordConfigBoundaryBatches() {
  for (ConfigId neighborConfig = 0; neighborConfig < configBoundaryBatches_.size();
       ++neighborConfig) {
    auto& batches = configBoundaryBatches_[neighborConfig];
    if (batches.toCanonical.empty()) {
      continue;
    }
    dispatchConfig(neighborConfig, [&](auto neighborCfg) {
      using NeighborReal = Real<decltype(neighborCfg)>;
      const auto typed = [](const std::vector<void*>& pointers) {
        std::vector<NeighborReal*> result(pointers.size());
        std::transform(pointers.begin(), pointers.end(), result.begin(), [](void* pointer) {
          return static_cast<NeighborReal*>(pointer);
        });
        return result;
      };

      for (const bool gts : {true, false}) {
        const auto relation = gts ? 0 : 1;
        if (!batches.integrals[relation].empty()) {
          const auto key = kernels::configboundary::timeKey(neighborConfig, gts);
          checkKey(key);
          auto derivatives = typed(batches.derivatives[relation]);
          auto integrals = typed(batches.integrals[relation]);
          (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, derivatives);
          (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, integrals);
        }
      }

      const auto key = kernels::configboundary::toCanonicalKey(neighborConfig);
      checkKey(key);
      auto integrals = typed(batches.toCanonical);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, integrals);
      (*currentTable_)[key].set(inner_keys::Wp::Id::CanonicalIdofs, batches.canonical);
    });
  }

  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    if (!converted_[side].empty()) {
      const auto key = kernels::configboundary::fromCanonicalKey(side);
      checkKey(key);
      (*currentTable_)[key].set(inner_keys::Wp::Id::CanonicalIdofs, fromCanonical_[side]);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, converted_[side]);
    }
  }
}

template <typename Cfg>
void NeighIntegrationRecorder<Cfg>::recordNeighborFluxIntegrals() {
  auto* faceNeighborsDevice = currentLayer_->var<LTS::FaceNeighborsDevice>();

  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicDofs {};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicIDofs {};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicAminusT {};

  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drDofs {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drGodunov {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drFluxSolver {};

  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regularDofsExt {};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> drDofsExt {};

  const auto* drMappingDevice = currentLayer_->var<LTS::DRMappingDevice>(Cfg());

  auto* dofsExt = currentLayer_->var<LTS::DofsExtScratch>(Cfg(), AllocationPlace::Device);

  const auto size = currentLayer_->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto data = currentLayer_->cellRef<Cfg>(cell, AllocationPlace::Device);
    auto dataHost = currentLayer_->cellRef<Cfg>(cell, AllocationPlace::Host);

    for (std::size_t face = 0; face < Cell::NumFaces; face++) {
      switch (dataHost.template get<LTS::CellInformation>().faceTypes[face]) {
      case FaceType::Regular: {
        // compute face type relation

        auto* neighborBufferPtr = static_cast<real*>(faceNeighborsDevice[cell][face]);
        // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
        if (neighborBufferPtr != nullptr) {
          const auto faceRelation =
              dataHost.template get<LTS::CellInformation>().faceRelations[face][0] + 4 * face;

          assert((*FaceRelations::Count) > faceRelation &&
                 "incorrect face relation count has been detected");

          regularPeriodicDofs[face][faceRelation].push_back(
              static_cast<real*>(data.template get<LTS::Dofs>()));
          const auto& cellInformation = dataHost.template get<LTS::CellInformation>();
          if (cellInformation.neighborConfigIds[face] != configIdOf<Cfg>()) {
            // the time integral of the neighbor, converted into the configuration of the layer
            regularPeriodicIDofs[face][faceRelation].push_back(
                convertedRegistry_.at({neighborBufferPtr, cellInformation.faceRelations[face][0]}));
          } else {
            regularPeriodicIDofs[face][faceRelation].push_back(
                idofsAddressRegistry_[neighborBufferPtr]);
          }
          regularPeriodicAminusT[face][faceRelation].push_back(
              reinterpret_cast<real*>(&data.template get<LTS::NeighboringIntegration>()));
          if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
            regularDofsExt[face][faceRelation].push_back(static_cast<real*>(dofsExt) +
                                                         kernels::size<tensor::Qext<Cfg>>() * cell);
          }
        }
        break;
      }
      case FaceType::DynamicRupture: {
        const std::size_t faceRelation =
            drMappingDevice[cell][face].side + 4 * drMappingDevice[cell][face].faceRelation;
        assert((*DrFaceRelations::Count) > faceRelation &&
               "incorrect face relation count in dyn. rupture has been detected");
        assert(drMappingDevice[cell][face].side == face &&
               "the batched neighbor integral only visits the side of the face itself");
        drDofs[face][faceRelation].push_back(static_cast<real*>(data.template get<LTS::Dofs>()));
        drGodunov[face][faceRelation].push_back(drMappingDevice[cell][face].godunov);
        drFluxSolver[face][faceRelation].push_back(drMappingDevice[cell][face].fluxSolver);
        if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
          drDofsExt[face][faceRelation].push_back(static_cast<real*>(dofsExt) +
                                                  kernels::size<tensor::Qext<Cfg>>() * cell);
        }
        break;
      }
      case FaceType::FreeSurface:
        [[fallthrough]];
      case FaceType::Outflow:
        [[fallthrough]];
      case FaceType::Analytical:
        [[fallthrough]];
      case FaceType::FreeSurfaceGravity:
        [[fallthrough]];
      case FaceType::Dirichlet: {
        // Do not need to compute anything in the neighboring macro-kernel
        // for most boundary conditions
        break;
      }
      default: {
        logError() << "unknown boundary condition type: "
                   << static_cast<int>(
                          dataHost.template get<LTS::CellInformation>().faceTypes[face]);
      }
      }
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // regular and periodic
    for (size_t faceRelation = 0; faceRelation < (*FaceRelations::Count); ++faceRelation) {
      if (!regularPeriodicDofs[face][faceRelation].empty()) {
        const ConditionalKey key(
            *KernelNames::NeighborFlux, *FaceKinds::Regular, face, faceRelation);
        checkKey(key);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs,
                                  regularPeriodicIDofs[face][faceRelation]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs,
                                  regularPeriodicDofs[face][faceRelation]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::NeighborIntegrationData,
                                  regularPeriodicAminusT[face][faceRelation]);
        if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
          (*currentTable_)[key].set(inner_keys::Wp::Id::DofsExt,
                                    regularDofsExt[face][faceRelation]);
        }
      }
    }

    // dynamic rupture
    for (std::size_t faceRelation = 0; faceRelation < (*DrFaceRelations::Count); ++faceRelation) {
      if (!drDofs[face][faceRelation].empty()) {
        const ConditionalKey key(
            *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);
        checkKey(key);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, drDofs[face][faceRelation]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Godunov, drGodunov[face][faceRelation]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::FluxSolver, drFluxSolver[face][faceRelation]);
        if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
          (*currentTable_)[key].set(inner_keys::Wp::Id::DofsExt, drDofsExt[face][faceRelation]);
        }
      }
    }
  }
}

#define SEISSOL_INSTANTIATE(Cfg) template class seissol::recording::NeighIntegrationRecorder<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

// NOLINTEND (-misc-const-correctness)
