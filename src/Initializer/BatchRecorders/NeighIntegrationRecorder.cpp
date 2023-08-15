#include "Recorders.h"
#include "utils/logger.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>
#include <Common/cellconfigconv.hpp>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

template<typename Config>
void NeighIntegrationRecorder<Config>::record(LTS<Config>& handler, Layer& layer) {
  typename kernels::NeighborData<Config>::Loader loader;
  loader.load(handler, layer);
  setUpContext(handler, layer, loader);
  idofsAddressRegistry.clear();

  recordDofsTimeEvaluation();
  recordNeighbourFluxIntegrals();
}

template<typename Config>
void NeighIntegrationRecorder<Config>::recordDofsTimeEvaluation() {
  RealT*(*faceNeighbors)[4] = currentLayer->var(currentHandler->faceNeighbors);
  RealT* integratedDofsScratch =
      static_cast<RealT*>(currentLayer->getScratchpadMemory(currentHandler->integratedDofsScratch));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<RealT*> ltsIDofsPtrs{};
    std::vector<RealT*> ltsDerivativesPtrs{};
    std::vector<RealT*> gtsDerivativesPtrs{};
    std::vector<RealT*> gtsIDofsPtrs{};

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);

      for (unsigned face = 0; face < 4; ++face) {
        RealT* neighbourBuffer = faceNeighbors[cell][face];

        // check whether a neighbour element idofs has not been counted twice
        if ((idofsAddressRegistry.find(neighbourBuffer) == idofsAddressRegistry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighbourBuffer != nullptr) {

            if (data.cellInformation.faceTypes[face] != FaceType::outflow &&
                data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {

              bool isNeighbProvidesDerivatives = ((data.cellInformation.ltsSetup >> face) % 2) == 1;

              if (isNeighbProvidesDerivatives) {
                RealT* NextTempIDofsPtr = &integratedDofsScratch[integratedDofsAddressCounter];

                bool isGtsNeigbour = ((data.cellInformation.ltsSetup >> (face + 4)) % 2) == 1;
                if (isGtsNeigbour) {

                  idofsAddressRegistry[neighbourBuffer] = NextTempIDofsPtr;
                  gtsIDofsPtrs.push_back(NextTempIDofsPtr);
                  gtsDerivativesPtrs.push_back(neighbourBuffer);

                } else {
                  idofsAddressRegistry[neighbourBuffer] = NextTempIDofsPtr;
                  ltsIDofsPtrs.push_back(NextTempIDofsPtr);
                  ltsDerivativesPtrs.push_back(neighbourBuffer);
                }
                integratedDofsAddressCounter += Yateto<Config>::Tensor::I::size();
              } else {
                idofsAddressRegistry[neighbourBuffer] = neighbourBuffer;
              }
            }
          }
        }
      }
    }

    if (!gtsIDofsPtrs.empty()) {
      ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Derivatives, gtsDerivativesPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, gtsIDofsPtrs);
    }

    if (!ltsIDofsPtrs.empty()) {
      ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Derivatives, ltsDerivativesPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, ltsIDofsPtrs);
    }
  }
}


template<typename Config>
void NeighIntegrationRecorder<Config>::recordNeighbourFluxIntegrals() {
  RealT*(*faceNeighbors)[4] = currentLayer->var(currentHandler->faceNeighbors);

  std::array<std::vector<RealT*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicDofs {};
  std::array<std::vector<RealT*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicIDofs {};
  std::array<std::vector<RealT*>[*FaceRelations::Count], *FaceId::Count> regularPeriodicAminusT {};

  std::array<std::vector<RealT*>[*DrFaceRelations::Count], *FaceId::Count> drDofs {};
  std::array<std::vector<RealT*>[*DrFaceRelations::Count], *FaceId::Count> drGodunov {};
  std::array<std::vector<RealT*>[*DrFaceRelations::Count], *FaceId::Count> drFluxSolver {};

  CellDRMapping(*drMapping)[4] = currentLayer->var(currentHandler->drMapping);

  const auto size = currentLayer->getNumberOfCells();
  for (unsigned cell = 0; cell < size; ++cell) {
    auto data = currentLoader->entry(cell);
    for (unsigned int face = 0; face < 4; face++) {
      switch (data.cellInformation.faceTypes[face]) {
      case FaceType::regular:
        [[fallthrough]];
      case FaceType::periodic: {
        // compute face type relation

        RealT* neighbourBufferPtr = faceNeighbors[cell][face];
        // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
        if (neighbourBufferPtr != nullptr) {
          unsigned faceRelation = data.cellInformation.faceRelations[face][1] +
                                  3 * data.cellInformation.faceRelations[face][0] + 12 * face;

          assert((*FaceRelations::Count) > faceRelation &&
                 "incorrect face relation count has been detected");

          regularPeriodicDofs[face][faceRelation].push_back(static_cast<RealT*>(data.dofs));
          regularPeriodicIDofs[face][faceRelation].push_back(
              idofsAddressRegistry[neighbourBufferPtr]);
          regularPeriodicAminusT[face][faceRelation].push_back(
              static_cast<RealT*>(data.neighIntegrationOnDevice.nAmNm1[face]));
        }
        break;
      }
      case FaceType::freeSurface: {
        break;
      }
      case FaceType::dynamicRupture: {
        unsigned faceRelation = drMapping[cell][face].side + 4 * drMapping[cell][face].faceRelation;
        assert((*DrFaceRelations::Count) > faceRelation &&
               "incorrect face relation count in dyn. rupture has been detected");
        drDofs[face][faceRelation].push_back(static_cast<RealT*>(data.dofs));
        drGodunov[face][faceRelation].push_back(drMapping[cell][face].godunov);
        drFluxSolver[face][faceRelation].push_back(drMapping[cell][face].fluxSolver);

        break;
      }
      case FaceType::outflow:
        [[fallthrough]];
      case FaceType::analytical:
        [[fallthrough]];
      case FaceType::freeSurfaceGravity:
        [[fallthrough]];
      case FaceType::dirichlet: {
        // Do not need to compute anything in the neighboring macro-kernel
        // for outflow, analytical, freeSurfaceGravity and dirichlet
        // boundary conditions
        break;
      }
      default: {
        logError() << "unknown boundary condition type: "
                   << static_cast<int>(data.cellInformation.faceTypes[face]);
      }
      }
    }
  }

  for (unsigned int face = 0; face < 4; ++face) {
    // regular and periodic
    for (size_t faceRelation = 0; faceRelation < (*FaceRelations::Count); ++faceRelation) {
      if (!regularPeriodicDofs[face][faceRelation].empty()) {
        ConditionalKey key(*KernelNames::NeighborFlux,
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
        ConditionalKey key(
            *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);
        checkKey(key);

        (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, drDofs[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Godunov, drGodunov[face][faceRelation]);
        (*currentTable)[key].set(inner_keys::Wp::Id::FluxSolver, drFluxSolver[face][faceRelation]);
      }
    }
  }
}

const DeclareForAllConfigs<NeighIntegrationRecorder> declNeighIntegrationRecorder;
