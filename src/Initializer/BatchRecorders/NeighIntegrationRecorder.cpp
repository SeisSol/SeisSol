#include "Recorders.h"
#include "utils/logger.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>


using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void NeighIntegrationRecorder::record(LTS &handler, Layer &layer) {
  kernels::NeighborData::Loader loader;
  loader.load(handler, layer);
  setUpContext(handler, layer, loader);
  idofsAddressRegistry.clear();

  recordDofsTimeEvaluation();
  recordNeighbourFluxIntegrals();
}


void NeighIntegrationRecorder::recordDofsTimeEvaluation() {
  real *(*faceNeighbors)[4] = currentLayer->var(currentHandler->faceNeighbors);
  real *idofsScratch = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsScratch));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real *> ltsIDofsPtrs{};
    std::vector<real *> ltsDerivativesPtrs{};
    std::vector<real *> gtsDerivativesPtrs{};
    std::vector<real *> gtsIDofsPtrs{};

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);

      for (unsigned face = 0; face < 4; ++face) {
        real *neighbourBuffer = faceNeighbors[cell][face];

        // check whether a neighbour element idofs has not been counted twice
        if ((idofsAddressRegistry.find(neighbourBuffer) == idofsAddressRegistry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighbourBuffer != nullptr) {

            if (data.cellInformation.faceTypes[face] != FaceType::outflow &&
                data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {

              bool isNeighbProvidesDerivatives = ((data.cellInformation.ltsSetup >> face) % 2) == 1;

              if (isNeighbProvidesDerivatives) {
                real *NextTempIDofsPtr = &idofsScratch[integratedDofsAddressCounter];

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
      ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
      checkKey(key);
      (*currentTable)[key].content[*EntityId::Derivatives] = new BatchPointers(gtsDerivativesPtrs);
      (*currentTable)[key].content[*EntityId::Idofs] = new BatchPointers(gtsIDofsPtrs);
    }

    if (!ltsIDofsPtrs.empty()) {
      ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
      checkKey(key);
      (*currentTable)[key].content[*EntityId::Derivatives] = new BatchPointers(ltsDerivativesPtrs);
      (*currentTable)[key].content[*EntityId::Idofs] = new BatchPointers(ltsIDofsPtrs);
    }
  }
}


void NeighIntegrationRecorder::recordNeighbourFluxIntegrals() {
  real *(*faceNeighbors)[4] = currentLayer->var(currentHandler->faceNeighbors);

  std::array<std::vector<real *>[*FaceRelations::Count], *FaceId::Count> regularPeriodicDofs {};
  std::array<std::vector<real *>[*FaceRelations::Count], *FaceId::Count> regularPeriodicIDofs {};
  std::array<std::vector<real *>[*FaceRelations::Count], *FaceId::Count> regularPeriodicAminusT {};

  std::array<std::vector<real *>[*DrFaceRelations::Count], *FaceId::Count> drDofs {};
  std::array<std::vector<real *>[*DrFaceRelations::Count], *FaceId::Count> drGodunov {};
  std::array<std::vector<real *>[*DrFaceRelations::Count], *FaceId::Count> drFluxSolver {};

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

          real *neighbourBufferPtr = faceNeighbors[cell][face];
          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighbourBufferPtr != nullptr) {
            unsigned faceRelation = data.cellInformation.faceRelations[face][1] +
                                    3 * data.cellInformation.faceRelations[face][0] + 12 * face;

            assert((*FaceRelations::Count) > faceRelation && "incorrect face relation count has been detected");

            regularPeriodicDofs[face][faceRelation].push_back(static_cast<real *>(data.dofs));
            regularPeriodicIDofs[face][faceRelation].push_back(idofsAddressRegistry[neighbourBufferPtr]);
            regularPeriodicAminusT[face][faceRelation].push_back(
                static_cast<real *>(data.neighIntegrationOnDevice.nAmNm1[face]));
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
          drDofs[face][faceRelation].push_back(static_cast<real *>(data.dofs));
          drGodunov[face][faceRelation].push_back(drMapping[cell][face].godunov);
          drFluxSolver[face][faceRelation].push_back(drMapping[cell][face].fluxSolver);

          break;
        }
        case FaceType::outflow:
          break;
        case FaceType::analytical: {
          logError() << "analytical boundary condition is not supported in batched computations";
          break;
        }
        case FaceType::freeSurfaceGravity: {
          logError() << "freeSurfaceGravity boundary condition is not supported in batched computations";
          break;
        }
        case FaceType::dirichlet: {
          logError() << "dirichlet boundary condition is not supported in batched computations";
          break;
        }
        default: {
          logError() << "unknown boundary condition type: " << static_cast<int>(data.cellInformation.faceTypes[face]);
        }
      }
    }
  }

  for (unsigned int face = 0; face < 4; ++face) {
    // regular and periodic
    for (size_t faceRelation = 0; faceRelation < (*FaceRelations::Count); ++faceRelation) {
      if (!regularPeriodicDofs[face][faceRelation].empty()) {
        ConditionalKey key(*KernelNames::NeighborFlux, (FaceKinds::Regular || FaceKinds::Periodic), face, faceRelation);
        checkKey(key);

        (*currentTable)[key].content[*EntityId::Idofs] = new BatchPointers(regularPeriodicIDofs[face][faceRelation]);
        (*currentTable)[key].content[*EntityId::Dofs] = new BatchPointers(regularPeriodicDofs[face][faceRelation]);
        (*currentTable)[key].content[*EntityId::AminusT] =
            new BatchPointers(regularPeriodicAminusT[face][faceRelation]);
      }
    }

    // dynamic rupture
    for (unsigned faceRelation = 0; faceRelation < (*DrFaceRelations::Count); ++faceRelation) {
      if (!drDofs[face][faceRelation].empty()) {
        ConditionalKey key(*KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);
        checkKey(key);

        (*currentTable)[key].content[*EntityId::Dofs] = new BatchPointers(drDofs[face][faceRelation]);
        (*currentTable)[key].content[*EntityId::Godunov] = new BatchPointers(drGodunov[face][faceRelation]);
        (*currentTable)[key].content[*EntityId::FluxSolver] = new BatchPointers(drFluxSolver[face][faceRelation]);
      }
    }
  }
}
