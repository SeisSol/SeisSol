#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

#include "DataTypes/Table.hpp"
#include "DataTypes/Condition.hpp"
#include "DataTypes/ConditionalTable.hpp"
#include "DataTypes/EncodedConstants.hpp"

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void LocalIntegrationRecorder::record(LTS& handler, Layer& layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  setUpContext(handler, layer, loader);
  idofsAddressRegistry.clear();

  recordTimeAndVolumeIntegrals();
  recordFreeSurfaceBc();
  recordLocalFluxIntegral();
  recordDisplacements();
}

void LocalIntegrationRecorder::recordTimeAndVolumeIntegrals() {
  real* idofsScratch =
      static_cast<real*>(currentLayer->getScratchpadMemory(currentHandler->idofsScratch));
  real* derivativesScratch =
      static_cast<real*>(currentLayer->getScratchpadMemory(currentHandler->derivativesScratch));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real*> dofsPtrs(size, nullptr);
    std::vector<real*> starPtrs(size, nullptr);
    std::vector<real*> idofsPtrs{};
    std::vector<real*> dQPtrs(size, nullptr);
    std::vector<real*> ltsBuffers{};
    std::vector<real*> idofsForLtsBuffers{};

    idofsPtrs.reserve(size);

    real** derivatives = currentLayer->var(currentHandler->derivatives);
    real** buffers = currentLayer->var(currentHandler->buffers);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);

      // dofs
      dofsPtrs[cell] = static_cast<real*>(data.dofs);

      // idofs
      real* nextIdofPtr = &idofsScratch[integratedDofsAddressCounter];
      bool isBuffersProvided = ((data.cellInformation.ltsSetup >> 8) % 2) == 1;
      bool isLtsBuffers = ((data.cellInformation.ltsSetup >> 10) % 2) == 1;

      if (isBuffersProvided) {
        if (isLtsBuffers) {
          // lts buffers may require either accumulation or overriding (in case of reset command)
          idofsPtrs.push_back(nextIdofPtr);

          idofsForLtsBuffers.push_back(nextIdofPtr);
          ltsBuffers.push_back(buffers[cell]);

          idofsAddressRegistry[cell] = nextIdofPtr;
          integratedDofsAddressCounter += tensor::I::size();
        } else {
          // gts buffers have to be always overridden
          idofsPtrs.push_back(buffers[cell]);
          idofsAddressRegistry[cell] = buffers[cell];
        }
      } else {
        idofsPtrs.push_back(nextIdofPtr);
        idofsAddressRegistry[cell] = nextIdofPtr;
        integratedDofsAddressCounter += tensor::I::size();
      }

      // stars
      starPtrs[cell] = static_cast<real*>(data.localIntegrationOnDevice.starMatrices[0]);

      // derivatives
      bool isDerivativesProvided = ((data.cellInformation.ltsSetup >> 9) % 2) == 1;
      if (isDerivativesProvided) {
        dQPtrs[cell] = derivatives[cell];

      } else {
        dQPtrs[cell] = &derivativesScratch[derivativesAddressCounter];
        derivativesAddressCounter += yateto::computeFamilySize<tensor::dQ>();
      }
    }
    // just to be sure that we took all branches while filling in idofsPtrs vector
    assert(dofsPtrs.size() == idofsPtrs.size());

    ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    checkKey(key);

    (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Star, starPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Derivatives, dQPtrs);

    if (!idofsForLtsBuffers.empty()) {
      ConditionalKey key(*KernelNames::Time, *ComputationKind::WithLtsBuffers);

      (*currentTable)[key].set(inner_keys::Wp::Id::Buffers, ltsBuffers);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsForLtsBuffers);
    }
  }
}

void LocalIntegrationRecorder::recordLocalFluxIntegral() {
  const auto size = currentLayer->getNumberOfCells();
  for (unsigned face = 0; face < 4; ++face) {
    std::vector<real*> idofsPtrs{};
    std::vector<real*> dofsPtrs{};
    std::vector<real*> aplusTPtrs{};

    idofsPtrs.reserve(size);
    dofsPtrs.reserve(size);
    aplusTPtrs.reserve(size);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);

      // no element local contribution in the case of dynamic rupture boundary conditions
      if (data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {
        idofsPtrs.push_back(idofsAddressRegistry[cell]);
        dofsPtrs.push_back(static_cast<real*>(data.dofs));
        aplusTPtrs.push_back(static_cast<real*>(data.localIntegrationOnDevice.nApNm1[face]));
      }
    }

    // NOTE: we can check any container, but we must check that a set is not empty!
    if (!dofsPtrs.empty()) {
      ConditionalKey key(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::AplusT, aplusTPtrs);
    }
  }
}

void LocalIntegrationRecorder::recordDisplacements() {
  real*(*faceDisplacements)[4] = currentLayer->var(currentHandler->faceDisplacements);
  std::array<std::vector<real*>, 4> iVelocitiesPtrs{{}};
  std::array<std::vector<real*>, 4> displacementsPtrs{};

  const auto size = currentLayer->getNumberOfCells();
  for (unsigned cell = 0; cell < size; ++cell) {
    auto data = currentLoader->entry(cell);

    for (unsigned face = 0; face < 4; ++face) {
      auto isRequired = faceDisplacements[cell][face] != nullptr;
      auto notFreeSurfaceGravity =
          data.cellInformation.faceTypes[face] != FaceType::freeSurfaceGravity;

      if (isRequired && notFreeSurfaceGravity) {
        auto Iview = init::I::view::create(idofsAddressRegistry[cell]);
        // NOTE: velocity components are between 6th and 8th columns
        constexpr unsigned firstVelocityComponent{6};
        iVelocitiesPtrs[face].push_back(&Iview(0, firstVelocityComponent));
        displacementsPtrs[face].push_back(faceDisplacements[cell][face]);
      }
    }
  }

  for (unsigned face = 0; face < 4; ++face) {
    if (!displacementsPtrs[face].empty()) {
      ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Ivelocities, iVelocitiesPtrs[face]);
      (*currentTable)[key].set(inner_keys::Wp::Id::FaceDisplacement, displacementsPtrs[face]);
    }
  }
}

void LocalIntegrationRecorder::recordFreeSurfaceBc() {
  const auto size = currentLayer->getNumberOfCells();
  constexpr size_t nodalAvgDisplacementsSize = tensor::averageNormalDisplacement::size();

  real* nodalAvgDisplacements =
      static_cast<real*>(currentLayer->getScratchpadMemory(currentHandler->nodalAvgDisplacements));

  if (size > 0) {
    std::array<std::vector<unsigned>, 4> freeSurfaceGravityCells{};
    std::array<std::vector<real*>, 4> nodalAvgDisplacementsPtrs{};

    std::array<std::vector<unsigned>, 4> dirichletCells{};
    std::array<std::vector<unsigned>, 4> analiticalSolutionCells{};

    size_t nodalAvgDisplacementsCounter{0};

    kernels::LocalTmp tmp;
    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);

      for (unsigned face = 0; face < 4; ++face) {
        if (data.cellInformation.faceTypes[face] == FaceType::freeSurfaceGravity) {
          assert(data.faceDisplacements[face] != nullptr);
          freeSurfaceGravityCells[face].push_back(cell);

          real* displ{&nodalAvgDisplacements[nodalAvgDisplacementsCounter]};
          nodalAvgDisplacementsPtrs[face].push_back(displ);
          nodalAvgDisplacementsCounter += nodalAvgDisplacementsSize;
        }

        if (data.cellInformation.faceTypes[face] == FaceType::dirichlet) {
          dirichletCells[face].push_back(cell);
        }

        if (data.cellInformation.faceTypes[face] == FaceType::analytical) {
          analiticalSolutionCells[face].push_back(cell);
        }
      }
    }

    for (unsigned face = 0; face < 4; ++face) {
      if (!freeSurfaceGravityCells[face].empty()) {
        ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, face);
        checkKey(key);
        (*currentIndicesTable)[key].set(inner_keys::Indices::Id::Cells,
                                        freeSurfaceGravityCells[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::NodalAvgDisplacements,
                                 nodalAvgDisplacementsPtrs[face]);
      }

      if (!dirichletCells[face].empty()) {
        ConditionalKey key(*KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
        checkKey(key);
        (*currentIndicesTable)[key].set(inner_keys::Indices::Id::Cells, dirichletCells[face]);
      }

      if (!analiticalSolutionCells[face].empty()) {
        ConditionalKey key(*KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
        checkKey(key);
        (*currentIndicesTable)[key].set(inner_keys::Indices::Id::Cells,
                                        analiticalSolutionCells[face]);
      }
    }
  }
}
