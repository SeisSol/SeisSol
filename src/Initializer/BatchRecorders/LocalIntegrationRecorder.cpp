#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

#include "DataTypes/BatchTable.hpp"
#include "DataTypes/Condition.hpp"
#include "DataTypes/ConditionalTable.hpp"
#include "DataTypes/EncodedConstants.hpp"

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;


void LocalIntegrationRecorder::record(LTS &handler, Layer &layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  setUpContext(handler, layer, loader);
  idofsAddressRegistry.clear();

  recordTimeAndVolumeIntegrals();
  recordLocalFluxIntegral();
  recordDisplacements();
}


void LocalIntegrationRecorder::recordTimeAndVolumeIntegrals() {
  real *idofsScratch = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsScratch));
  real *derivativesScratch = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->derivativesScratch));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real *> dofsPtrs(size, nullptr);
    std::vector<real *> starPtrs(size, nullptr);
    std::vector<real *> idofsPtrs{};
    std::vector<real *> dQPtrs(size, nullptr);
    std::vector<real *> ltsBuffers{};
    std::vector<real *> idofsForLtsBuffers{};

    idofsPtrs.reserve(size);

    real **derivatives = currentLayer->var(currentHandler->derivatives);
    real **buffers = currentLayer->var(currentHandler->buffers);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);

      // dofs
      dofsPtrs[cell] = static_cast<real *>(data.dofs);

      // idofs
      real *nextIdofPtr = &idofsScratch[integratedDofsAddressCounter];
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
      starPtrs[cell] = static_cast<real *>(data.localIntegrationOnDevice.starMatrices[0]);

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

    (*currentTable)[key].content[*EntityId::Dofs] = new BatchPointers(dofsPtrs);
    (*currentTable)[key].content[*EntityId::Star] = new BatchPointers(starPtrs);
    (*currentTable)[key].content[*EntityId::Idofs] = new BatchPointers(idofsPtrs);
    (*currentTable)[key].content[*EntityId::Derivatives] = new BatchPointers(dQPtrs);


    if (!idofsForLtsBuffers.empty()) {
      ConditionalKey key(*KernelNames::Time, *ComputationKind::WithLtsBuffers);

      (*currentTable)[key].content[*EntityId::Buffers] = new BatchPointers(ltsBuffers);
      (*currentTable)[key].content[*EntityId::Idofs] = new BatchPointers(idofsForLtsBuffers);
    }
  }
}


void LocalIntegrationRecorder::recordLocalFluxIntegral() {
  const auto size = currentLayer->getNumberOfCells();
  for (unsigned face = 0; face < 4; ++face) {
    std::vector<real *> idofsPtrs{};
    std::vector<real *> dofsPtrs{};
    std::vector<real *> aplusTPtrs{};

    idofsPtrs.reserve(size);
    dofsPtrs.reserve(size);
    aplusTPtrs.reserve(size);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);

      // no element local contribution in the case of dynamic rupture boundary conditions
      if (data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {
        idofsPtrs.push_back(idofsAddressRegistry[cell]);
        dofsPtrs.push_back(static_cast<real *>(data.dofs));
        aplusTPtrs.push_back(static_cast<real *>(data.localIntegrationOnDevice.nApNm1[face]));
      }
    }

    // NOTE: we can check any container, but we must check that a set is not empty!
    if (!dofsPtrs.empty()) {
      ConditionalKey key(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);
      checkKey(key);
      (*currentTable)[key].content[*EntityId::Idofs] = new BatchPointers(idofsPtrs);
      (*currentTable)[key].content[*EntityId::Dofs] = new BatchPointers(dofsPtrs);
      (*currentTable)[key].content[*EntityId::AplusT] = new BatchPointers(aplusTPtrs);
    }
  }
}


void LocalIntegrationRecorder::recordDisplacements() {
  real **displacements = currentLayer->var(currentHandler->displacements);
  std::vector<real *> iVelocitiesPtrs{};
  std::vector<real *> displacementsPtrs{};

  // NOTE: velocity components are between 6th and 8th columns
  constexpr unsigned OFFSET_TO_VELOCITIES = 6 * seissol::tensor::I::Shape[0];
  const auto size = currentLayer->getNumberOfCells();
  for (unsigned cell = 0; cell < size; ++cell) {
    if (displacements[cell] != nullptr) {
      real *iVelocity = &idofsAddressRegistry[cell][OFFSET_TO_VELOCITIES];
      iVelocitiesPtrs.push_back(iVelocity);
      displacementsPtrs.push_back(displacements[cell]);
    }
  }
  if (!displacementsPtrs.empty()) {
    ConditionalKey key(*KernelNames::Displacements);
    checkKey(key);
    (*currentTable)[key].content[*EntityId::Ivelocities] = new BatchPointers(iVelocitiesPtrs);
    (*currentTable)[key].content[*EntityId::Displacements] = new BatchPointers(displacementsPtrs);
  }
}
