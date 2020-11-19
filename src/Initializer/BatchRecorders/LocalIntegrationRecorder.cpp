#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>
#include <iostream>

#include "DataTypes/Condition.hpp"
#include "DataTypes/ConditionalTable.hpp"
#include "DataTypes/EncodedConstants.hpp"
#include "DataTypes/BatchTable.hpp"

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void LocalIntegrationRecorder::record(LTS &handler, Layer &layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  loadContext(handler, layer, loader);

  recordTimeIntegral();
  recordVolumeIntegral();
  recordLocalFluxIntegral();
  recordDisplacements();
}


void LocalIntegrationRecorder::recordTimeIntegral() {
  std::vector<real *> dofsBatch{};
  dofsBatch.resize(currentLayer->getNumberOfCells(), nullptr);

  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  checkKey(*currentTable, key);

  (*currentTable)[key].content[*EntityId::Dofs] = new BatchPointers(dofsBatch);
}


void LocalIntegrationRecorder::recordVolumeIntegral() {

}


void LocalIntegrationRecorder::recordLocalFluxIntegral() {

}


void LocalIntegrationRecorder::recordDisplacements() {

}