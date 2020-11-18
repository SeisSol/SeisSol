#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>
#include <iostream>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void LocalIntegrationRecorder::record(LTS &handler, Layer &layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  auto& conditionalTable = layer.getCondBatchTable();

  recordTimeIntegral();
  recordVolumeIntegral();
  recordLocalFluxIntegral();
}


void LocalIntegrationRecorder::recordTimeIntegral() {

}


void LocalIntegrationRecorder::recordVolumeIntegral() {

}


void LocalIntegrationRecorder::recordLocalFluxIntegral() {

}