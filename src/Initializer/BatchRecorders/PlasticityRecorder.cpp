#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void PlasticityRecorder::record(LTS &handler, Layer &layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  auto& conditionalTable = layer.getCondBatchTable();
}