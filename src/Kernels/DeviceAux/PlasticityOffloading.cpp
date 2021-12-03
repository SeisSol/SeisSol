#include <Kernels/precision.hpp>
#include <Kernels/Plasticity.h>

#include <cstring>
#include <algorithm>
#include <cmath>
#include <generated_code/kernel.h>
#include <generated_code/init.h>
#include <Kernels/common.hpp>

/*
#ifdef ACL_DEVICE
#include "device.h"
#include "DeviceAux/PlasticityAux.h"
using namespace device;
#endif
*/

#include <omp.h>
#include <iostream>

namespace seissol::kernels {
void testOmpOffloading() {
  bool canOffload = false;
#pragma omp target map(tofrom: canOffload)
  {
    if (!omp_is_initial_device()) {
      canOffload = true;
    }
  }

  if (canOffload) {
    std::cout << "Device can perform OMP Offloading" << std::endl;
  } else {
    std::cout << "Device cannot perform OMP Offloading" << std::endl;
  }
}



unsigned Plasticity::computePlasticityBatched(double oneMinusIntegratingFactor,
                                              double timeStepWidth,
                                              double T_v,
                                              GlobalData const *global,
                                              initializers::recording::ConditionalBatchTableT &table,
                                              PlasticityData *plasticity) {

#ifdef ACL_DEVICE
  seissol::kernels::testOmpOffloading();
  return 0;
  /*
  static_assert(tensor::Q::Shape[0] == tensor::QStressNodal::Shape[0],
                  "modal and nodal dofs must have the same leading dimensions");
    static_assert(tensor::Q::Shape[0] == tensor::v::Shape[0],
                  "modal dofs and vandermonde matrix must hage the same leading dimensions");

    DeviceInstance &device = DeviceInstance::getInstance();
    ConditionalKey key(*KernelNames::Plasticity);
    auto defaultStream = device.api->getDefaultStream();

    if (table.find(key) != table.end()) {
      unsigned stackMemCounter{0};
      BatchTable& entry = table[key];
      const size_t numElements = (entry.content[*EntityId::Dofs])->getSize();
      real** modalStressTensors = (entry.content[*EntityId::Dofs])->getPointers();
      real** nodalStressTensors = (entry.content[*EntityId::NodalStressTensor])->getPointers();

      const size_t firsModesSize = NUM_STRESS_COMPONENTS * numElements * sizeof(real);
      real *firsModes = reinterpret_cast<real*>(device.api->getStackMemory(firsModesSize));
      ++stackMemCounter;

      device::aux::plasticity::saveFirstModes(firsModes,
                                              const_cast<const real**>(modalStressTensors),
                                              numElements,
                                              defaultStream);

      auto defaultStream = device.api->getDefaultStream();

      assert(global->replicateStresses != nullptr && "replicateStresses has not been initialized");
      real** initLoad = (entry.content[*EntityId::InitialLoad])->getPointers();
      kernel::gpu_plConvertToNodal m2nKrnl;
      m2nKrnl.v = global->vandermondeMatrix;
      m2nKrnl.QStress = const_cast<const real**>(modalStressTensors);
      m2nKrnl.QStressNodal = nodalStressTensors;
      m2nKrnl.replicateInitialLoadingM = global->replicateStresses;
      m2nKrnl.initialLoadingM = const_cast<const real**>(initLoad);
      m2nKrnl.streamPtr = defaultStream;
      m2nKrnl.numElements = numElements;

      const size_t MAX_TMP_MEM = m2nKrnl.TmpMaxMemRequiredInBytes * numElements;
      real *tmpMem = (real*)(device.api->getStackMemory(MAX_TMP_MEM));
      ++stackMemCounter;

      m2nKrnl.linearAllocator.initialize(tmpMem);
      m2nKrnl.execute();
      m2nKrnl.linearAllocator.free();

      int *isAdjustableVector =
          reinterpret_cast<int*>(device.api->getStackMemory(numElements * sizeof(int)));
      ++stackMemCounter;

      device::aux::plasticity::adjustDeviatoricTensors(nodalStressTensors,
                                                       isAdjustableVector,
                                                       plasticity,
                                                       oneMinusIntegratingFactor,
                                                       numElements,
                                                       defaultStream);

      unsigned numAdjustedDofs = device.algorithms.reduceVector(isAdjustableVector,
                                                                numElements,
                                                                ::device::ReductionType::Add,
                                                                defaultStream);

      // apply stress adjustment
      device::aux::plasticity::adjustModalStresses(modalStressTensors,
                                                   const_cast<const real **>(nodalStressTensors),
                                                   global->vandermondeMatrixInverse,
                                                   isAdjustableVector,
                                                   numElements,
                                                   defaultStream);

      // compute Pstrains
      real **pstrains = entry.content[*EntityId::Pstrains]->getPointers();
      device::aux::plasticity::computePstrains(pstrains,
                                               isAdjustableVector,
                                               const_cast<const real **>(modalStressTensors),
                                               firsModes,
                                               plasticity,
                                               oneMinusIntegratingFactor,
                                               timeStepWidth,
                                               T_v,
                                               numElements,
                                               defaultStream);

      // NOTE: Temp memory must be properly clean after using negative signed integers
      // This kind of memory is mainly used for floating-point numbers. Negative signed ints might corrupt
      // the most significant bits. We came to this conclusion by our first-hand experience
      device.algorithms.fillArray(reinterpret_cast<char*>(isAdjustableVector),
                                  static_cast<char>(0),
                                  numElements * sizeof(int),
                                  defaultStream);

      for (unsigned i = 0; i < stackMemCounter; ++i) {
        device.api->popStackMemory();
      }
      return numAdjustedDofs;
    }
    else {
      return 0;
    }
  */
#else
  assert(false && "no implementation provided");
  return 0;
#endif // ACL_DEVICE
}
}  // namespace seissol::kernels