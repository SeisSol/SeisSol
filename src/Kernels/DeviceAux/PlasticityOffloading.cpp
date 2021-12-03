#include <Kernels/precision.hpp>
#include <Kernels/Plasticity.h>

#include <cstring>
#include <algorithm>
#include <cmath>
#include <generated_code/kernel.h>
#include <generated_code/init.h>
#include <Kernels/common.hpp>

#include "device.h"
#include <omp.h>

#include <iostream>  // TODO: remove



namespace seissol::kernels {
#pragma omp declare target
template<typename Tensor>
constexpr size_t leadDim() {
  return Tensor::Stop[0] - Tensor::Start[0];
}

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
squareRoot(T x) {
  return std::is_same<T, double>::value ? sqrt(x) : sqrtf(x);
}
#pragma omp end declare target

size_t align(size_t number, size_t min) {
    return min * ((number + min - 1) / min);
}



#define NUM_STRESS_COMPONENTS 6
using namespace device;

unsigned Plasticity::computePlasticityBatched(double oneMinusIntegratingFactor,
                                              double timeStepWidth,
                                              double T_v,
                                              GlobalData const *global,
                                              initializers::recording::ConditionalBatchTableT &table,
                                              PlasticityData *plasticity) {

#ifdef ACL_DEVICE

  static_assert(tensor::Q::Shape[0] == tensor::QStressNodal::Shape[0],
                "modal and nodal dofs must have the same leading dimensions");
  static_assert(tensor::Q::Shape[0] == tensor::v::Shape[0],
                "modal dofs and vandermonde matrix must hage the same leading dimensions");

  DeviceInstance &device = DeviceInstance::getInstance();
  ConditionalKey key(*KernelNames::Plasticity);
  if (table.find(key) != table.end()) {

      unsigned stackMemCounter{0};
      BatchTable& entry = table[key];
      const size_t numElements = (entry.content[*EntityId::Dofs])->getSize();

      auto laneSize = device.api->getLaneSize();
      auto dev = device.api->getDeviceId();
      
      // Note: there are some problems with unified memory and OpenMP
      // Thus, making a copy to a normal memory
      constexpr unsigned dofsSize = tensor::Q::Size;
      const size_t dofsPtrsCopySize = dofsSize * numElements * sizeof(real);
      real *dofsCopy = reinterpret_cast<real*>(device.api->getStackMemory(dofsPtrsCopySize));
      ++stackMemCounter;

      auto stream = device.api->getDefaultStream();
      real** dofsPtrs = (entry.content[*EntityId::Dofs])->getPointers();
      device.algorithms.copyScatterToUniform(dofsPtrs, dofsCopy, dofsSize, numElements, stream);
      device.api->fastStreamsSync();

      const size_t firsModesSize = NUM_STRESS_COMPONENTS * numElements * sizeof(real);
      real *firstModes = reinterpret_cast<real*>(device.api->getStackMemory(firsModesSize));
      ++stackMemCounter;


      #pragma omp target teams num_teams(numElements) \
        is_device_ptr(firstModes, dofsCopy) device(dev)
      #pragma omp distribute
      for (size_t element = 0; element < numElements; ++element) {
        #pragma omp parallel num_threads(laneSize) 
        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
          constexpr auto next = leadDim<init::Q>();
          firstModes[i + NUM_STRESS_COMPONENTS * element] = dofsCopy[i * next + element * dofsSize];
        }
      }

      real** nodalStressTensors = (entry.content[*EntityId::NodalStressTensor])->getPointers();
      assert(global->replicateStresses != nullptr && "replicateStresses has not been initialized");
      real** initLoad = (entry.content[*EntityId::InitialLoad])->getPointers();

      kernel::gpu_plConvertToNodal m2nKrnl;
      m2nKrnl.v = global->vandermondeMatrix;
      m2nKrnl.QStress = const_cast<const real**>(dofsPtrs);
      m2nKrnl.QStressNodal = nodalStressTensors;
      m2nKrnl.replicateInitialLoadingM = global->replicateStresses;
      m2nKrnl.initialLoadingM = const_cast<const real**>(initLoad);
      m2nKrnl.streamPtr = stream;
      m2nKrnl.numElements = numElements;

      const size_t MAX_TMP_MEM = m2nKrnl.TmpMaxMemRequiredInBytes * numElements;
      real *tmpMem = (real*)(device.api->getStackMemory(MAX_TMP_MEM));
      ++stackMemCounter;

      m2nKrnl.linearAllocator.initialize(tmpMem);
      m2nKrnl.execute();
      m2nKrnl.linearAllocator.free();
      device.api->fastStreamsSync();
      
      int *isAdjustableVector =
          reinterpret_cast<int*>(device.api->getStackMemory(numElements * sizeof(int)));
      ++stackMemCounter;

      {
        constexpr unsigned numNodesPerElement = tensor::QStressNodal::Shape[0];
        constexpr auto next = leadDim<init::QStressNodal>();

        #pragma omp target teams num_teams(numElements) \
          is_device_ptr(nodalStressTensors, plasticity, isAdjustableVector) device(dev)
        #pragma omp distribute
        for (size_t element = 0; element < numElements; ++element) {
         
          bool isAdjusted{false};
          real factors[numNodesPerElement];
          real stresses[numNodesPerElement][NUM_STRESS_COMPONENTS];

          #pragma omp parallel num_threads(align(numNodesPerElement, laneSize)) \
            shared(isAdjusted, factors, stresses)
          #pragma omp for schedule(static, 1) 
          for (size_t node = 0; node < numNodesPerElement; ++node) {

            real *localStresses = &stresses[node][0];
            real *elementTensors = nodalStressTensors[element];
              for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
                 localStresses[i] = elementTensors[node + next * i];
              }

            // 2. Compute the mean stress for each node
            real meanStress = localStresses[0] 
                            + localStresses[1] 
                            + localStresses[2];
            meanStress /= 3.0f;

            // 3. Compute deviatoric stress tensor
            localStresses[0] -= meanStress;
            localStresses[1] -= meanStress;
            localStresses[2] -= meanStress;

            // 4. Compute the second invariant for each node
            real tau = 0.5 * (localStresses[0] * localStresses[0] +
                              localStresses[1] * localStresses[1] +
                              localStresses[2] * localStresses[2]);

            tau += (localStresses[3] * localStresses[3] +
                    localStresses[4] * localStresses[4] +
                    localStresses[5] * localStresses[5]);
            tau = squareRoot(tau);

            // 5. Compute the plasticity criteria
            const real cohesionTimesCosAngularFriction 
                = plasticity[element].cohesionTimesCosAngularFriction;
            const real sinAngularFriction 
                = plasticity[element].sinAngularFriction;
            real taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
            taulim = taulim > static_cast<real>(0.0) ? taulim : static_cast<real>(0.0);

            factors[node] = 0.0;
            if (tau > taulim) {
              isAdjusted = true;
              factors[node] = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
            }
          }        
          #pragma omp barrier

          #pragma omp parallel num_threads(align(numNodesPerElement, laneSize)) \
            shared(isAdjusted, factors, stresses)
          #pragma omp for schedule(static, 1) 
          for (size_t node = 0; node < numNodesPerElement; ++node) {
            real *localStresses = &stresses[node][0];
            real *elementTensors = nodalStressTensors[element];

            if (isAdjusted) {
              for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
                elementTensors[node + next * i] = localStresses[i] * factors[node];
              }
            }

            if (node == 0) {
              isAdjustableVector[element] = static_cast<int>(isAdjusted);
            }
          }
        }
        
        real *vandermondeMatrixInverse = global->vandermondeMatrixInverse;

        #pragma omp target teams num_teams(numElements) \
          is_device_ptr(nodalStressTensors, isAdjustableVector) \
          is_device_ptr(dofsCopy, vandermondeMatrixInverse) device(dev)
        #pragma omp distribute
        for (size_t element = 0; element < numElements; ++element) {
         
          real shrMem[numNodesPerElement][NUM_STRESS_COMPONENTS];
          int toAdjust = isAdjustableVector[element];

          #pragma omp parallel num_threads(align(numNodesPerElement, laneSize)) \
            shared(shrMem) firstprivate(toAdjust)
          #pragma omp for schedule(static, 1) 
          for (size_t node = 0; node < numNodesPerElement; ++node) {
            if (toAdjust) {
              // load to shr. mem
              real *elementTensors = nodalStressTensors[element];
              for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
                shrMem[node][i] = elementTensors[node + next * i];
              }
            }
          }


          #pragma omp parallel num_threads(align(numNodesPerElement, laneSize)) \
            shared(shrMem) firstprivate(toAdjust)
          #pragma omp for schedule(static, 1) 
          for (size_t node = 0; node < numNodesPerElement; ++node) {
            constexpr real zero{0.0};
            real accumulator[NUM_STRESS_COMPONENTS] = {zero, zero, zero, zero, zero, zero};
            real value{zero};

            constexpr auto ivmColumn = leadDim<init::vInv>();
            if (toAdjust) {
              for (size_t k = 0; k < numNodesPerElement; ++k) {
                value = vandermondeMatrixInverse[node + ivmColumn * k];
                for (size_t i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
                  accumulator[i] += value * shrMem[k][i];
                }
              }
              constexpr auto dofsColumn = leadDim<init::Q>();
              for (size_t i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
                dofsCopy[node + dofsColumn * i] += accumulator[i];
              }
            }
          }
        }
      }

      device.algorithms.copyUniformToScatter(dofsCopy, dofsPtrs, dofsSize, numElements, stream);
      device.api->fastStreamsSync();

/*
      constexpr unsigned numNodesPerElement = tensor::QStressNodal::Shape[0];
      dim3 block(numNodesPerElement, 1, 1);
      dim3 grid(numElements, 1, 1);
      auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
      kernel_adjustDeviatoricTensors<<<grid, block, 0, stream>>>(nodalStressTensors,
                                                                 isAdjustableVector,
                                                                 plasticity,
                                                                 oneMinusIntegratingFactor);
*/
      for (unsigned i = 0; i < stackMemCounter; ++i) {
        device.api->popStackMemory();
      }
  }


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
