#include <Kernels/DeviceAux/PlasticityAux.h>
#include <init.h>
#include <cmath>
#include <CL/sycl.hpp>


namespace seissol::kernels::device::aux::plasticity {

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value,T>::type
squareRoot(T x) {
  return std::is_same<T, double>::value ? sqrt(x) : sqrtf(x);
}

//--------------------------------------------------------------------------------------------------
void saveFirstModes(real *firstModes,
                    const real **modalStressTensors,
                    const size_t numElements,
                    void *queuePtr) {
  cl::sycl::range<3> groupCount(NUM_STRESS_COMPONENTS, 1, 1);
  cl::sycl::range<3> groupSize(numElements, 1, 1);

  assert(queuePtr != nullptr && "a pointer to a SYCL queue must be a valid one");
  auto queue = reinterpret_cast<cl::sycl::queue*>(queuePtr);

  queue->submit([&](cl::sycl::handler &cgh) {
    cgh.parallel_for(cl::sycl::nd_range < 3 >
                     {{groupCount.get(0) * groupSize.get(0), groupCount.get(1) * groupSize.get(1),
                       groupCount.get(2) * groupSize.get(2)}, groupSize}, [=](cl::sycl::nd_item<3> item) {
      constexpr unsigned numModesPerElement = tensor::Q::Shape[0];
      firstModes[item.get_local_id(0) + item.get_local_range(0) * item.get_group().get_id(0)] =
          modalStressTensors[item.get_group().get_id(0)][item.get_local_id(0) * numModesPerElement];
    });
  });
}


//--------------------------------------------------------------------------------------------------
void adjustDeviatoricTensors(real **nodalStressTensors,
                             int *isAdjustableVector,
                             const PlasticityData *plasticity,
                             const double oneMinusIntegratingFactor,
                             const size_t numElements,
                             void *queuePtr) {
  constexpr unsigned numNodesPerElement = tensor::QStressNodal::Shape[0];
  cl::sycl::range<3> groupCount(numNodesPerElement, 1, 1);
  cl::sycl::range<3> groupSize(numElements, 1, 1);

  assert(queuePtr != nullptr && "a pointer to a SYCL queue must be a valid one");
  auto queue = reinterpret_cast<cl::sycl::queue*>(queuePtr);

  queue->submit([&](cl::sycl::handler &cgh) {

    cl::sycl::accessor<int, 1, cl::sycl::access::mode::read_write, cl::sycl::access::target::local> isAdjusted(1, cgh);

    cgh.parallel_for(cl::sycl::nd_range < 3 >
                     {{groupCount.get(0) * groupSize.get(0), groupCount.get(1) * groupSize.get(1),
                       groupCount.get(2) * groupSize.get(2)}, groupSize}, [=](cl::sycl::nd_item<3> item) {
      real *elementTensors = nodalStressTensors[item.get_group().get_id(0)];
      real localStresses[NUM_STRESS_COMPONENTS];


      // NOTE: item.get_local_range(0) == tensor::QStressNodal::Shape[0] i.e., num nodes
      constexpr unsigned numNodesPerElement = tensor::QStressNodal::Shape[0];
      #pragma unroll
      for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
        localStresses[i] = elementTensors[item.get_local_id(0) + numNodesPerElement * i];
      }

      // 2. Compute the mean stress for each node
      real meanStress = (localStresses[0] + localStresses[1] + localStresses[2]) / 3.0f;

      // 3. Compute deviatoric stress tensor
      #pragma unroll
      for (int i = 0; i < 3; ++i) {
        localStresses[i] -= meanStress;
      }

      // 4. Compute the second invariant for each node
      real tau = 0.5 * (localStresses[0] * localStresses[0] +
                        localStresses[1] * localStresses[1] +
                        localStresses[2] * localStresses[2]);
      tau += (localStresses[3] * localStresses[3] +
              localStresses[4] * localStresses[4] +
              localStresses[5] * localStresses[5]);
      tau = squareRoot(tau);

      // 5. Compute the plasticity criteria
      const real cohesionTimesCosAngularFriction = plasticity[item.get_group().get_id(0)].cohesionTimesCosAngularFriction;
      const real sinAngularFriction = plasticity[item.get_group().get_id(0)].sinAngularFriction;
      real taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
      taulim = cl::sycl::fmax(static_cast<real>(0.0), taulim);

      if (item.get_local_id(0) == 0) { isAdjusted[0] = static_cast<int>(false); }
      item.barrier();

      // 6. Compute the yield factor
      real factor = 0.0;
      if (tau > taulim) {
        isAdjusted[0] = static_cast<int>(true);
        factor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
      }

      // 7. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
      item.barrier();
      if (isAdjusted[0]) {
        #pragma unroll
        for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
          elementTensors[item.get_local_id(0) + item.get_local_range(0) * i] = localStresses[i] * factor;
        }
      }

      if (item.get_local_id(0) == 0) {
        isAdjustableVector[item.get_group().get_id(0)] = isAdjusted[0];
      }
    });
  });
}

//--------------------------------------------------------------------------------------------------
void adjustModalStresses(real **modalStressTensors,
                         const real **nodalStressTensors,
                         const real *inverseVandermondeMatrix,
                         const int *isAdjustableVector,
                         const size_t numElements,
                         void *queuePtr) {
  constexpr unsigned numNodesPerElement = init::vInv::Shape[0];
  cl::sycl::range<3> groupCount(numNodesPerElement, 1, 1);
  cl::sycl::range<3> groupSize(numElements, 1, 1);

  assert(queuePtr != nullptr && "a pointer to a SYCL queue must be a valid one");
  auto queue = reinterpret_cast<cl::sycl::queue*>(queuePtr);

  // NOTE: item.get_local_range(0) == init::QStressNodal::Shape[0]
  constexpr int numNodes = init::QStressNodal::Shape[0];
  constexpr size_t nodalTensorSize = numNodes * NUM_STRESS_COMPONENTS;

  queue->submit([&](cl::sycl::handler &cgh) {
    cl::sycl::accessor<real, 1, cl::sycl::access::mode::read_write, cl::sycl::access::target::local> shrMem(
        nodalTensorSize, cgh);

    cgh.parallel_for(cl::sycl::nd_range < 3 >
                     {{groupCount.get(0) * groupSize.get(0), groupCount.get(1) * groupSize.get(1),
                       groupCount.get(2) * groupSize.get(2)}, groupSize}, [=](cl::sycl::nd_item<3> item) {
      if (isAdjustableVector[item.get_group().get_id(0)]) {


        real *modalTensor = modalStressTensors[item.get_group().get_id(0)];
        const real *nodalTensor = nodalStressTensors[item.get_group().get_id(0)];

        for (int n = 0; n < NUM_STRESS_COMPONENTS; ++n) {
          shrMem[item.get_local_id(0) + numNodes * n] = nodalTensor[item.get_local_id(0) + numNodes * n];
        }
        item.barrier();

        // matrix multiply: (numNodes x numNodes) * (numNodes x NUM_STRESS_COMPONENTS)
        real accumulator[NUM_STRESS_COMPONENTS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        real value = 0.0;
        // inverseVandermondeMatrix - square matrix
        for (int k = 0; k < numNodes; ++k) {
          value = inverseVandermondeMatrix[item.get_local_id(0) + numNodes * k];

          #pragma unroll
          for (int n = 0; n < NUM_STRESS_COMPONENTS; ++n) {
            accumulator[n] += value * shrMem[k + numNodes * n];
          }
        }

        constexpr unsigned numModesPerElement = init::Q::Shape[0];
        #pragma unroll
        for (int n = 0; n < NUM_STRESS_COMPONENTS; ++n) {
          modalTensor[item.get_local_id(0) + numModesPerElement * n] += accumulator[n];
        }
      }
    });
  });
}

//--------------------------------------------------------------------------------------------------
void computePstrains(real **pstrains,
                     const int *isAdjustableVector,
                     const real **modalStressTensors,
                     const real *firsModes,
                     const PlasticityData *plasticity,
                     const double oneMinusIntegratingFactor,
                     const double timeStepWidth,
                     const double T_v,
                     const size_t numElements,
                     void *queuePtr) {

  cl::sycl::range<3> groupCount(NUM_STRESS_COMPONENTS, 32, 1);
  size_t numBlocks = (numElements + groupCount.get(1) - 1) / groupCount.get(1);
  cl::sycl::range<3> groupSize(numBlocks, 1, 1);

  assert(queuePtr != nullptr && "a pointer to a SYCL queue must be a valid one");
  auto queue = reinterpret_cast<cl::sycl::queue*>(queuePtr);

  queue->submit([&](cl::sycl::handler &cgh) {
    cl::sycl::accessor<real, 1, cl::sycl::access::mode::read_write, cl::sycl::access::target::local> squaredDuDtPstrains(
        NUM_STRESS_COMPONENTS, cgh);

    cgh.parallel_for(cl::sycl::nd_range < 3 >
                     {{groupCount.get(0) * groupSize.get(0), groupCount.get(1) * groupSize.get(1),
                       groupCount.get(2) * groupSize.get(2)}, groupSize}, [=](cl::sycl::nd_item<3> item) {

      // compute element id
      size_t index = item.get_local_id(1) + item.get_group().get_id(0) * item.get_local_range(1);
      if ((isAdjustableVector[index]) && (index < numElements)) {
        // NOTE: Six threads (x-dimension) work on the same element.

        // get local data
        real *localPstrains = pstrains[index];
        const real *localModalTensor = modalStressTensors[index];
        const real *localFirstMode = &firsModes[NUM_STRESS_COMPONENTS * index];
        const PlasticityData *localData = &plasticity[index];

        constexpr unsigned numModesPerElement = init::Q::Shape[0];
        real factor = localData->mufactor / (T_v * oneMinusIntegratingFactor);
        real duDtPstrain = factor * (localFirstMode[item.get_local_id(0)] -
                                     localModalTensor[item.get_local_id(0) * numModesPerElement]);
        localPstrains[item.get_local_id(0)] += timeStepWidth * duDtPstrain;

        real coefficient = item.get_local_id(0) < 3 ? 0.5f : 1.0f;
        squaredDuDtPstrains[item.get_local_id(0)] = coefficient * duDtPstrain * duDtPstrain;
        item.barrier();

        if (item.get_local_id(0) == 0) {
          real sum = 0.0;

          #pragma unroll
          for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
            sum += squaredDuDtPstrains[i];
          }
          localPstrains[6] += (timeStepWidth * squareRoot(duDtPstrain));
        }
      }
    });
  });
}
} // namespace seissol::kernels::device::aux::plasticity
