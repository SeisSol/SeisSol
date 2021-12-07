#include <Kernels/DeviceAux/PlasticityAux.h>
#include <init.h>
#include <cmath>
#include <type_traits>

// NOTE: using c++14 because of cuda@10
namespace seissol {
namespace kernels {
namespace device {
namespace aux {
namespace plasticity {

template<typename T>
__forceinline__ __device__ typename std::enable_if<std::is_floating_point<T>::value, T>::type
squareRoot(T x) {
  return std::is_same<T, double>::value ? sqrt(x) : sqrtf(x);
}

template<typename T>
__forceinline__ __device__ typename std::enable_if<std::is_floating_point<T>::value, T>::type
maxValue(T x, T y) {
  return std::is_same<T, double>::value ? fmax(x, y) : fmaxf(x, y);
}

template<typename Tensor>
__forceinline__  __device__
constexpr size_t leadDim() {
  return Tensor::Stop[0] - Tensor::Start[0];
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_saveFirstMode(real *firstModes,
                                     const real **modalStressTensors) {
  constexpr auto modalStressTensorsColumn = leadDim<init::Q>();
  firstModes[threadIdx.x + blockDim.x * blockIdx.x] =
      modalStressTensors[blockIdx.x][threadIdx.x * modalStressTensorsColumn];
}

void saveFirstModes(real *firstModes,
                    const real **modalStressTensors,
                    const size_t numElements,
                    void *streamPtr) {
  dim3 block(NUM_STRESS_COMPONENTS, 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_saveFirstMode<<<grid, block, 0, stream>>>(firstModes, modalStressTensors);
}


//--------------------------------------------------------------------------------------------------
__global__ void kernel_adjustDeviatoricTensors(real **nodalStressTensors,
                                               int *isAdjustableVector,
                                               const PlasticityData *plasticity,
                                               const double oneMinusIntegratingFactor) {
  real *elementTensors = nodalStressTensors[blockIdx.x];
  real localStresses[NUM_STRESS_COMPONENTS];


  constexpr auto elementTensorsColumn = leadDim<init::QStressNodal>();
  #pragma unroll
  for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
    localStresses[i] = elementTensors[threadIdx.x + elementTensorsColumn * i];
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
  const real cohesionTimesCosAngularFriction = plasticity[blockIdx.x].cohesionTimesCosAngularFriction;
  const real sinAngularFriction = plasticity[blockIdx.x].sinAngularFriction;
  real taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
  taulim = maxValue(static_cast<real>(0.0), taulim);

  __shared__ int isAdjusted;
  if (threadIdx.x == 0) { isAdjusted = static_cast<int>(false); }
  __syncthreads();

  // 6. Compute the yield factor
  real factor = 0.0;
  if (tau > taulim) {
    isAdjusted = static_cast<int>(true);
    factor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
  }

  // 7. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
  __syncthreads();
  if (isAdjusted) {
    #pragma unroll
    for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
      elementTensors[threadIdx.x + elementTensorsColumn * i] = localStresses[i] * factor;
    }
  }

  if (threadIdx.x == 0) {
    isAdjustableVector[blockIdx.x] = isAdjusted;
  }
}

void adjustDeviatoricTensors(real **nodalStressTensors,
                             int *isAdjustableVector,
                             const PlasticityData *plasticity,
                             const double oneMinusIntegratingFactor,
                             const size_t numElements,
                             void *streamPtr) {
  constexpr unsigned numNodes = tensor::QStressNodal::Shape[0];
  dim3 block(numNodes, 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_adjustDeviatoricTensors<<<grid, block, 0, stream>>>(nodalStressTensors,
                                                             isAdjustableVector,
                                                             plasticity,
                                                             oneMinusIntegratingFactor);
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_adjustModalStresses(real **modalStressTensors,
                                           const real **nodalStressTensors,
                                           const real *inverseVandermondeMatrix,
                                           const int *isAdjustableVector) {
  if (isAdjustableVector[blockIdx.x]) {

    constexpr int numNodes = init::QStressNodal::Shape[0];
    __shared__ real shrMem[NUM_STRESS_COMPONENTS][numNodes];

    real *modalTensor = modalStressTensors[blockIdx.x];
    const real *nodalTensor = nodalStressTensors[blockIdx.x];

    constexpr size_t nodalTensorColumn = leadDim<init::QStressNodal>();
    for (int n = 0; n < NUM_STRESS_COMPONENTS; ++n) {
      shrMem[n][threadIdx.x] = nodalTensor[threadIdx.x + nodalTensorColumn * n];
    }
    __syncthreads();

    real accumulator[NUM_STRESS_COMPONENTS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    constexpr auto vInvColumn = leadDim<init::vInv>();
    for (int k = 0; k < numNodes; ++k) {
      real value = inverseVandermondeMatrix[threadIdx.x + vInvColumn * k];

      #pragma unroll
      for (int n = 0; n < NUM_STRESS_COMPONENTS; ++n) {
        accumulator[n] += value * shrMem[n][k];
      }
    }

    constexpr size_t modalTensorColumn = leadDim<init::Q>();
    #pragma unroll
    for (int n = 0; n < NUM_STRESS_COMPONENTS; ++n) {
      modalTensor[threadIdx.x + modalTensorColumn * n] += accumulator[n];
    }
  }
}

void adjustModalStresses(real **modalStressTensors,
                         const real **nodalStressTensors,
                         const real *inverseVandermondeMatrix,
                         const int *isAdjustableVector,
                         const size_t numElements,
                         void *streamPtr) {
  constexpr unsigned numNodesPerElement = init::QStressNodal::Shape[0];
  dim3 block(numNodesPerElement, 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);

  kernel_adjustModalStresses<<<grid, block, 0, stream>>>(modalStressTensors,
                                                         nodalStressTensors,
                                                         inverseVandermondeMatrix,
                                                         isAdjustableVector);
}

//--------------------------------------------------------------------------------------------------
__global__ void kernel_computePstrains(real **pstrains,
                                       const int *isAdjustableVector,
                                       const real **modalStressTensors,
                                       const real *firsModes,
                                       const PlasticityData *plasticity,
                                       const double oneMinusIntegratingFactor,
                                       const double timeStepWidth,
                                       const double T_v,
                                       const size_t numElements) {
  // compute element id
  size_t index = threadIdx.y + blockIdx.x * blockDim.y;
  if ((isAdjustableVector[index]) && (index < numElements)) {
    // NOTE: Six threads (x-dimension) work on the same element.

    // get local data
    real *localPstrains = pstrains[index];
    const real *localModalTensor = modalStressTensors[index];
    const real *localFirstMode = &firsModes[NUM_STRESS_COMPONENTS * index];
    const PlasticityData *localData = &plasticity[index];

    constexpr auto elementTensorsColumn = leadDim<init::QStressNodal>();
    real factor = localData->mufactor / (T_v * oneMinusIntegratingFactor);
    real duDtPstrain = factor * (localFirstMode[threadIdx.x] - localModalTensor[threadIdx.x * elementTensorsColumn]);
    localPstrains[threadIdx.x] += timeStepWidth * duDtPstrain;

    __shared__ real squaredDuDtPstrains[NUM_STRESS_COMPONENTS];
    real coefficient = threadIdx.x < 3 ? static_cast<real>(0.5) : static_cast<real>(1.0);
    squaredDuDtPstrains[threadIdx.x] = coefficient * duDtPstrain * duDtPstrain;
    __syncthreads();

    if (threadIdx.x == 0) {
      real sum = 0.0;

      #pragma unroll
      for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
        sum += squaredDuDtPstrains[i];
      }
      localPstrains[6] += (timeStepWidth * squareRoot(duDtPstrain));
    }
  }
}


void computePstrains(real **pstrains,
                     const int *isAdjustableVector,
                     const real **modalStressTensors,
                     const real *firsModes,
                     const PlasticityData *plasticity,
                     const double oneMinusIntegratingFactor,
                     const double timeStepWidth,
                     const double T_v,
                     const size_t numElements,
                     void *streamPtr) {
  dim3 block(NUM_STRESS_COMPONENTS, 32, 1);
  size_t numBlocks = (numElements + block.y - 1) / block.y;
  dim3 grid(numBlocks, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(streamPtr);
  kernel_computePstrains<<<grid, block, 0, stream>>>(pstrains,
                                                     isAdjustableVector,
                                                     modalStressTensors,
                                                     firsModes,
                                                     plasticity,
                                                     oneMinusIntegratingFactor,
                                                     timeStepWidth,
                                                     T_v,
                                                     numElements);
}
} // namespace plasticity
} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol
