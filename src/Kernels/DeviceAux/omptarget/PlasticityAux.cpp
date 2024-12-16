#include <Kernels/DeviceAux/PlasticityAux.h>
#include <init.h>
#include <cmath>
#include <omp.h>


namespace seissol::kernels::device::aux::plasticity {

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
squareRoot(T x) {
  return std::is_same<T, double>::value ? sqrt(x) : sqrtf(x);
}

template<typename Tensor>
constexpr size_t leadDim() {
  return Tensor::Stop[0] - Tensor::Start[0];
}


void adjustDeviatoricTensors(real **nodalStressTensors,
                             unsigned *isAdjustableVector,
                             const model::PlasticityData *plasticity,
                             const double oneMinusIntegratingFactor,
                             const size_t numElements,
                             void *queuePtr) {
  constexpr unsigned numNodes = tensor::QStressNodal::Shape[0];
  auto stream = reinterpret_cast<int*>(queuePtr);

  #pragma omp target teams nowait depend(inout: stream[0]) num_teams(numElements)
  {
    int isAdjusted[1];
    #pragma omp parallel num_threads(numNodes)
    {
        auto wid = omp_get_team_num();
        auto tid = omp_get_thread_num();

        real *elementTensors = nodalStressTensors[wid];
      real localStresses[NUM_STRESS_COMPONENTS];


      constexpr auto elementTensorsColumn = leadDim<init::QStressNodal>();
      #pragma unroll
      for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
        localStresses[i] = elementTensors[tid + elementTensorsColumn * i];
      }

      // 1. Compute the mean stress for each node
      real meanStress = (localStresses[0] + localStresses[1] + localStresses[2]) / 3.0f;

      // 2. Compute deviatoric stress tensor
      #pragma unroll
      for (int i = 0; i < 3; ++i) {
        localStresses[i] -= meanStress;
      }

      // 3. Compute the second invariant for each node
      real tau = 0.5 * (localStresses[0] * localStresses[0] +
                        localStresses[1] * localStresses[1] +
                        localStresses[2] * localStresses[2]);
      tau += (localStresses[3] * localStresses[3] +
              localStresses[4] * localStresses[4] +
              localStresses[5] * localStresses[5]);
      tau = squareRoot(tau);

      // 4. Compute the plasticity criteria
      const real cohesionTimesCosAngularFriction = plasticity[wid].cohesionTimesCosAngularFriction;
      const real sinAngularFriction = plasticity[wid].sinAngularFriction;
      real taulim = cohesionTimesCosAngularFriction - meanStress * sinAngularFriction;
      taulim = std::max(static_cast<real>(0.0), taulim);

      if (tid == 0) { isAdjusted[0] = static_cast<unsigned>(false); }
      #pragma omp barrier

      // 5. Compute the yield factor
      real factor = 0.0;
      if (tau > taulim) {
        isAdjusted[0] = static_cast<unsigned>(true);
        factor = ((taulim / tau) - 1.0) * oneMinusIntegratingFactor;
      }

      // 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
      #pragma omp barrier
      if (isAdjusted[0]) {
      #pragma unroll
        for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
          elementTensors[tid + elementTensorsColumn * i] = localStresses[i] * factor;
        }
      }
      if (tid == 0) {
        isAdjustableVector[wid] = isAdjusted[0];
      }
    }
  }
}


void adjustPointers(real *QEtaNodal,
                    real **QEtaNodalPtrs,
                    real *QEtaModal,
                    real **QEtaModalPtrs,
                    real *dUdTpstrain,
                    real **dUdTpstrainPtrs,
                    size_t numElements,
                    void *queuePtr) {
  auto stream = reinterpret_cast<int*>(queuePtr);

  #pragma omp target teams distribute parallel for depend(inout: stream[0]) nowait
  for (size_t tid = 0; tid < numElements; ++tid) {
    QEtaNodalPtrs[tid] = &QEtaNodal[tensor::QEtaNodal::Size * tid];
    QEtaModalPtrs[tid] = &QEtaModal[tensor::QEtaModal::Size * tid];
    dUdTpstrainPtrs[tid] = &dUdTpstrain[tensor::QStressNodal::Size * tid];
  }
}


void computePstrains(real **pstrains,
                     const model::PlasticityData *plasticityData,
                     real **dofs,
                     real *prevDofs,
                     real **dUdTpstrain,
                     double T_v,
                     double oneMinusIntegratingFactor,
                     double timeStepWidth,
                     unsigned *isAdjustableVector,
                     size_t numElements,
                     void *queuePtr) {
  constexpr unsigned numNodes = tensor::Q::Shape[0];

  auto stream = reinterpret_cast<int*>(queuePtr);

  #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t wid = 0; wid < numElements; ++wid) {
    #pragma omp parallel for
    for (size_t lid = 0; lid < numNodes; ++lid) {
        if (isAdjustableVector[wid]) {
        real *localDofs = dofs[wid];
        real *localPrevDofs = &prevDofs[tensor::Q::Size * wid];
        const model::PlasticityData *localData = &plasticityData[wid];
        real *localPstrain = pstrains[wid];
        real *localDuDtPstrain = dUdTpstrain[wid];

        #pragma unroll
        for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
          int q = lid + i * leadDim<init::Q>();
          real factor = localData->mufactor / (T_v * oneMinusIntegratingFactor);
          real nodeDuDtPstrain = factor * (localPrevDofs[q] - localDofs[q]);

          static_assert(leadDim<init::QStress>() == leadDim<init::Q>());
          localPstrain[q] += timeStepWidth * nodeDuDtPstrain;
          localDuDtPstrain[q] = nodeDuDtPstrain;
        }
      }
    }
  }
}


void pstrainToQEtaModal(real **pstrains,
                       real **QEtaModalPtrs,
                       unsigned *isAdjustableVector,
                       size_t numElements,
                       void *queuePtr) {

  auto stream = reinterpret_cast<int*>(queuePtr);

  #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t wid = 0; wid < numElements; ++wid) {
    #pragma omp parallel for
    for (size_t lid = 0; lid < tensor::QEtaModal::Size; ++lid) {
    if (isAdjustableVector[wid]) {
        real *localQEtaModal = QEtaModalPtrs[wid];
        real *localPstrain = pstrains[wid];
        localQEtaModal[lid] = localPstrain[NUM_STRESS_COMPONENTS * leadDim<init::QStressNodal>() + lid];
      }
    }
  }
}


void qEtaModalToPstrain(real **QEtaModalPtrs,
                        real **pstrains,
                        unsigned *isAdjustableVector,
                        size_t numElements,
                        void *queuePtr) {
                            auto stream = reinterpret_cast<int*>(queuePtr);

  #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t wid = 0; wid < numElements; ++wid) {
    #pragma omp parallel for
    for (size_t lid = 0; lid < tensor::QEtaModal::Size; ++lid) {
   if (isAdjustableVector[wid]) {
        real *localQEtaModal = QEtaModalPtrs[wid];
        real *localPstrain = pstrains[wid];
        localPstrain[NUM_STRESS_COMPONENTS * leadDim<init::QStressNodal>() + lid] = localQEtaModal[lid];
      }
    }
  }
}


void updateQEtaNodal(real **QEtaNodalPtrs,
                     real **QStressNodalPtrs,
                     double timeStepWidth,
                     unsigned *isAdjustableVector,
                     size_t numElements,
                     void *queuePtr) {
                        auto stream = reinterpret_cast<int*>(queuePtr);
                        #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t wid = 0; wid < numElements; ++wid) {
    #pragma omp parallel for
    for (size_t lid = 0; lid < tensor::QStressNodal::Shape[0]; ++lid) {
   if (isAdjustableVector[wid]) {
        real *localQEtaNodal = QEtaNodalPtrs[wid];
        real *localQStressNodal = QStressNodalPtrs[wid];
        real factor{0.0};

        constexpr auto ld = leadDim<init::QStressNodal>();
        #pragma unroll
        for (int i = 0; i < NUM_STRESS_COMPONENTS; ++i) {
          factor += localQStressNodal[lid + i * ld] * localQStressNodal[lid + i * ld];
        }

        localQEtaNodal[lid] = std::max(static_cast<real>(0.0), localQEtaNodal[lid])
                              + timeStepWidth * squareRoot(static_cast<real>(0.5) * factor);
      }
    }
  }
}

} // namespace seissol::kernels::device::aux::plasticity
