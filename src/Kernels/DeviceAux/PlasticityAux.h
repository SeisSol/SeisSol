#ifndef SEISSOL_DEVICEAUX_PLASTICITY_H
#define SEISSOL_DEVICEAUX_PLASTICITY_H

#include "Initializer/BasicTypedefs.hpp"
#include <stddef.h>

#define NUM_STRESS_COMPONENTS 6

// NOTE: using c++14 because of cuda@10
namespace seissol {
namespace kernels {
namespace device {
namespace aux {
namespace plasticity {
void adjustDeviatoricTensors(real **nodalStressTensors,
                             unsigned *isAdjustableVector,
                             const PlasticityData *plasticity,
                             double oneMinusIntegratingFactor,
                             size_t numElements,
                             void* streamPtr);

void adjustPointers(real *QEtaNodal,
                    real **QEtaNodalPtrs,
                    real *QEtaModal,
                    real **QEtaModalPtrs,
                    real *dUdTpstrain,
                    real **dUdTpstrainPtrs,
                    size_t numElements,
                    void *streamPtr);

void computePstrains(real **pstrains,
                     const PlasticityData *plasticityData,
                     real **dofs,
                     real *prevDofs,
                     real **dUdTpstrain,
                     double T_v,
                     double oneMinusIntegratingFactor,
                     double timeStepWidth,
                     unsigned *isAdjustableVector,
                     size_t numElements,
                     void *streamPtr);

void pstrainToQEtaModal(real **pstrains,
                        real **QEtaModalPtrs,
                        unsigned *isAdjustableVector,
                        size_t numElements,
                        void *streamPtr);

void qEtaModalToPstrain(real **QEtaModalPtrs,
                        real **pstrains,
                        unsigned *isAdjustableVector,
                        size_t numElements,
                        void *streamPtr);

void updateQEtaNodal(real **QEtaNodalPtrs,
                     real **QStressNodalPtrs,
                     double timeStepWidth,
                     unsigned *isAdjustableVector,
                     size_t numElements,
                     void *streamPtr);

} // namespace plasticity
} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol


#endif // SEISSOL_DEVICEAUX_PLASTICITY_H
