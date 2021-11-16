#ifndef SEISSOL_DEVICEAUX_PLASTICITY_H
#define SEISSOL_DEVICEAUX_PLASTICITY_H

#include <Initializer/BasicTypedefs.hpp>
#include <stddef.h>

#define NUM_STRESS_COMPONENTS 6

// NOTE: using c++14 because of cuda@10
namespace seissol {
namespace kernels {
namespace device {
namespace aux {
namespace plasticity {
void saveFirstModes(real *firstModes,
                    const real **modalStressTensors,
                    size_t numElements);

void adjustDeviatoricTensors(real **nodalStressTensors,
                             int *isAdjustableVector,
                             const PlasticityData *plasticity,
                             double oneMinusIntegratingFactor,
                             size_t numElements);

void adjustModalStresses(real **modalStressTensors,
                         const real **nodalStressTensors,
                         const real *inverseVandermondeMatrix,
                         const int *isAdjustableVector,
                         size_t numElements);

void computePstrains(real **pstrains,
                     const int *isAdjustableVector,
                     const real **modalStressTensors,
                     const real *firsModes,
                     const PlasticityData *plasticity,
                     double oneMinusIntegratingFactor,
                     double timeStepWidth,
                     double T_v,
                     size_t numElements);
} // namespace plasticity
} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol


#endif // SEISSOL_DEVICEAUX_PLASTICITY_H
