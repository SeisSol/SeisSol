#ifndef SEISSOL_DEVICEAUX_PLASTICITY_H
#define SEISSOL_DEVICEAUX_PLASTICITY_H

#include <Initializer/BasicTypedefs.hpp>

#define NUM_STREESS_COMPONENTS 6

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
} // namespace algorithms
} // namespace seissol

#endif // SEISSOL_DEVICEAUX_PLASTICITY_H
