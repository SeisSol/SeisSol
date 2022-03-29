#ifndef ENERGYOUTPUT_H
#define ENERGYOUTPUT_H

#include <Initializer/typedefs.hpp>
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/Lut.hpp>

namespace seissol {
namespace writer {

real computePlasticMoment(MeshReader const& i_meshReader,
                          seissol::initializers::LTSTree* i_ltsTree,
                          seissol::initializers::LTS* i_lts,
                          seissol::initializers::Lut* i_ltsLut);

real computeStaticWork(GlobalData const* global,
                       real* degreesOfFreedomPlus,
                       real* degreesOfFreedomMinus,
                       DRFaceInformation const& faceInfo,
                       DRGodunovData const& godunovData,
                       real slip[seissol::tensor::slipInterpolated::size()]);

void printDynamicRuptureEnergies(GlobalData const* global,
                                 seissol::initializers::DynamicRupture* dynRup,
                                 seissol::initializers::LTSTree* dynRupTree,
                                 MeshReader const& i_meshReader,
                                 seissol::initializers::LTSTree* i_ltsTree,
                                 seissol::initializers::LTS* i_lts,
                                 seissol::initializers::Lut* i_ltsLut,
                                 bool usePlasticity);

void printEnergies(GlobalData const* global,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::initializers::LTSTree* dynRupTree,
                   MeshReader const& meshReader,
                   seissol::initializers::LTSTree* ltsTree,
                   seissol::initializers::LTS* lts,
                   seissol::initializers::Lut* ltsLut,
                   bool usePlasticity);

} // namespace writer
} // namespace seissol

#endif // ENERGYOUTPUT_H
