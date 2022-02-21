#ifndef ENERGYOUTPUT_H
#define ENERGYOUTPUT_H

#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>
#include <Initializer/typedefs.hpp>

namespace seissol {
namespace writer {

void printPlasticMoment(MeshReader const& i_meshReader, seissol::initializers::LTSTree* i_ltsTree,
                        seissol::initializers::LTS* i_lts, seissol::initializers::Lut* i_ltsLut);
}
} // namespace seissol

#endif // ENERGYOUTPUT_H
