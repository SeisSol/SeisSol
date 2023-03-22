
#ifndef INIT_HPP_
#define INIT_HPP_

#include <vector>

#include "Initializer/tree/Lut.hpp"
#include "Initializer/typedefs.hpp"

namespace seissol::initializer::initprocedure {

struct LtsInfo {
    unsigned* ltsMeshToFace = nullptr;
    MeshStructure* meshStructure = nullptr;
    TimeStepping timeStepping;

    ~LtsInfo() {
        // TODO: refactor LtsLayout, so that the following checks can be removed entirely.
        if (ltsMeshToFace != nullptr) {
            delete[] ltsMeshToFace;
            ltsMeshToFace = nullptr;
        }
        if (meshStructure != nullptr) {
            delete[] meshStructure;
            meshStructure = nullptr;
        }
    }
};

void seissolMain();

}

#endif
