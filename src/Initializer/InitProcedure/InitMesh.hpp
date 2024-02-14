#ifndef INITPROC_MESH_H
#define INITPROC_MESH_H

#include "Initializer/InitProcedure/Init.hpp"

namespace seissol {
class SeisSol;
namespace initializer::initprocedure {
void initMesh(seissol::SeisSol& seissolInstance);
}
} // namespace seissol

#endif
