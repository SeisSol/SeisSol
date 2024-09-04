#ifndef INITPROC_MESH_H
#define INITPROC_MESH_H

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initMesh(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif
