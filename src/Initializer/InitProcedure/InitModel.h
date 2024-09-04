#ifndef INITPROC_CELLS_H
#define INITPROC_CELLS_H

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initModel(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif
