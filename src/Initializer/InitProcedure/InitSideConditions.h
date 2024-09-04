#ifndef INITPROC_SIDECONDITIONS_H
#define INITPROC_SIDECONDITIONS_H

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initSideConditions(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif
