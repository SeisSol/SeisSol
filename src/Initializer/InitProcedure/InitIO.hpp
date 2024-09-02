#ifndef INITPROC_IO_H
#define INITPROC_IO_H

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initIO(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif
