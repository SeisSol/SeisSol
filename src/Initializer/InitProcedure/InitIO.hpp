#ifndef INITPROC_IO_H
#define INITPROC_IO_H

#include "Initializer/InitProcedure/Init.hpp"

namespace seissol {
  class SeisSol;
  namespace initializers::initprocedure {
void initIO(seissol::SeisSol& seissolInstance);
  }
} // namespace seissol::initializers::initprocedure

#endif
