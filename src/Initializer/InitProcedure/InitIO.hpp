#ifndef INITPROC_IO_H
#define INITPROC_IO_H

#include "Initializer/InitProcedure/Init.hpp"

namespace seissol {
class SeisSol;
namespace initializer::initprocedure {
void initIO(seissol::SeisSol& seissolInstance);
} // namespace initializer::initprocedure
} // namespace seissol

#endif
