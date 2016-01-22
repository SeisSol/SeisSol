#ifndef GENERATED_SIZES_H_
#define GENERATED_SIZES_H_

#include <generated_code/initialization/bind.h>

namespace seissol {
  namespace model {
    namespace AminusT {
      unsigned const rows = NUMBER_OF_QUANTITIES;
      unsigned const cols = NUMBER_OF_QUANTITIES;
      unsigned const reals = NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES;
    }
    namespace AplusT {
      unsigned const rows = NUMBER_OF_QUANTITIES;
      unsigned const cols = NUMBER_OF_QUANTITIES;
      unsigned const reals = NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES;
    }
    namespace AstarT {
      unsigned const reals = STAR_NNZ;
    }
    namespace BstarT {
      unsigned const reals = STAR_NNZ;
    }
    namespace CstarT {
      unsigned const reals = STAR_NNZ;
    }
  }
}
#endif
