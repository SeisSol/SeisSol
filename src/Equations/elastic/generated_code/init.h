#ifndef GENERATED_INIT_H_
#define GENERATED_INIT_H_

#include <Numerical_aux/MatrixView.h>
#include <generated_code/initialization/bind.h>

namespace seissol {
  namespace model {
    namespace AminusT {
      unsigned const reals = NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES;
      int (* const index)(unsigned, unsigned) = &colMjrIndex<NUMBER_OF_QUANTITIES>;
    }
    namespace AplusT {
      unsigned const reals = NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES;
      int (* const index)(unsigned, unsigned) = &colMjrIndex<NUMBER_OF_QUANTITIES>;
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
