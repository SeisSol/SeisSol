#ifndef GENERATED_INIT_H_
#define GENERATED_INIT_H_

#include <Numerical_aux/MatrixView.h>
#include <generated_code/initialization/bind.h>

namespace seissol {
  namespace model {
    namespace AminusT {
      int (* const index)(unsigned, unsigned) = &colMjrIndex<NUMBER_OF_QUANTITIES>;
    }
    namespace AplusT {
      int (* const index)(unsigned, unsigned) = &colMjrIndex<NUMBER_OF_QUANTITIES>;
    }
  }
}
#endif
