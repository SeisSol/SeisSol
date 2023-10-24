#ifndef MULTIPLE_SIMULATIONS_H
#define MULTIPLE_SIMULATIONS_H

namespace seissol::multipleSimulations {
#ifdef MULTIPLE_SIMULATIONS
  constexpr unsigned int numberOfSimulations = MULTIPLE_SIMULATIONS;
  constexpr unsigned int basisFunctionDimension = 1;
#else
  constexpr unsigned int numberOfSimulations = 1;
  constexpr unsigned int basisFunctionDimension = 0;
#endif

}

#endif
