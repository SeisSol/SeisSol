#ifndef MULTIPLE_SIMULATIONS_H
#define MULTIPLE_SIMULATIONS_H

namespace seissol {
#ifdef MULTIPLE_SIMULATIONS
  constexpr unsigned int numberOfMultipleSimulations = MULTIPLE_SIMULATIONS;
#else
  constexpr unsigned int numberOfMultipleSimulations = 1;
#endif
}

#endif
