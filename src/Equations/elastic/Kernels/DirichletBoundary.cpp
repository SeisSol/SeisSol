// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "Kernels/DirichletBoundary.h"
#pragma GCC diagnostic pop

namespace seissol::kernels {

void computeAverageDisplacement(double deltaT,
					 const real* timeDerivatives,
					 const unsigned int derivativesOffsets[ConvergenceOrder],
					 real timeIntegrated[tensor::I::size()] 
					 ) {
  // TODO(Lukas) Only compute integral for displacement, not for all vars.
  assert(reinterpret_cast<uintptr_t>(timeDerivatives) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(timeIntegrated) % Alignment == 0);
  assert(deltaT > 0);
  
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeIntegrated;
  for (size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + derivativesOffsets[i];
  }
  
  real factorial = 2.0;
  double power = deltaT * deltaT;
  
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    intKrnl.power(der) = power / factorial;

    factorial *= der + 2.0;
    power *= deltaT;
  }
  intKrnl.execute();
}

} // namespace seissol::kernels

