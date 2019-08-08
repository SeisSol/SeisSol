#include "Kernels/DirichletBoundary.h"
void seissol::kernels::computeAverageDisplacement(double deltaT,
					 const real* timeDerivatives,
					 const unsigned int derivativesOffsets[CONVERGENCE_ORDER],
					 real timeIntegrated[tensor::I::size()] 
					 ) {
  // TODO(Lukas) Verify & optimize this.
  // TODO(Lukas) Only compute integral for displacement, not for all vars.
  assert((uintptr_t)(timeDerivatives) % ALIGNMENT == 0);
  assert((uintptr_t)(timeIntegrated) % ALIGNMENT == 0);
  assert(deltaT > 0);
  
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + derivativesOffsets[i];
  }
  
  auto factorial = static_cast<real>(2.0 * 1.0);
  auto power = deltaT * deltaT;
  
  for (int der = 0; der < CONVERGENCE_ORDER; ++der) {
    intKrnl.power = power / factorial;
    intKrnl.execute(der);
    
    factorial *= der + 1;
    power *= deltaT;
  }
}

