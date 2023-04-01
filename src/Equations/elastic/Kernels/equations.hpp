#ifndef SEISSOL_EQUATIONS
#define SEISSOL_EQUATIONS


constexpr unsigned numberOfBasisFunctions(unsigned O) {
  return O * (O + 1) * (O + 2) / 6;
}

constexpr unsigned numberOfAlignedBasisFunctions(unsigned O) {
  return (numberOfBasisFunctions(O) * sizeof(real) + (ALIGNMENT - (numberOfBasisFunctions(O) * sizeof(real)) % ALIGNMENT) % ALIGNMENT) / sizeof(real);
}

constexpr unsigned numberOfAlignedDerBasisFunctions(unsigned O) {
  return (O > 0) ? numberOfAlignedBasisFunctions(O) + numberOfAlignedDerBasisFunctions(O-1) : 0;
}

#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS numberOfAlignedBasisFunctions(CONVERGENCE_ORDER)
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS numberOfAlignedDerBasisFunctions(CONVERGENCE_ORDER)

#endif  // SEISSOL_EQUATIONS

