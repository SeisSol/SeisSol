#pragma once

#include "Kernels/precision.hpp"
#include "generated_code/init.h"

namespace seissol::kernels::local_flux::aux::details {
void launchFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                              real** displacementDataPtrs,
                              double* rhos,
                              double g,
                              size_t numElements,
                              void* deviceStream);

void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream);
} // namespace seissol::kernels::local_flux::aux::details

namespace seissol::kernels::local_flux::aux {
template <typename Derived>
struct DirichletBoundaryAux {
  void evaluate(real** dofsFaceBoundaryNodalPtrs, size_t numElements, void* deviceStream) {
    static_cast<Derived*>(this)->dispatch(dofsFaceBoundaryNodalPtrs, numElements, deviceStream);
  }
};

struct FreeSurfaceGravity : public DirichletBoundaryAux<FreeSurfaceGravity> {
  real** displacementDataPtrs{};
  double* rhos;
  double g{};

  void dispatch(real** dofsFaceBoundaryNodalPtrs, size_t numElements, void* deviceStream) {

    assert(displacementDataPtrs != nullptr);
    assert(rhos != nullptr);
    details::launchFreeSurfaceGravity(
        dofsFaceBoundaryNodalPtrs, displacementDataPtrs, rhos, g, numElements, deviceStream);
  }
};

struct EasiBoundary : public DirichletBoundaryAux<EasiBoundary> {
  real** easiBoundaryMapPtrs{};
  real** easiBoundaryConstantPtrs{};

  void dispatch(real** dofsFaceBoundaryNodalPtrs, size_t numElements, void* deviceStream) {

    assert(easiBoundaryMapPtrs != nullptr);
    assert(easiBoundaryConstantPtrs != nullptr);
    details::launchEasiBoundary(dofsFaceBoundaryNodalPtrs,
                                easiBoundaryMapPtrs,
                                easiBoundaryConstantPtrs,
                                numElements,
                                deviceStream);
  }
};

} // namespace seissol::kernels::local_flux::aux

namespace seissol::kernels::time::aux {
void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                             real** displacementToGlobalDataPtrs,
                             real** TPtrs,
                             real** TinvPtrs,
                             size_t numElements,
                             void* deviceStream);

void initializeTaylorSeriesForGravitationalBoundary(real** prevCoefficientsPtrs,
                                                    real** integratedDisplacementNodalPtrs,
                                                    real** rotatedFaceDisplacementPtrs,
                                                    double deltaTInt,
                                                    size_t numElements,
                                                    void* deviceStream);

void computeInvAcousticImpedance(
    double* invImpedances, double* rhos, double* lambdas, size_t numElements, void* deviceStream);

void updateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                   real** prevCoefficientsPtrs,
                                   real** integratedDisplacementNodalPtrs,
                                   real** dofsFaceNodalPtrs,
                                   double* invImpedances,
                                   double* rhos,
                                   double g,
                                   double factorEvaluated,
                                   double factorInt,
                                   size_t numElements,
                                   void* deviceStream);
} // namespace seissol::kernels::time::aux