// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_DEVICEAUX_KERNELSAUX_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_DEVICEAUX_KERNELSAUX_H_

#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"

namespace seissol::kernels::time::aux {
#ifdef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
void taylorSum(
    std::size_t count, real** target, const real** source, const real* coeffs, void* stream);
#endif
} // namespace seissol::kernels::time::aux

namespace seissol::kernels::local_flux::aux {
template <typename Cfg, typename Derived>
struct DirichletBoundaryAux {
  using real = Real<Cfg>;
  void evaluate(real** dofsFaceBoundaryNodalPtrs, size_t numElements, void* deviceStream) {
    static_cast<Derived*>(this)->dispatch(dofsFaceBoundaryNodalPtrs, numElements, deviceStream);
  }
};

template <typename Cfg>
struct FreeSurfaceGravity : public DirichletBoundaryAux<Cfg, FreeSurfaceGravity<Cfg>> {
  using real = Real<Cfg>;
  real** displacementDataPtrs{};
  double* rhos{};
  double g{};

  void dispatch(real** dofsFaceBoundaryNodalPtrs, size_t numElements, void* deviceStream);
};

template <typename Cfg>
struct EasiBoundary : public DirichletBoundaryAux<Cfg, EasiBoundary<Cfg>> {
  using real = Real<Cfg>;
  real** easiBoundaryMapPtrs{};
  real** easiBoundaryConstantPtrs{};

  void dispatch(real** dofsFaceBoundaryNodalPtrs, size_t numElements, void* deviceStream);
};

} // namespace seissol::kernels::local_flux::aux

namespace seissol::kernels::time::aux {
template <typename Cfg>
struct TimeAux {
  using real = Real<Cfg>;
  static void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                                      real** displacementToGlobalDataPtrs,
                                      real** TPtrs,
                                      real** TinvPtrs,
                                      size_t numElements,
                                      void* deviceStream);

  static void initializeTaylorSeriesForGravitationalBoundary(real** prevCoefficientsPtrs,
                                                             real** integratedDisplacementNodalPtrs,
                                                             real** rotatedFaceDisplacementPtrs,
                                                             double deltaTInt,
                                                             size_t numElements,
                                                             void* deviceStream);

  static void computeInvAcousticImpedance(
      double* invImpedances, double* rhos, double* lambdas, size_t numElements, void* deviceStream);

  static void updateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
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
};
} // namespace seissol::kernels::time::aux

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DEVICEAUX_KERNELSAUX_H_
