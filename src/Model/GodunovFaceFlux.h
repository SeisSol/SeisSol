// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_GODUNOVFACEFLUX_H_
#define SEISSOL_SRC_MODEL_GODUNOVFACEFLUX_H_

#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/Precision.h"

namespace seissol::model::detail {

/// The pair of matrices a face applies, assembled from a Godunov state.
///
/// In a header of its own because the three solvers that assemble this way
/// share it word for word, and because the kernels it names are generated
/// only for a build whose flux is linear -- a solver that transports more
/// than its state has no Godunov state to assemble from.
inline void assembleGodunovFaceFlux(bool rusanov,
                                    double fluxScale,
                                    FaceType faceType,
                                    real* aPlusT,
                                    real* aMinusT,
                                    const real* centralFluxData,
                                    const real* rusanovPlusData,
                                    const real* rusanovMinusData,
                                    const real* qGodLocalData,
                                    const real* qGodNeighborData,
                                    const real* rusanovPlusNull,
                                    const real* rusanovMinusNull,
                                    const real* matTData,
                                    const real* matTinvData,
                                    const real* matATtildeData) {

  kernel::computeFluxSolverLocal localKrnl;
  localKrnl.fluxScale = fluxScale;
  localKrnl.AplusT = aPlusT;
  if (faceType == FaceType::DynamicRupture) {
    localKrnl.fluxScale = 0;
  }
  if (rusanov) {
    localKrnl.QgodLocal = centralFluxData;
    localKrnl.QcorrLocal = rusanovPlusData;
  } else {
    localKrnl.QgodLocal = qGodLocalData;
    localKrnl.QcorrLocal = rusanovPlusNull;
  }
  localKrnl.T = matTData;
  localKrnl.Tinv = matTinvData;
  localKrnl.star(0) = matATtildeData;
  localKrnl.execute();

  kernel::computeFluxSolverNeighbor neighKrnl;
  neighKrnl.fluxScale = fluxScale;
  neighKrnl.AminusT = aMinusT;
  if (rusanov) {
    neighKrnl.QgodNeighbor = centralFluxData;
    neighKrnl.QcorrNeighbor = rusanovMinusData;
  } else {
    neighKrnl.QgodNeighbor = qGodNeighborData;
    neighKrnl.QcorrNeighbor = rusanovMinusNull;
  }
  neighKrnl.T = matTData;
  neighKrnl.Tinv = matTinvData;
  neighKrnl.star(0) = matATtildeData;
  if (faceType == FaceType::Dirichlet || faceType == FaceType::FreeSurfaceGravity) {
    // already rotated
    neighKrnl.Tinv = init::identityT::Values;
  }
  neighKrnl.execute();
}

/// The boundary conditions a face folds into the pair it applies.
///
/// A Dirichlet datum and the free surface with gravity are constant over the
/// face, so what they impose is a map on the ghost state and becomes part of
/// the local flux solver; what is left of either condition is the offset the
/// local kernel adds. Here for the same two reasons as the assembly above: the
/// three solvers that fold do it word for word, and the kernels are generated
/// only where the pair's ghost-state index runs over the quantities of the
/// Riemann problem.
inline void foldBoundaryIntoFaceFlux(FaceType faceType,
                                     real* aPlusT,
                                     real* aMinusT,
                                     const real* matTinvData,
                                     const real* dirichletMap) {
  if (faceType == FaceType::Dirichlet) {
    kernel::foldDirichlet foldKrnl;
    foldKrnl.AplusT = aPlusT;
    foldKrnl.AminusT = aMinusT;
    foldKrnl.Tinv = matTinvData;
    foldKrnl.dirichletMap = dirichletMap;
    foldKrnl.execute();
  }

  if (faceType == FaceType::FreeSurfaceGravity) {
    kernel::foldFreeSurfaceGravity foldKrnl;
    // fsgMap is a constant; only the pool holds it
    foldKrnl.bindGlobals(seissol::Pool::host());
    foldKrnl.AplusT = aPlusT;
    foldKrnl.AminusT = aMinusT;
    foldKrnl.Tinv = matTinvData;
    foldKrnl.execute();
  }
}

} // namespace seissol::model::detail

#endif // SEISSOL_SRC_MODEL_GODUNOVFACEFLUX_H_
