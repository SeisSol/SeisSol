// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_PHYSICS_DRUCKERPRAGER_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_PHYSICS_DRUCKERPRAGER_H_

#include "Core/Common/Marker.h"

#include <algorithm>
#include <cmath>

/**
 * The plastic correction at a node, once for the host and the device: the Drucker-Prager yield
 * criterion, and the relaxation of the stress beyond it towards the yield surface (Wollherr et al.,
 * 2018). The host evaluates the nodes of a cell in a loop, the device one node per thread; the
 * conversions between the modal and the nodal stresses stay with them.
 */
namespace seissol::model::plasticity {

/// the stress components that the correction acts on: xx, yy, zz, xy, yz, xz
constexpr int StressComponents = 6;

/// what the yield criterion finds at a node
template <typename RealT>
struct NodeYield {
  /// the deviatoric stress s_ij = sigma_ij + sigma0_ij - m delta_ij
  RealT deviator[StressComponents]{};
  /// whether the stress lies beyond the yield surface, tau > tau_c
  bool yields{false};
  /// (tau_c / tau - 1) (1 - exp(-dt / Tv)) where the node yields, else 0
  RealT factor{};
};

/**
 * The yield criterion at a node whose total stress (sigma_ij + sigma0_ij) is `stress`:
 *   m = s_ii / 3,  s_ij = stress_ij - m delta_ij,  tau = sqrt(0.5 s_ij s_ji),
 *   tau_c = max(0, c cos(phi) - m sin(phi)),
 * where a shear component stands for the two entries ij and ji of the tensor.
 * `oneMinusIntegratingFactor` is 1 - exp(-dt / Tv), see Plasticity::computeRelaxTime.
 */
template <typename RealT>
SEISSOL_HOSTDEVICE inline NodeYield<RealT> yieldOfNode(const RealT (&stress)[StressComponents],
                                                       RealT cohesionTimesCosAngularFriction,
                                                       RealT sinAngularFriction,
                                                       RealT oneMinusIntegratingFactor) {
  NodeYield<RealT> node{};

  const RealT meanStress = (stress[0] + stress[1] + stress[2]) / static_cast<RealT>(3);
  for (int i = 0; i < StressComponents; ++i) {
    node.deviator[i] = i < 3 ? stress[i] - meanStress : stress[i];
  }

  const auto* s = node.deviator;
  RealT tau = static_cast<RealT>(0.5) * (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
  tau += (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
  tau = std::sqrt(tau);

  const RealT taulim = std::max(static_cast<RealT>(0.0),
                                cohesionTimesCosAngularFriction - meanStress * sinAngularFriction);

  node.yields = tau > taulim;
  node.factor = node.yields ? ((taulim / tau) - static_cast<RealT>(1.0)) * oneMinusIntegratingFactor
                            : static_cast<RealT>(0.0);
  return node;
}

/**
 * The relaxation of a node with the yield `node`. The relaxation changes the stress at a fixed
 * total strain by sigmaNew_ij - sigma_ij = yield s_ij (`stressUpdate`), where, by equation (13) of
 * Wollherr et al., sigmaNew_ij = f^* s_ij + m delta_ij - sigma0_ij with f^* = 1 + yield, so that it
 * turns the stress it releases into plastic strain, by equation (10):
 *   dstrain_ij = (sigma_ij - sigmaNew_ij) / (2 mu) = -yield s_ij / (2 mu)   (`pstrainUpdate`),
 * whatever the time step is relative to Tv, and for Tv = 0 as well; `mufactor` is 1 / (2 mu).
 * Returns the growth of the accumulated plastic strain eta = int sqrt(0.5 dstrain_ij/dt
 * dstrain_ij/dt) dt: the relaxation scales the deviatoric stress of the node without turning it, so
 * that its plastic strain keeps one direction over the step, and eta grows by sqrt(0.5 dstrain_ij
 * dstrain_ij) exactly. A node that does not yield (factor 0) changes nothing.
 */
template <typename RealT>
SEISSOL_HOSTDEVICE inline RealT relaxNode(const NodeYield<RealT>& node,
                                          RealT mufactor,
                                          RealT (&stressUpdate)[StressComponents],
                                          RealT (&pstrainUpdate)[StressComponents]) {
  RealT pstrainSq = 0;
  for (int i = 0; i < StressComponents; ++i) {
    stressUpdate[i] = node.deviator[i] * node.factor;
    pstrainUpdate[i] = -mufactor * stressUpdate[i];
    pstrainSq += (i < 3 ? static_cast<RealT>(0.5) : static_cast<RealT>(1.0)) * pstrainUpdate[i] *
                 pstrainUpdate[i];
  }
  return std::sqrt(pstrainSq);
}

} // namespace seissol::model::plasticity

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_PHYSICS_DRUCKERPRAGER_H_
