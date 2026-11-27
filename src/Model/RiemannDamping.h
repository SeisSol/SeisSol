// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_RIEMANNDAMPING_H_
#define SEISSOL_SRC_MODEL_RIEMANNDAMPING_H_

#include "Model/Quantities.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace seissol::model {

/**
 * The part of the Riemann problem a component of the face-aligned quantity vector belongs to (the
 * normal of the face is the x axis): the traction on the face or the velocity, each with its normal
 * component and its tangential ones. A scalar traction (a pressure) and a scalar velocity act along
 * the normal. The stress components that do not act on the face (yy, zz, yz) and the quantities
 * outside the Riemann problem belong to none.
 */
enum class FaceComponent : std::uint8_t {
  None,
  TractionNormal,
  TractionTangential,
  VelocityNormal,
  VelocityTangential,
};

template <std::size_t N>
constexpr FaceComponent faceComponentOf(const std::array<QuantityGroup, N>& groups,
                                        std::size_t index) {
  std::size_t offset = 0;
  for (const auto& group : groups) {
    const auto extent = group.extent();
    if (index < offset + extent) {
      const auto component = index - offset;
      const bool traction =
          group.role == FaceRole::Traction || group.role == FaceRole::ExtraTraction;
      const bool velocity =
          group.role == FaceRole::Velocity || group.role == FaceRole::ExtraVelocity;
      if (!traction && !velocity) {
        return FaceComponent::None;
      }
      bool normal = false;
      bool tangential = false;
      switch (group.kind) {
      case QuantityKind::Scalar:
        normal = true;
        break;
      case QuantityKind::Vector:
        normal = component == 0;
        tangential = component != 0;
        break;
      case QuantityKind::SymTensor2:
        // Voigt order xx yy zz xy yz xz: the traction on the face is (xx, xy, xz)
        normal = component == 0;
        tangential = component == 3 || component == 5;
        break;
      }
      if (normal) {
        return traction ? FaceComponent::TractionNormal : FaceComponent::VelocityNormal;
      }
      if (tangential) {
        return traction ? FaceComponent::TractionTangential : FaceComponent::VelocityTangential;
      }
      return FaceComponent::None;
    }
    offset += extent;
  }
  return FaceComponent::None;
}

/**
 * Whether a component of the face-aligned quantity vector belongs to the extra pair of the Riemann
 * problem (FaceRole::ExtraTraction and FaceRole::ExtraVelocity: the pressure and the velocity of
 * the fluid of a poroelastic solid) rather than to the traction and velocity of the material
 * itself.
 */
template <std::size_t N>
constexpr bool isExtraComponent(const std::array<QuantityGroup, N>& groups, std::size_t index) {
  std::size_t offset = 0;
  for (const auto& group : groups) {
    if (index < offset + group.extent()) {
      return group.role == FaceRole::ExtraTraction || group.role == FaceRole::ExtraVelocity;
    }
    offset += group.extent();
  }
  return false;
}

constexpr bool isTraction(FaceComponent component) {
  return component == FaceComponent::TractionNormal ||
         component == FaceComponent::TractionTangential;
}

constexpr bool isVelocity(FaceComponent component) {
  return component == FaceComponent::VelocityNormal ||
         component == FaceComponent::VelocityTangential;
}

constexpr bool isNormal(FaceComponent component) {
  return component == FaceComponent::TractionNormal || component == FaceComponent::VelocityNormal;
}

/**
 * The shares of the upwind damping in a numerical flux (fluxdamping), one per pair of the Riemann
 * problem, kind of jump and direction: of the velocity jumps of the material (normal, tangential),
 * which act through the traction of the Riemann state, of its traction jumps (normal, tangential),
 * which act through its velocity, and of the jumps of the extra pair (the fluid of a poroelastic
 * solid), which has a normal direction only: of its velocity and of its traction (the pressure).
 * A component without a partner of its own direction in its pair (the tangential velocity of the
 * fluid) takes the share of that direction in the pair of the material. All 1 is the Godunov flux,
 * all 0 its part without damping (on a face within one material the centered flux).
 */
using DampingShares = std::array<double, 6>;

/// The shares of the Godunov flux itself.
inline constexpr DampingShares UpwindDamping{1, 1, 1, 1, 1, 1};

constexpr bool isUpwindDamping(const DampingShares& shares) {
  for (const auto share : shares) {
    if (share != 1.0) {
      return false;
    }
  }
  return true;
}

/// The index in DampingShares of the share of a component of the jumps of the velocity
/// (tractionJumps = false) or of the traction.
constexpr std::size_t dampingSlot(FaceComponent component, bool extra, bool tractionJumps) {
  if (extra && isNormal(component)) {
    return tractionJumps ? 5 : 4;
  }
  return (tractionJumps ? 2 : 0) + (isNormal(component) ? 0 : 1);
}

/**
 * Scales the upwind damping of the Godunov state of a face (transposed: row = the quantity the
 * state is taken from, column = the component of the state). The damping of the velocity jumps is
 * the part of the traction of the state that depends on the velocities, the damping of the traction
 * jumps the part of its velocity that depends on the tractions; both are symmetric positive
 * semi-definite operators on the face components (also with a material contrast). Each is scaled
 * symmetrically, entry (a, b) by sqrt(s_a s_b) with the shares (dampingSlot) of the component a of
 * the state and b of the quantity, so that it stays positive semi-definite where components of
 * other directions or pairs couple (anisotropic or poroelastic solids). What the state takes from
 * the same kind of quantity, the (impedance-weighted) averages, is kept. The scaled entries
 * cancel in qGodLocal + qGodNeighbor, which therefore stays the identity.
 */
template <typename MaterialT, typename LocalViewT, typename NeighborViewT>
void scaleRiemannDamping(LocalViewT& qGodLocal,
                         NeighborViewT& qGodNeighbor,
                         const DampingShares& shares) {
  constexpr auto Groups = MaterialT::PrimaryGroups;
  constexpr std::size_t Count = totalExtent(Groups);
  const auto share = [&](std::size_t index, bool tractionJumps) {
    return shares[dampingSlot(
        faceComponentOf(Groups, index), isExtraComponent(Groups, index), tractionJumps)];
  };
  for (std::size_t from = 0; from < Count; ++from) {
    const auto quantity = faceComponentOf(Groups, from);
    for (std::size_t to = 0; to < Count; ++to) {
      const auto state = faceComponentOf(Groups, to);
      double factor = 1.0;
      if (isTraction(state) && isVelocity(quantity)) {
        factor = std::sqrt(share(to, false) * share(from, false));
      } else if (isVelocity(state) && isTraction(quantity)) {
        factor = std::sqrt(share(to, true) * share(from, true));
      } else {
        continue;
      }
      if (qGodLocal.isInRange(from, to)) {
        using RealT = std::decay_t<decltype(qGodLocal(from, to))>;
        qGodLocal(from, to) = static_cast<RealT>(qGodLocal(from, to) * factor);
      }
      if (qGodNeighbor.isInRange(from, to)) {
        using RealT = std::decay_t<decltype(qGodNeighbor(from, to))>;
        qGodNeighbor(from, to) = static_cast<RealT>(qGodNeighbor(from, to) * factor);
      }
    }
  }
}

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_RIEMANNDAMPING_H_
