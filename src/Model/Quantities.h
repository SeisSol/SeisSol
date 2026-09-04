// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_QUANTITIES_H_
#define SEISSOL_SRC_MODEL_QUANTITIES_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace seissol::model {

/// How a group of quantities transforms under the face rotation.
enum class QuantityKind : std::uint8_t {
  /// One rotation invariant component, e.g. a pressure.
  Scalar,
  /// Three components transforming as a first-order tensor.
  Vector,
  /// Six components in Voigt order (xx, yy, zz, xy, yz, xz).
  SymTensor2,
};

constexpr std::size_t extentOf(QuantityKind kind) {
  switch (kind) {
  case QuantityKind::Scalar:
    return 1;
  case QuantityKind::Vector:
    return 3;
  case QuantityKind::SymTensor2:
    return 6;
  }
  return 0;
}

/// Whether and how a group takes part in the face-local quantities.
///
/// `Traction` and `Velocity` mark the groups carrying the mechanical traction
/// and velocity vectors. `Extra*` groups contribute one further row each, below
/// those vectors, and only through their normal component.
enum class FaceRole : std::uint8_t {
  None,
  Traction,
  Velocity,
  ExtraTraction,
  ExtraVelocity,
};

/// A contiguous run of quantities sharing a transformation behaviour.
struct QuantityGroup {
  std::string_view name;
  QuantityKind kind{QuantityKind::Scalar};
  FaceRole role{FaceRole::None};

  [[nodiscard]] constexpr std::size_t extent() const { return extentOf(kind); }
};

namespace detail {

template <std::size_t A, std::size_t B>
constexpr std::array<QuantityGroup, A + B> concat(const std::array<QuantityGroup, A>& first,
                                                  const std::array<QuantityGroup, B>& second) {
  std::array<QuantityGroup, A + B> result{};
  for (std::size_t i = 0; i < A; ++i) {
    result[i] = first[i];
  }
  for (std::size_t i = 0; i < B; ++i) {
    result[A + i] = second[i];
  }
  return result;
}

template <std::size_t A, std::size_t B, typename... Rest>
constexpr auto concat(const std::array<QuantityGroup, A>& first,
                      const std::array<QuantityGroup, B>& second,
                      const Rest&... rest) {
  return concat(concat(first, second), rest...);
}

} // namespace detail

/// Places `Repetitions` copies of a mechanism block behind the primary groups.
///
/// How often the block appears on the quantity axis is a property of the
/// solver, not of the material: the fused layout repeats it once per
/// relaxation mechanism, the split layout carries a single copy and keeps the
/// mechanism index as a separate tensor dimension.
template <std::size_t Repetitions, std::size_t N, std::size_t M>
constexpr std::array<QuantityGroup, N + M * Repetitions>
    withMechanisms(const std::array<QuantityGroup, N>& primary,
                   const std::array<QuantityGroup, M>& mechanism) {
  std::array<QuantityGroup, N + M * Repetitions> result{};
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = primary[i];
  }
  for (std::size_t rep = 0; rep < Repetitions; ++rep) {
    for (std::size_t i = 0; i < M; ++i) {
      result[N + rep * M + i] = mechanism[i];
    }
  }
  return result;
}

/// Number of quantities the groups occupy.
template <std::size_t N>
constexpr std::size_t totalExtent(const std::array<QuantityGroup, N>& groups) {
  std::size_t sum = 0;
  for (const auto& group : groups) {
    sum += group.extent();
  }
  return sum;
}

/// Quantity index at which the first group with `role` starts.
template <std::size_t N>
constexpr std::size_t roleOffset(const std::array<QuantityGroup, N>& groups, FaceRole role) {
  std::size_t offset = 0;
  for (const auto& group : groups) {
    if (group.role == role) {
      return offset;
    }
    offset += group.extent();
  }
  return totalExtent(groups);
}

/// Weight of a component of `kind` in the trace of the tensor it represents.
///
/// A scalar stands for an isotropic tensor, so its one component enters the
/// trace three times over; a symmetric second-order tensor contributes its
/// three diagonal components once each. Energy expressions written for one of
/// them therefore carry over to the other by these weights alone.
constexpr double traceWeight(QuantityKind kind) {
  return kind == QuantityKind::Scalar ? 3.0 : 1.0;
}

/// Components of `kind` that enter the trace.
constexpr std::size_t traceComponents(QuantityKind kind) {
  return kind == QuantityKind::Scalar ? 1 : 3;
}

/// Kind of the first group with `role`.
template <std::size_t N>
constexpr QuantityKind roleKind(const std::array<QuantityGroup, N>& groups, FaceRole role) {
  for (const auto& group : groups) {
    if (group.role == role) {
      return group.kind;
    }
  }
  return QuantityKind::Scalar;
}

/// Number of components in the first group with `role`, or zero if absent.
template <std::size_t N>
constexpr std::size_t roleExtent(const std::array<QuantityGroup, N>& groups, FaceRole role) {
  for (const auto& group : groups) {
    if (group.role == role) {
      return group.extent();
    }
  }
  return 0;
}

/// True if the groups can describe a layout of `numQuantities` quantities.
template <std::size_t N>
constexpr bool quantitiesWellFormed(const std::array<QuantityGroup, N>& groups,
                                    std::size_t numQuantities) {
  std::size_t tractions = 0;
  std::size_t velocities = 0;
  for (const auto& group : groups) {
    if (group.name.empty()) {
      return false;
    }
    tractions += static_cast<std::size_t>(group.role == FaceRole::Traction);
    velocities += static_cast<std::size_t>(group.role == FaceRole::Velocity);
  }
  return tractions <= 1 && velocities <= 1 && totalExtent(groups) == numQuantities;
}

// Reusable fragments. The materials compose their layout from these rather
// than spelling out the same groups again.

inline constexpr std::array AcousticQuantities{
    QuantityGroup{"pprime", QuantityKind::Scalar, FaceRole::Traction},
    QuantityGroup{"v", QuantityKind::Vector, FaceRole::Velocity},
};

inline constexpr std::array ElasticQuantities{
    QuantityGroup{"s", QuantityKind::SymTensor2, FaceRole::Traction},
    QuantityGroup{"v", QuantityKind::Vector, FaceRole::Velocity},
};

inline constexpr std::array PoroelasticExtraQuantities{
    QuantityGroup{"p", QuantityKind::Scalar, FaceRole::ExtraTraction},
    QuantityGroup{"vf", QuantityKind::Vector, FaceRole::ExtraVelocity},
};

/// One relaxation mechanism of a material whose traction is a stress tensor.
inline constexpr std::array ElasticMechanismQuantities{
    QuantityGroup{"theta", QuantityKind::SymTensor2},
};

/// One relaxation mechanism of a material whose traction is a scalar.
inline constexpr std::array AcousticMechanismQuantities{
    QuantityGroup{"theta", QuantityKind::Scalar},
};

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_QUANTITIES_H_
