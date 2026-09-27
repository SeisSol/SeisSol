// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_DUAL_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_DUAL_H_

#include "Common/Marker.h"

#include <cmath>

namespace seissol::dr::friction_law {

/**
 * A forward-mode dual number: a value and the derivative of that value with respect to one
 * chosen variable, carried along every operation.
 *
 * A function written once against a generic scalar type yields its value when instantiated with
 * `real` and its value together with the exact derivative when instantiated with `Dual<real>`.
 * Seeding is `Dual<real>(x, 1)` for the variable to differentiate by; every other quantity enters
 * as a constant and carries a zero derivative.
 *
 * The friction laws use this to hand the slip-rate inversion a derivative that follows the whole
 * chain -- state variable, friction coefficient, effective normal stress, pore pressure -- without
 * a second set of functions that has to stay in step with the first. Differentiating by hand would
 * need one derivative function per layer; a finite difference is not an option, since the state
 * variable is a small perturbation on a large constant and differencing it loses every digit in
 * single precision.
 *
 * The type is trivially copyable and every operation is inlined, so it goes through device code
 * and vectorised loops as a pair of scalars.
 */
template <typename T>
struct Dual {
  T value{};
  T derivative{};

  Dual() = default;

  /// a constant: it carries no dependence on the differentiation variable
  SEISSOL_HOSTDEVICE constexpr Dual(T constantValue) // NOLINT(google-explicit-constructor)
      : value(constantValue), derivative(0) {}

  SEISSOL_HOSTDEVICE constexpr Dual(T nodeValue, T nodeDerivative)
      : value(nodeValue), derivative(nodeDerivative) {}
};

/**
 * The value of a scalar or of a dual number. Used wherever a formula has to branch on the value:
 * the branch picks a formula, and the chosen formula is what gets differentiated.
 */
SEISSOL_HOSTDEVICE constexpr float valueOf(float x) { return x; }
SEISSOL_HOSTDEVICE constexpr double valueOf(double x) { return x; }
template <typename T>
SEISSOL_HOSTDEVICE constexpr T valueOf(Dual<T> x) {
  return x.value;
}

/// The derivative carried by a scalar is zero; the one carried by a dual number is its own.
SEISSOL_HOSTDEVICE constexpr float derivativeOf(float /*x*/) { return 0; }
SEISSOL_HOSTDEVICE constexpr double derivativeOf(double /*x*/) { return 0; }
template <typename T>
SEISSOL_HOSTDEVICE constexpr T derivativeOf(Dual<T> x) {
  return x.derivative;
}

// arithmetic

template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator+(Dual<T> a, Dual<T> b) {
  return {a.value + b.value, a.derivative + b.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator-(Dual<T> a, Dual<T> b) {
  return {a.value - b.value, a.derivative - b.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator-(Dual<T> a) {
  return {-a.value, -a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator*(Dual<T> a, Dual<T> b) {
  return {a.value * b.value, a.derivative * b.value + a.value * b.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator/(Dual<T> a, Dual<T> b) {
  const T inverse = T(1) / b.value;
  return {a.value * inverse, (a.derivative - a.value * inverse * b.derivative) * inverse};
}

// elementary functions; each evaluates its transcendental once and reuses it for the derivative

template <typename T>
SEISSOL_HOSTDEVICE Dual<T> exp(Dual<T> a) {
  const T e = std::exp(a.value);
  return {e, e * a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> expm1(Dual<T> a) {
  const T e = std::expm1(a.value);
  return {e, (e + T(1)) * a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> log(Dual<T> a) {
  return {std::log(a.value), a.derivative / a.value};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> log1p(Dual<T> a) {
  return {std::log1p(a.value), a.derivative / (T(1) + a.value)};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> sqrt(Dual<T> a) {
  const T s = std::sqrt(a.value);
  return {s, a.derivative / (T(2) * s)};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> asinh(Dual<T> a) {
  return {std::asinh(a.value), a.derivative / std::sqrt(T(1) + a.value * a.value)};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> abs(Dual<T> a) {
  return a.value >= T(0) ? a : -a;
}
/// only the exponent is a plain number here, which is what the friction laws need
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> pow(Dual<T> a, T exponent) {
  const T p = std::pow(a.value, exponent);
  return {p, exponent * p / a.value * a.derivative};
}

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_DUAL_H_
