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
#include <type_traits>

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
 *
 * First order only: nesting a dual inside a dual would give second derivatives, but the elementary
 * functions below reach for the standard library by qualified name, which a nested value would not
 * find. Nothing here needs a second derivative.
 *
 * Comparisons are deliberately absent. A formula that has to branch does so on `valueOf`, which
 * says at the call site that the branch picks a formula and that the picked formula is what
 * carries the derivative -- a comparison spelled `a < b` would hide that decision.
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

  SEISSOL_HOSTDEVICE constexpr Dual& operator+=(Dual b) {
    value += b.value;
    derivative += b.derivative;
    return *this;
  }
  SEISSOL_HOSTDEVICE constexpr Dual& operator-=(Dual b) {
    value -= b.value;
    derivative -= b.derivative;
    return *this;
  }
  SEISSOL_HOSTDEVICE constexpr Dual& operator*=(Dual b) {
    derivative = derivative * b.value + value * b.derivative;
    value *= b.value;
    return *this;
  }
  SEISSOL_HOSTDEVICE constexpr Dual& operator/=(Dual b) {
    const T inverse = T(1) / b.value;
    derivative = (derivative - value * inverse * b.derivative) * inverse;
    value *= inverse;
    return *this;
  }
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

// ---------------------------------------------------------------------------
// arithmetic
// ---------------------------------------------------------------------------

template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator+(Dual<T> a, Dual<T> b) {
  return {a.value + b.value, a.derivative + b.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator-(Dual<T> a, Dual<T> b) {
  return {a.value - b.value, a.derivative - b.derivative};
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
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator+(Dual<T> a) {
  return a;
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator-(Dual<T> a) {
  return {-a.value, -a.derivative};
}

/// A constant on either side, without the detour through a dual that carries a zero derivative.
/// The scalar has to be the dual's own, so that a mixed-precision expression stays visible as the
/// cast it is.
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator+(Dual<T> a, T b) {
  return {a.value + b, a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator+(T a, Dual<T> b) {
  return {a + b.value, b.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator-(Dual<T> a, T b) {
  return {a.value - b, a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator-(T a, Dual<T> b) {
  return {a - b.value, -b.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator*(Dual<T> a, T b) {
  return {a.value * b, a.derivative * b};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator*(T a, Dual<T> b) {
  return {a * b.value, a * b.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator/(Dual<T> a, T b) {
  const T inverse = T(1) / b;
  return {a.value * inverse, a.derivative * inverse};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> operator/(T a, Dual<T> b) {
  const T inverse = T(1) / b.value;
  const T quotient = a * inverse;
  return {quotient, -quotient * inverse * b.derivative};
}

// ---------------------------------------------------------------------------
// classification
// ---------------------------------------------------------------------------

/// A dual number is finite when both halves are: a finite value beside an infinite derivative is
/// what a Newton step that has left the domain looks like, and the assertions want to catch it.
template <typename T>
SEISSOL_HOSTDEVICE bool isfinite(Dual<T> a) {
  return std::isfinite(a.value) && std::isfinite(a.derivative);
}
template <typename T>
SEISSOL_HOSTDEVICE bool isnan(Dual<T> a) {
  return std::isnan(a.value) || std::isnan(a.derivative);
}
template <typename T>
SEISSOL_HOSTDEVICE bool isinf(Dual<T> a) {
  return std::isinf(a.value) || std::isinf(a.derivative);
}

// ---------------------------------------------------------------------------
// elementary functions
//
// Each evaluates its transcendental once and builds the derivative from the result, and each
// derivative is arranged so that it keeps the range of the value beside it: a formula that
// overflows an intermediate would hand the inversion a zero or infinite slope where the function
// itself is perfectly well behaved.
// ---------------------------------------------------------------------------

template <typename T>
SEISSOL_HOSTDEVICE Dual<T> exp(Dual<T> a) {
  const T e = std::exp(a.value);
  return {e, e * a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> exp2(Dual<T> a) {
  constexpr T Ln2 = 0.69314718055994530942;
  const T e = std::exp2(a.value);
  return {e, e * Ln2 * a.derivative};
}
/// exp(v) = expm1(v) + 1 reuses the call, at the price of the low digits of the derivative once
/// expm1 approaches -1 -- from v = -37 down there are none left. That is where the derivative is
/// exponentially small and, in the state variable it integrates, multiplies a bounded steady state;
/// a second exp for the sake of it would cost one transcendental per Newton step.
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
SEISSOL_HOSTDEVICE Dual<T> log2(Dual<T> a) {
  constexpr T Ln2 = 0.69314718055994530942;
  return {std::log2(a.value), a.derivative / (a.value * Ln2)};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> log10(Dual<T> a) {
  constexpr T Ln10 = 2.30258509299404568402;
  return {std::log10(a.value), a.derivative / (a.value * Ln10)};
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
SEISSOL_HOSTDEVICE Dual<T> cbrt(Dual<T> a) {
  const T c = std::cbrt(a.value);
  return {c, a.derivative / (T(3) * c * c)};
}
/// hypot divides through the result first, so that a large leg cannot overflow the numerator
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> hypot(Dual<T> a, Dual<T> b) {
  const T h = std::hypot(a.value, b.value);
  const T inverse = T(1) / h;
  return {h, (a.value * inverse) * a.derivative + (b.value * inverse) * b.derivative};
}

/// a constant exponent, which is the common case in the friction laws
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> pow(Dual<T> a, T exponent) {
  const T p = std::pow(a.value, exponent);
  return {p, exponent * p / a.value * a.derivative};
}
/// both the base and the exponent may depend on the variable: d(a^b) = a^b (b' ln a + b a'/a)
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> pow(Dual<T> a, Dual<T> b) {
  const T p = std::pow(a.value, b.value);
  return {p, p * (b.derivative * std::log(a.value) + b.value * a.derivative / a.value)};
}
/// a constant base: d(a^b) = a^b b' ln a
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> pow(T base, Dual<T> b) {
  const T p = std::pow(base, b.value);
  return {p, p * std::log(base) * b.derivative};
}

// ---------------------------------------------------------------------------
// trigonometric
// ---------------------------------------------------------------------------

template <typename T>
SEISSOL_HOSTDEVICE Dual<T> sin(Dual<T> a) {
  return {std::sin(a.value), std::cos(a.value) * a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> cos(Dual<T> a) {
  return {std::cos(a.value), -std::sin(a.value) * a.derivative};
}
/// the derivative is built from the tangent rather than from a second call to cos
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> tan(Dual<T> a) {
  const T t = std::tan(a.value);
  return {t, (T(1) + t * t) * a.derivative};
}
/// (1-v)(1+v) rather than 1-v*v: the factored form keeps its digits as |v| approaches one, which
/// is where the derivative matters most
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> asin(Dual<T> a) {
  return {std::asin(a.value), a.derivative / std::sqrt((T(1) - a.value) * (T(1) + a.value))};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> acos(Dual<T> a) {
  return {std::acos(a.value), -a.derivative / std::sqrt((T(1) - a.value) * (T(1) + a.value))};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> atan(Dual<T> a) {
  return {std::atan(a.value), a.derivative / (T(1) + a.value * a.value)};
}
/// d atan2(y, x) = (x y' - y x') / (x^2 + y^2), scaled by the hypotenuse so that neither the
/// numerator nor the denominator leaves the range the arguments themselves occupy
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> atan2(Dual<T> y, Dual<T> x) {
  const T h = std::hypot(y.value, x.value);
  const T inverse = T(1) / h;
  return {std::atan2(y.value, x.value),
          ((x.value * inverse) * y.derivative - (y.value * inverse) * x.derivative) * inverse};
}

// ---------------------------------------------------------------------------
// hyperbolic
// ---------------------------------------------------------------------------

template <typename T>
SEISSOL_HOSTDEVICE Dual<T> sinh(Dual<T> a) {
  return {std::sinh(a.value), std::cosh(a.value) * a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> cosh(Dual<T> a) {
  return {std::cosh(a.value), std::sinh(a.value) * a.derivative};
}
/// the one derivative here built from a second call rather than from the value beside it:
/// 1 - tanh(v)^2 cancels as tanh approaches one, and by |v| = 19 there is nothing left of it, while
/// 1 / cosh(v)^2 stays exact until it underflows
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> tanh(Dual<T> a) {
  const T c = std::cosh(a.value);
  return {std::tanh(a.value), a.derivative / (c * c)};
}
/// hypot rather than sqrt(1 + v*v): asinh is evaluated at arguments reaching to the edge of the
/// range, where 1 + v*v has long overflowed and would hand the inversion a zero derivative -- and
/// with it an infinite step. hypot costs about a fifth of the asinh it accompanies.
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> asinh(Dual<T> a) {
  return {std::asinh(a.value), a.derivative / std::hypot(T(1), a.value)};
}
/// two square roots rather than sqrt(v*v - 1): the product form neither overflows at a large
/// argument nor cancels as the argument approaches one
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> acosh(Dual<T> a) {
  return {std::acosh(a.value),
          a.derivative / (std::sqrt(a.value - T(1)) * std::sqrt(a.value + T(1)))};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> atanh(Dual<T> a) {
  return {std::atanh(a.value), a.derivative / ((T(1) - a.value) * (T(1) + a.value))};
}

// ---------------------------------------------------------------------------
// error function
// ---------------------------------------------------------------------------

template <typename T>
SEISSOL_HOSTDEVICE Dual<T> erf(Dual<T> a) {
  constexpr T TwoOverSqrtPi = 1.12837916709551257390;
  return {std::erf(a.value), TwoOverSqrtPi * std::exp(-a.value * a.value) * a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> erfc(Dual<T> a) {
  constexpr T TwoOverSqrtPi = 1.12837916709551257390;
  return {std::erfc(a.value), -TwoOverSqrtPi * std::exp(-a.value * a.value) * a.derivative};
}

// ---------------------------------------------------------------------------
// sign and selection
//
// Each of these picks one of its arguments by value; the picked one keeps its own derivative,
// which is the one-sided derivative of the piecewise function. At a kink -- abs at zero, fmax of
// two equal arguments -- that is the derivative from the right.
// ---------------------------------------------------------------------------

template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> abs(Dual<T> a) {
  return a.value >= T(0) ? a : -a;
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> fabs(Dual<T> a) {
  return abs(a);
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> copysign(Dual<T> a, Dual<T> b) {
  const T sign = std::copysign(T(1), a.value) * std::copysign(T(1), b.value);
  return {std::copysign(a.value, b.value), sign * a.derivative};
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> fmax(Dual<T> a, Dual<T> b) {
  return a.value >= b.value ? a : b;
}
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> fmin(Dual<T> a, Dual<T> b) {
  return a.value <= b.value ? a : b;
}
/// the positive difference: zero, with a zero derivative, once the arguments cross
template <typename T>
SEISSOL_HOSTDEVICE constexpr Dual<T> fdim(Dual<T> a, Dual<T> b) {
  return a.value > b.value ? a - b : Dual<T>(0);
}

// ---------------------------------------------------------------------------
// piecewise constant
//
// These are constant between their steps, so the derivative is zero -- which is the derivative
// everywhere except on the steps themselves, where none exists.
// ---------------------------------------------------------------------------

template <typename T>
SEISSOL_HOSTDEVICE Dual<T> floor(Dual<T> a) {
  return {std::floor(a.value), T(0)};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> ceil(Dual<T> a) {
  return {std::ceil(a.value), T(0)};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> round(Dual<T> a) {
  return {std::round(a.value), T(0)};
}
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> trunc(Dual<T> a) {
  return {std::trunc(a.value), T(0)};
}

// ---------------------------------------------------------------------------
// composites
// ---------------------------------------------------------------------------

/// a b + c in one rounding, and the product rule on top of it
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> fma(Dual<T> a, Dual<T> b, Dual<T> c) {
  return {std::fma(a.value, b.value, c.value),
          a.derivative * b.value + a.value * b.derivative + c.derivative};
}
/// a scaling by a power of two, which is exact and therefore scales the derivative with it
template <typename T>
SEISSOL_HOSTDEVICE Dual<T> ldexp(Dual<T> a, int exponent) {
  return {std::ldexp(a.value, exponent), std::ldexp(a.derivative, exponent)};
}

// Deliberately absent: the gamma functions, whose derivatives need the digamma function that the
// standard library does not provide, and fmod and remainder, which are piecewise in their second
// argument. Nothing in the friction laws asks for them.

// ---------------------------------------------------------------------------
// precision
// ---------------------------------------------------------------------------

/// Carry a value from one precision to another. A dual number takes its derivative along, a plain
/// scalar is simply converted, so a formula can move between the precisions its parts are stated
/// in without knowing whether it is being differentiated.
template <typename T, typename U>
SEISSOL_HOSTDEVICE constexpr Dual<T> dualCast(Dual<U> x) {
  return {static_cast<T>(x.value), static_cast<T>(x.derivative)};
}
template <typename T, typename U, std::enable_if_t<std::is_floating_point_v<U>, int> = 0>
SEISSOL_HOSTDEVICE constexpr T dualCast(U x) {
  return static_cast<T>(x);
}

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_DUAL_H_
