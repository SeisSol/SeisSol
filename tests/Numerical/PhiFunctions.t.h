// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Numerical/PhiFunctions.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace seissol::unit_test {

namespace phifunctionstest {

/**
 * phi_n(x) at 120 decimal digits, rounded once to the nearest double. Produced with
 * mpmath as hyp1f1(1, n+1, x) / n!, which is an evaluation path independent of the one
 * under test. The sample points bracket the switch at |x| = n+1 from both sides and cover
 * the range in which the defining difference exp(x) - sum_{k<n} x^k/k! loses every digit.
 */
constexpr std::size_t MaxOrder = 8;
constexpr std::size_t SampleCount = 17;

constexpr double ReferenceArguments[MaxOrder + 1][SampleCount] = {
    // order 0
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.ffffffffffc7bp-1,
        -0x1.ffffffffffc7bp-1,
        0x1.00000000001c2p+0,
        -0x1.00000000001c2p+0,
        0x1.4000000000000p+1,
        -0x1.4000000000000p+1,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 1
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.ffffffffffc7bp+0,
        -0x1.ffffffffffc7bp+0,
        0x1.00000000001c2p+1,
        -0x1.00000000001c2p+1,
        0x1.4000000000000p+2,
        -0x1.4000000000000p+2,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 2
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.7fffffffffd5cp+1,
        -0x1.7fffffffffd5cp+1,
        0x1.80000000002a3p+1,
        -0x1.80000000002a3p+1,
        0x1.e000000000000p+2,
        -0x1.e000000000000p+2,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 3
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.ffffffffffc7bp+1,
        -0x1.ffffffffffc7bp+1,
        0x1.00000000001c2p+2,
        -0x1.00000000001c2p+2,
        0x1.4000000000000p+3,
        -0x1.4000000000000p+3,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 4
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.3fffffffffdcdp+2,
        -0x1.3fffffffffdcdp+2,
        0x1.4000000000232p+2,
        -0x1.4000000000232p+2,
        0x1.9000000000000p+3,
        -0x1.9000000000000p+3,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 5
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.7fffffffffd5cp+2,
        -0x1.7fffffffffd5cp+2,
        0x1.80000000002a3p+2,
        -0x1.80000000002a3p+2,
        0x1.e000000000000p+3,
        -0x1.e000000000000p+3,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 6
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.bfffffffffcecp+2,
        -0x1.bfffffffffcecp+2,
        0x1.c000000000314p+2,
        -0x1.c000000000314p+2,
        0x1.1800000000000p+4,
        -0x1.1800000000000p+4,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 7
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.ffffffffffc7bp+2,
        -0x1.ffffffffffc7bp+2,
        0x1.00000000001c2p+3,
        -0x1.00000000001c2p+3,
        0x1.4000000000000p+4,
        -0x1.4000000000000p+4,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
    // order 8
    {
        0x0.0p+0,
        0x1.19799812dea11p-40,
        -0x1.19799812dea11p-40,
        0x1.0c6f7a0b5ed8dp-20,
        -0x1.0c6f7a0b5ed8dp-20,
        0x1.0000000000000p-2,
        -0x1.0000000000000p-2,
        0x1.0000000000000p+0,
        -0x1.0000000000000p+0,
        0x1.1fffffffffe05p+3,
        -0x1.1fffffffffe05p+3,
        0x1.20000000001fap+3,
        -0x1.20000000001fap+3,
        0x1.6800000000000p+4,
        -0x1.6800000000000p+4,
        0x1.2800000000000p+5,
        -0x1.2800000000000p+5,
    },
};

constexpr double ReferenceValues[MaxOrder + 1][SampleCount] = {
    // order 0
    {
        0x1.0000000000000p+0,
        0x1.0000000001198p+0,
        0x1.fffffffffdcd1p-1,
        0x1.000010c6f82d7p+0,
        0x1.ffffde7211d81p-1,
        0x1.48b5e3c3e8186p+0,
        0x1.8ebef9eac820bp-1,
        0x1.5bf0a8b145769p+1,
        0x1.78b56362cef38p-2,
        0x1.5bf0a8b145505p+1,
        0x1.78b56362cf1cfp-2,
        0x1.5bf0a8b1459cdp+1,
        0x1.78b56362ceca2p-2,
        0x1.85d6fd931e0bbp+3,
        0x1.50385c094f425p-4,
        0x1.4d13fbb1a001ap+53,
        0x1.898471fca6055p-54,
    },
    // order 1
    {
        0x1.0000000000000p+0,
        0x1.00000000008ccp+0,
        0x1.fffffffffee68p-1,
        0x1.000008637bff4p+0,
        0x1.ffffef3908bd2p-1,
        0x1.22d78f0fa061ap+0,
        0x1.c5041854df7d4p-1,
        0x1.b7e151628aed3p+0,
        0x1.43a54e4e98864p-1,
        0x1.98e64b8d4d9fdp+1,
        0x1.bab55571021a4p-2,
        0x1.98e64b8d4e15ep+1,
        0x1.bab5557101d76p-2,
        0x1.d7b8dc24d1f4bp+4,
        0x1.96d7133665114p-3,
        0x1.201148624531dp+48,
        0x1.bacf914c1bacfp-6,
    },
    // order 2
    {
        0x1.0000000000000p-1,
        0x1.00000000005ddp-1,
        0x1.ffffffffff446p-2,
        0x1.00000597a7f7bp-1,
        0x1.fffff4d0b06e7p-2,
        0x1.16bc787d030cdp-1,
        0x1.d7df3d590415dp-2,
        0x1.6fc2a2c515da5p-1,
        0x1.78b56362cef38p-2,
        0x1.c98b4e28da794p+0,
        0x1.d270c05b06a6bp-3,
        0x1.c98b4e28db163p+0,
        0x1.d270c05b06683p-3,
        0x1.ffdeadfcb9eecp+4,
        0x1.d95b17ab99a24p-4,
        0x1.f2476872a131dp+42,
        0x1.aed7cba3ff408p-6,
    },
    // order 3
    {
        0x1.5555555555555p-3,
        0x1.5555555555b33p-3,
        0x1.5555555554f78p-3,
        0x1.55555aecfd485p-3,
        0x1.55554fbdad87ep-3,
        0x1.6bc787d030ccfp-3,
        0x1.41061537df518p-3,
        0x1.bf0a8b1457695p-3,
        0x1.0e95393a62190p-3,
        0x1.4cc902e27364cp-1,
        0x1.3ed3eaa47e006p-4,
        0x1.4cc902e273e63p-1,
        0x1.3ed3eaa47dd80p-4,
        0x1.5f728c42e19fcp+4,
        0x1.4fdf2304a2172p-5,
        0x1.aef1a670fa019p+37,
        0x1.a385a1f722a29p-7,
    },
    // order 4
    {
        0x1.5555555555555p-5,
        0x1.5555555555a06p-5,
        0x1.55555555550a4p-5,
        0x1.555559cea87bap-5,
        0x1.555550dc02481p-5,
        0x1.672327adb7798p-5,
        0x1.44f401d7603dcp-5,
        0x1.a6d4d6fc08500p-5,
        0x1.1b00706bccf14p-5,
        0x1.656eca3773b54p-3,
        0x1.437cd1081fd6bp-6,
        0x1.656eca37744edp-3,
        0x1.437cd1081faf0p-6,
        0x1.5f2b14d806d39p+3,
        0x1.5b7beea88e3e6p-7,
        0x1.74b551add0264p+32,
        0x1.1087c481a1e02p-8,
    },
    // order 5
    {
        0x1.1111111111111p-7,
        0x1.1111111111432p-7,
        0x1.1111111110df0p-7,
        0x1.1111140c9dd41p-7,
        0x1.11110e15845c6p-7,
        0x1.1cdd258622425p-7,
        0x1.061537df51798p-7,
        0x1.45fe069acbea9p-7,
        0x1.d2a7274c4320ep-8,
        0x1.2fdbf0a13bcf1p-5,
        0x1.053f2119d387ep-8,
        0x1.2fdbf0a13c5cfp-5,
        0x1.053f2119d3687p-8,
        0x1.1346b71171511p+2,
        0x1.1c0007efb63bcp-9,
        0x1.4257a09649ee5p+27,
        0x1.09be956c384d4p-10,
    },
    // order 6
    {
        0x1.6c16c16c16c17p-10,
        0x1.6c16c16c16faap-10,
        0x1.6c16c16c16884p-10,
        0x1.6c16c4d4b79f4p-10,
        0x1.6c16be0375f1ep-10,
        0x1.79828ea226278p-10,
        0x1.5f7b2637f2f28p-10,
        0x1.a767ac4dd6cbdp-10,
        0x1.3debeb577c04fp-10,
        0x1.ab2e37579ed39p-8,
        0x1.5ea57e8123439p-11,
        0x1.ab2e37579fa99p-8,
        0x1.5ea57e812319fp-11,
        0x1.62c8365bae85ap+0,
        0x1.80479180b31a0p-12,
        0x1.16c85388aaf3fp+22,
        0x1.9edf8e82a9904p-13,
    },
    // order 7
    {
        0x1.a01a01a01a01ap-13,
        0x1.a01a01a01a3adp-13,
        0x1.a01a01a019c87p-13,
        0x1.a01a0508badebp-13,
        0x1.a019fe3779315p-13,
        0x1.ad79a6c1ecc1bp-13,
        0x1.9373668479de6p-13,
        0x1.da87570e00531p-13,
        0x1.7156b0a4d5e3cp-13,
        0x1.ffb35b099710cp-11,
        0x1.92b5df68741afp-14,
        0x1.ffb35b0998215p-11,
        0x1.92b5df6873ebap-14,
        0x1.84087b8a73803p-2,
        0x1.bbeb1515e1f29p-15,
        0x1.e237ea6d65d1ap+16,
        0x1.0e09599b76f89p-15,
    },
    // order 8
    {
        0x1.a01a01a01a01ap-16,
        0x1.a01a01a01a347p-16,
        0x1.a01a01a019cedp-16,
        0x1.a01a04a7c5703p-16,
        0x1.a019fe986e9d3p-16,
        0x1.abf4a43a58018p-16,
        0x1.94d363740468fp-16,
        0x1.d36aab6f328b4p-16,
        0x1.761a87da20eeep-16,
        0x1.0ae71af1d13f8p-13,
        0x1.943c0d40c78acp-17,
        0x1.0ae71af1d1d5ap-13,
        0x1.943c0d40c75bap-17,
        0x1.70860dfb243abp-4,
        0x1.bf90e53b31e8fp-18,
        0x1.a10dc3c8fd693p+11,
        0x1.2d7c410ee8030p-18,
    },
};

/// Deviation from the reference, relative, in units of the epsilon of the tested type.
template <typename T>
double ulpDistance(T computed, double reference) {
  if (reference == 0.0) {
    return computed == T{0} ? 0.0 : std::numeric_limits<double>::infinity();
  }
  return std::abs((static_cast<double>(computed) - reference) / reference) /
         static_cast<double>(std::numeric_limits<T>::epsilon());
}

/// Held well above the observed maximum so that a rounding change does not fail the test,
/// and far below the point where any digit of the result would be in doubt.
constexpr double MaxUlp = 16.0;

template <std::size_t N, typename T>
void checkOrder() {
  for (std::size_t i = 0; i < SampleCount; ++i) {
    const double x = ReferenceArguments[N][i];
    const double reference = ReferenceValues[N][i];
    const T value = functions::phi<N, T>(static_cast<T>(x));
    REQUIRE(std::isfinite(static_cast<double>(value)));
    REQUIRE(ulpDistance(value, reference) < MaxUlp);
  }
}

template <typename T, std::size_t... Ns>
void checkOrders(std::index_sequence<Ns...> /*orders*/) {
  (checkOrder<Ns, T>(), ...);
}

template <std::size_t N, typename T>
void checkZero() {
  REQUIRE(functions::phi<N, T>(T{0}) == T{1} / functions::phifunctions::factorial<T>(N));
}

template <typename T, std::size_t... Ns>
void checkZeros(std::index_sequence<Ns...> /*orders*/) {
  (checkZero<Ns, T>(), ...);
}

/**
 * phi_n(x) = x phi_{n+1}(x) + 1/n! ties neighbouring orders together. For negative
 * arguments the two summands cancel, by a factor that the identity itself dictates; the
 * bound admits exactly that factor rather than hiding it behind a loose constant. Keeping
 * it tight is what pins down the range in which phiUpTo may use the recurrence.
 */
template <std::size_t N, typename T>
void checkRecurrence() {
  if constexpr (N > 0) {
    for (int i = -100; i <= 100; ++i) {
      const auto x = static_cast<T>(i * 0.1);
      const T lower = functions::phi<N - 1, T>(x);
      const T summand = x * functions::phi<N, T>(x);
      const T offset = T{1} / functions::phifunctions::factorial<T>(N - 1);
      const double amplification =
          (std::abs(static_cast<double>(summand)) + static_cast<double>(offset)) /
          std::abs(static_cast<double>(lower));
      REQUIRE(ulpDistance(static_cast<T>(summand + offset), static_cast<double>(lower)) <
              MaxUlp * amplification);
    }
  }
}

template <typename T, std::size_t... Ns>
void checkRecurrences(std::index_sequence<Ns...> /*orders*/) {
  (checkRecurrence<Ns, T>(), ...);
}

/**
 * phi_n(x) = 1/(n-1)! int_0^1 exp((1-s) x) s^{n-1} ds is positive for every real argument,
 * however far the evaluation strays from the origin.
 */
template <std::size_t N, typename T>
void checkPositivity() {
  for (double x : {-700.0, -80.0, -37.0, -8.0, -1.0, 0.0, 1.0, 8.0, 37.0, 80.0, 700.0}) {
    const T value = functions::phi<N, T>(static_cast<T>(x));
    if (!std::isfinite(static_cast<double>(value))) {
      continue;
    }
    REQUIRE(value >= T{0});
    // only exp itself decays fast enough to underflow; from the first order on the decay is
    // algebraic, phi_n(x) -> -1/((n-1)! x), and the value stays representable
    if constexpr (N > 0) {
      REQUIRE(value > T{0});
    }
  }
}

template <typename T, std::size_t... Ns>
void checkPositivities(std::index_sequence<Ns...> /*orders*/) {
  (checkPositivity<Ns, T>(), ...);
}

/// Both branches have to agree where the implementation switches between them.
template <std::size_t N, typename T>
void checkBranchAgreement() {
  if constexpr (N > 0) {
    for (int sign : {1, -1}) {
      const auto x = static_cast<T>(sign * static_cast<int>(N + 1));
      const T series =
          functions::phifunctions::phiSeries<N, functions::phifunctions::DefaultSeriesLength<N, T>>(
              x);
      const T difference = (std::exp(x) - functions::phifunctions::truncatedExponential<N, T>(x)) /
                           functions::phifunctions::integerPower(x, N);
      REQUIRE(ulpDistance(series, static_cast<double>(difference)) < MaxUlp);
    }
  }
}

template <typename T, std::size_t... Ns>
void checkBranchAgreements(std::index_sequence<Ns...> /*orders*/) {
  (checkBranchAgreement<Ns, T>(), ...);
}

template <std::size_t N, typename T>
void checkRemainder() {
  for (std::size_t i = 0; i < SampleCount; ++i) {
    const double x = ReferenceArguments[N][i];
    const double reference = ReferenceValues[N][i] * std::pow(x, static_cast<double>(N));
    const T value = functions::expRemainder<N, T>(static_cast<T>(x));
    // the factor x^n drives the remainder out of the exponent range long before phi_n is in
    // any trouble, which is the reason the normalized form is the one to build on
    const bool representable =
        std::abs(reference) > static_cast<double>(std::numeric_limits<T>::min()) &&
        std::isfinite(static_cast<double>(value));
    if (representable) {
      REQUIRE(ulpDistance(value, reference) < MaxUlp);
    }
  }
}

template <typename T, std::size_t... Ns>
void checkRemainders(std::index_sequence<Ns...> /*orders*/) {
  (checkRemainder<Ns, T>(), ...);
}

template <std::size_t N, typename T>
void checkBatch(T x) {
  const auto batch = functions::phiUpTo<N, T>(x);
  REQUIRE(ulpDistance(batch[N], static_cast<double>(functions::phi<N, T>(x))) < MaxUlp);
  if constexpr (N > 0) {
    checkBatch<N - 1, T>(x);
  }
}

} // namespace phifunctionstest

TEST_CASE_TEMPLATE("Phi functions reproduce a high-precision reference" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace phifunctionstest;
  checkOrders<RealT>(std::make_index_sequence<MaxOrder + 1>{});
}

TEST_CASE_TEMPLATE("Phi functions reach their limit at the origin exactly" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace phifunctionstest;
  checkZeros<RealT>(std::make_index_sequence<MaxOrder + 1>{});
}

TEST_CASE_TEMPLATE("Phi functions satisfy the recurrence between neighbouring orders" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace phifunctionstest;
  checkRecurrences<RealT>(std::make_index_sequence<MaxOrder + 1>{});
}

TEST_CASE_TEMPLATE("Phi functions stay positive over the whole range" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace phifunctionstest;
  checkPositivities<RealT>(std::make_index_sequence<MaxOrder + 1>{});
}

TEST_CASE_TEMPLATE("Both phi evaluation branches agree where the implementation switches" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace phifunctionstest;
  checkBranchAgreements<RealT>(std::make_index_sequence<MaxOrder + 1>{});
}

TEST_CASE_TEMPLATE("Exponential series remainders reproduce the reference" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace phifunctionstest;
  checkRemainders<RealT>(std::make_index_sequence<MaxOrder + 1>{});
}

TEST_CASE("The first order remainder is expm1" * doctest::test_suite("numerical")) {
  for (int i = -300; i <= 300; ++i) {
    const double x = i * 0.05;
    REQUIRE(functions::expRemainder<1, double>(x) == std::expm1(x));
  }
}

TEST_CASE_TEMPLATE("The batched phi evaluation agrees with the individual orders" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  // both sides of the switch between the downward recurrence and the direct evaluation
  for (const double x : {-30.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 30.0}) {
    phifunctionstest::checkBatch<8, RealT>(static_cast<RealT>(x));
  }
}

TEST_CASE("The truncated phi series holds within the bound it was sized for" *
          doctest::test_suite("numerical")) {
  constexpr double Bound = 0.25;
  constexpr auto Terms =
      functions::phifunctions::seriesLength(3, Bound, std::numeric_limits<double>::epsilon());
  static_assert(Terms < functions::phifunctions::DefaultSeriesLength<3, double>,
                "a smaller bound has to buy a shorter series");
  for (int i = -250; i <= 250; ++i) {
    const double x = i * (Bound / 250.0);
    REQUIRE(phifunctionstest::ulpDistance(functions::phiTruncated<3, Terms>(x),
                                          functions::phi<3, double>(x)) < phifunctionstest::MaxUlp);
  }
}

} // namespace seissol::unit_test
