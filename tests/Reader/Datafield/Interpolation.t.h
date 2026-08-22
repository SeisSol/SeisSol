// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Cross-checks the sampler against a transcription of easi's scalar
// Grid<Derived>::mapWith driver (davschneller/precommit). The oracle is the
// point of this file: the two code paths read the same .nc files in the same
// run, and a divergence between them is a wrong simulation, not a failed unit
// test. Everything else here is the easi tests/interpolation.cpp battery,
// carried over so that the ported kernel keeps the properties the original had.

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::reader::datafield;

namespace datafield_test {

// --------------------------------------------------------------------------
// Backing store: values in grid-axis order, axis 0 fastest, components last.
// --------------------------------------------------------------------------
struct Store {
  std::vector<double> f64;
  std::vector<float> f32;
  ArrayView view;
  std::vector<double> minv, deltav;
  std::vector<unsigned> numv;
};

Store makeStore(const std::vector<double>& min,
                const std::vector<double>& delta,
                const std::vector<unsigned>& num,
                unsigned comps,
                const std::function<double(const double*, unsigned)>& f) {
  Store s;
  s.minv = min;
  s.deltav = delta;
  s.numv = num;
  const unsigned dim = static_cast<unsigned>(num.size());
  std::size_t total = 1;
  std::vector<std::size_t> stride(dim);
  for (unsigned d = 0; d < dim; ++d) {
    stride[d] = total;
    total *= num[d];
  }
  s.f64.resize(total * comps);
  std::vector<double> pos(dim);
  for (std::size_t flat = 0; flat < total; ++flat) {
    std::size_t rest = flat;
    for (unsigned d = 0; d < dim; ++d) {
      pos[d] = min[d] + static_cast<double>(rest % num[d]) * delta[d];
      rest /= num[d];
    }
    for (unsigned c = 0; c < comps; ++c) {
      s.f64[flat * comps + c] = f(pos.data(), c);
    }
  }
  s.f32.assign(s.f64.begin(), s.f64.end());

  s.view.data = s.f64.data();
  s.view.type = ElementType::Double;
  s.view.dimensions = dim;
  s.view.components = comps;
  s.view.componentStride = 1;
  for (unsigned d = 0; d < dim; ++d) {
    s.view.min[d] = min[d];
    s.view.deltaInv[d] = 1.0 / delta[d];
    s.view.num[d] = num[d];
    s.view.stride[d] = stride[d] * comps;
  }
  return s;
}

// --------------------------------------------------------------------------
// Oracle: easi's scalar driver, AoS in / AoS out.
// --------------------------------------------------------------------------
void oracle(const ArrayView& v, Interpolation type, const double* x, std::size_t n, double* out) {
  const StencilKernel k = kernelOf(type);
  const unsigned width = k.width;
  const unsigned dim = v.dimensions;
  const unsigned comps = v.components;
  std::size_t stencilSize = 1;
  for (unsigned d = 0; d < dim; ++d) {
    stencilSize *= width;
  }
  std::vector<double> stencil(stencilSize * comps), work(stencilSize * comps);

  for (std::size_t i = 0; i < n; ++i) {
    std::vector<double> w(dim * MaxStencilWidth, 0.0);
    std::vector<long> start(dim);
    for (unsigned d = 0; d < dim; ++d) {
      const long num = static_cast<long>(v.num[d]);
      const double top = static_cast<double>(num - 1);
      double raw = (x[i * dim + d] - v.min[d]) * v.deltaInv[d];
      if (!(raw > 0.0)) {
        raw = 0.0;
      } else if (raw > top) {
        raw = top;
      }
      double* wd = w.data() + d * MaxStencilWidth;
      if (k.roundToNearest) {
        start[d] = std::lround(raw);
        wd[0] = 1.0;
        continue;
      }
      const long maxBase = (num >= 2) ? num - 2 : 0;
      const double lower = std::floor(raw);
      long index = static_cast<long>(lower);
      double s = 0.0;
      if (index > maxBase) {
        index = maxBase;
        s = (num >= 2) ? 1.0 : 0.0;
      } else {
        s = raw - lower;
      }
      if (width <= 2) {
        start[d] = index;
        interpolationWeights(type, s, wd);
      } else {
        const long maxStart = num - static_cast<long>(width);
        long ws = (num >= static_cast<long>(width)) ? index + k.offset : 0;
        if (ws < 0) {
          ws = 0;
        } else if (ws > maxStart) {
          ws = maxStart > 0 ? maxStart : 0;
        }
        const unsigned slot = static_cast<unsigned>(index - ws);
        if (num >= static_cast<long>(width) && slot == static_cast<unsigned>(-k.offset)) {
          interpolationWeights(type, s, wd);
        } else {
          linearFallbackWeights(width, slot, s, wd);
        }
        start[d] = ws;
      }
    }

    std::vector<unsigned> si(dim, 0);
    for (std::size_t f = 0; f < stencilSize; ++f) {
      std::size_t off = 0;
      for (unsigned d = 0; d < dim; ++d) {
        long node = start[d] + static_cast<long>(si[d]);
        const long limit = static_cast<long>(v.num[d]) - 1;
        if (node > limit) {
          node = limit;
        }
        off += static_cast<std::size_t>(node) * v.stride[d];
      }
      for (unsigned c = 0; c < comps; ++c) {
        const std::size_t e = off + c * v.componentStride;
        stencil[f * comps + c] = (v.type == ElementType::Float)
                                     ? static_cast<double>(static_cast<const float*>(v.data)[e])
                                     : static_cast<const double*>(v.data)[e];
      }
      for (unsigned d = 0; d < dim; ++d) {
        if (++si[d] < width) {
          break;
        }
        si[d] = 0;
      }
    }

    const double* src = stencil.data();
    std::size_t cnt = stencilSize;
    for (int d = static_cast<int>(dim) - 1; d >= 0; --d) {
      cnt /= width;
      const double* wd = w.data() + static_cast<unsigned>(d) * MaxStencilWidth;
      for (std::size_t p = 0; p < cnt; ++p) {
        for (unsigned c = 0; c < comps; ++c) {
          double acc = 0.0;
          for (unsigned kk = 0; kk < width; ++kk) {
            acc += wd[kk] * src[(p + kk * cnt) * comps + c];
          }
          work[p * comps + c] = acc;
        }
      }
      src = work.data();
    }
    for (unsigned c = 0; c < comps; ++c) {
      out[i * comps + c] = work[c];
    }
  }
}

// The kernel takes a pointer per axis and per component now (see
// Interpolation.h); these tests were written against the flat blocks and the
// two describe the same memory, so the adapter is address arithmetic only. The
// pointer-array behaviour that the flat form CANNOT express -- skipping a
// component, aliasing out onto x -- gets its own cases below.
template <typename Out>
std::size_t sampleBatchFlat(const ArrayView& v,
                            Interpolation t,
                            const double* x,
                            std::size_t n,
                            Out* out,
                            SampleScratch& sc) {
  std::vector<const double*> xs(v.dimensions);
  std::vector<Out*> outs(v.components);
  for (unsigned d = 0; d < v.dimensions; ++d) {
    xs[d] = x + static_cast<std::size_t>(d) * n;
  }
  for (unsigned c = 0; c < v.components; ++c) {
    outs[c] = out + static_cast<std::size_t>(c) * n;
  }
  return sampleBatch(v, t, xs.data(), n, outs.data(), sc);
}

// AoS helper over the SoA kernel
std::vector<double>
    interp(const ArrayView& v, Interpolation t, const std::vector<std::vector<double>>& pts) {
  const std::size_t n = pts.size();
  std::vector<double> xs(n * v.dimensions), ys(n * v.components), r(n * v.components);
  for (std::size_t i = 0; i < n; ++i) {
    for (unsigned d = 0; d < v.dimensions; ++d) {
      xs[d * n + i] = pts[i][d];
    }
  }
  SampleScratch sc;
  sampleBatchFlat(v, t, xs.data(), n, ys.data(), sc);
  for (std::size_t i = 0; i < n; ++i) {
    for (unsigned c = 0; c < v.components; ++c) {
      r[i * v.components + c] = ys[c * n + i];
    }
  }
  return r;
}

const char* nameOf(Interpolation t) {
  switch (t) {
  case Interpolation::Nearest:
    return "nearest";
  case Interpolation::Linear:
    return "linear";
  case Interpolation::Cubic:
    return "cubic";
  }
  return "?";
}

} // namespace datafield_test

TEST_CASE("Datafield interpolation matches the easi reference") {
  using namespace datafield_test;

  SUBCASE("cubic weights are Keys a = -1/2 and sum to one") {
    const std::vector<double> nodes = {3.0, -1.0, 4.0, 1.0, -5.0, 9.0};
    auto store = makeStore({0.0}, {1.0}, {6}, 1, [&](const double* p, unsigned) {
      return nodes[static_cast<std::size_t>(std::lround(p[0]))];
    });
    const double expected = -1.0 / 16.0 * nodes[1] + 9.0 / 16.0 * nodes[2] + 9.0 / 16.0 * nodes[3] -
                            1.0 / 16.0 * nodes[4];
    REQUIRE(interp(store.view, Interpolation::Cubic, {{2.5}})[0] ==
            doctest::Approx(expected).epsilon(1e-14));
    for (int i = 0; i <= 20; ++i) {
      double w[MaxStencilWidth];
      interpolationWeights(Interpolation::Cubic, i / 20.0, w);
      REQUIRE(std::accumulate(w, w + 4, 0.0) == doctest::Approx(1.0).epsilon(1e-14));
    }
  }

  SUBCASE("every scheme reproduces the node values") {
    auto f = [](const double* p, unsigned) { return std::sin(1.7 * p[0]) + 0.3 * p[0]; };
    for (auto type : {Interpolation::Nearest, Interpolation::Linear, Interpolation::Cubic}) {
      auto store = makeStore({0.0}, {0.25}, {12}, 1, f);
      for (unsigned i = 0; i < 12; ++i) {
        const double x = 0.25 * i;
        REQUIRE(interp(store.view, type, {{x}})[0] == doctest::Approx(f(&x, 0)).epsilon(1e-13));
      }
    }
  }

  SUBCASE("polynomial reproduction is exact up to the scheme order, and no further") {
    struct Case {
      Interpolation type;
      unsigned maxDegree;
    };
    for (const Case& c : {Case{Interpolation::Linear, 1}, Case{Interpolation::Cubic, 2}}) {
      for (unsigned degree = 0; degree <= 3; ++degree) {
        auto store = makeStore({-2.0}, {0.5}, {16}, 1, [degree](const double* p, unsigned) {
          return std::pow(p[0], degree);
        });
        double maxError = 0.0;
        for (int i = 0; i <= 40; ++i) {
          const double x = -1.0 + 4.0 * i / 40.0;
          maxError = std::max(
              maxError, std::fabs(interp(store.view, c.type, {{x}})[0] - std::pow(x, degree)));
        }
        if (degree <= c.maxDegree) {
          REQUIRE(maxError < 1e-12);
        } else {
          REQUIRE(maxError > 1e-6);
        }
      }
    }
  }

  SUBCASE("convergence order is second for linear and third for cubic") {
    auto f = [](const double* p, unsigned) { return std::sin(3.1 * p[0]) * std::exp(0.7 * p[1]); };
    struct Case {
      Interpolation type;
      double order;
    };
    for (const Case& c : {Case{Interpolation::Linear, 2.0}, Case{Interpolation::Cubic, 3.0}}) {
      std::vector<double> errors;
      for (unsigned n : {33u, 65u, 129u}) {
        const double h = 1.0 / (n - 1);
        auto store = makeStore({0.0, 0.0}, {h, h}, {n, n}, 1, f);
        std::vector<std::vector<double>> points;
        for (int i = 0; i < 25; ++i) {
          for (int j = 0; j < 25; ++j) {
            points.push_back({0.2 + 0.6 * i / 24.0, 0.2 + 0.6 * j / 24.0});
          }
        }
        const auto y = interp(store.view, c.type, points);
        double maxError = 0.0;
        for (std::size_t k = 0; k < points.size(); ++k) {
          maxError = std::max(maxError, std::fabs(y[k] - f(points[k].data(), 0)));
        }
        errors.push_back(maxError);
      }
      const double order = std::log2(errors[errors.size() - 2] / errors.back());
      REQUIRE(std::fabs(order - c.order) < 0.3);
    }
  }

  SUBCASE("the stencil degrades to linear near an edge instead of extrapolating") {
    auto f = [](const double* p, unsigned) { return std::sin(1.3 * p[0]) + 0.5 * p[0] * p[0]; };
    auto store = makeStore({0.0}, {0.5}, {10}, 1, f);
    for (double x : {0.05, 0.25, 0.49, 4.05, 4.3, 4.49}) {
      REQUIRE(interp(store.view, Interpolation::Cubic, {{x}})[0] ==
              doctest::Approx(interp(store.view, Interpolation::Linear, {{x}})[0]).epsilon(1e-15));
    }
    REQUIRE(std::fabs(interp(store.view, Interpolation::Cubic, {{0.75}})[0] -
                      interp(store.view, Interpolation::Linear, {{0.75}})[0]) > 1e-9);
  }

  SUBCASE("out-of-domain and non-finite coordinates clamp instead of indexing wild") {
    auto store =
        makeStore({0.0}, {1.0}, {5}, 1, [](const double* p, unsigned) { return 2.0 * p[0] + 1.0; });
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE(interp(store.view, Interpolation::Cubic, {{-10.0}})[0] == doctest::Approx(1.0));
    REQUIRE(interp(store.view, Interpolation::Cubic, {{10.0}})[0] == doctest::Approx(9.0));
    REQUIRE(interp(store.view, Interpolation::Cubic, {{nan}})[0] == doctest::Approx(1.0));
    REQUIRE(interp(store.view, Interpolation::Cubic, {{-inf}})[0] == doctest::Approx(1.0));
    REQUIRE(interp(store.view, Interpolation::Cubic, {{inf}})[0] == doctest::Approx(9.0));

    const double xs[] = {-3.0, 1.5, 2.5, 99.0, nan};
    double out[5];
    SampleScratch scratch;
    REQUIRE(sampleBatchFlat(store.view, Interpolation::Linear, xs, 5, out, scratch) == 3);
  }

  SUBCASE("a degenerate axis always samples its single slice") {
    auto store = makeStore(
        {0.0, 5.0}, {1.0, 1.0}, {6, 1}, 1, [](const double* p, unsigned) { return 3.0 + p[0]; });
    store.view.deltaInv[1] = 0.0;
    REQUIRE(interp(store.view, Interpolation::Cubic, {{2.0, 5.0}})[0] == doctest::Approx(5.0));
    REQUIRE(interp(store.view, Interpolation::Cubic, {{2.0, 42.0}})[0] == doctest::Approx(5.0));
  }

  SUBCASE("results do not depend on lane order or batch size") {
    auto f = [](const double* p, unsigned) { return std::sin(2.0 * p[0]) * std::cos(1.4 * p[1]); };
    auto store = makeStore({0.0, 0.0}, {0.25, 0.25}, {20, 20}, 1, f);
    std::vector<std::vector<double>> clustered;
    std::vector<std::vector<double>> scrambled;
    for (int i = 0; i < 64; ++i) {
      clustered.push_back({1.05 + 0.001 * (i % 8), 2.05 + 0.001 * (i / 8)});
    }
    for (int i = 0; i < 64; ++i) {
      scrambled.push_back(clustered[(i * 37) % 64]);
    }
    const auto a = interp(store.view, Interpolation::Cubic, clustered);
    const auto b = interp(store.view, Interpolation::Cubic, scrambled);
    for (int i = 0; i < 64; ++i) {
      REQUIRE(a[(i * 37) % 64] == b[i]);
      REQUIRE(interp(store.view, Interpolation::Cubic, {clustered[i]})[0] == a[i]);
    }
  }

  SUBCASE("the batch kernel agrees with the easi oracle bit for bit") {
    std::mt19937_64 rng(9001);
    std::uniform_real_distribution<double> uniform(-1.5, 6.5);
    auto f = [](const double* p, unsigned c) {
      return std::sin(1.9 * p[0] + 0.3) * std::cos(0.7 * p[1]) + 0.11 * p[2] + 0.5 * c;
    };
    for (unsigned dim = 1; dim <= 3; ++dim) {
      std::vector<double> min(dim, -1.0);
      std::vector<double> delta(dim);
      std::vector<unsigned> num(dim);
      for (unsigned d = 0; d < dim; ++d) {
        delta[d] = 0.3 + 0.05 * d;
        num[d] = 9 + 2 * d;
      }
      auto store = makeStore(min, delta, num, 2, [&](const double* p, unsigned c) {
        const double q[3] = {p[0], dim > 1 ? p[1] : 0.0, dim > 2 ? p[2] : 0.0};
        return f(q, c);
      });
      const std::size_t n = 257; // not a multiple of any vector width
      std::vector<double> aos(n * dim);
      std::vector<double> soa(n * dim);
      for (std::size_t i = 0; i < n; ++i) {
        for (unsigned d = 0; d < dim; ++d) {
          const double x = uniform(rng);
          aos[i * dim + d] = x;
          soa[d * n + i] = x;
        }
      }
      for (auto type : {Interpolation::Nearest, Interpolation::Linear, Interpolation::Cubic}) {
        std::vector<double> reference(n * 2);
        std::vector<double> actual(n * 2);
        oracle(store.view, type, aos.data(), n, reference.data());
        SampleScratch scratch;
        sampleBatchFlat(store.view, type, soa.data(), n, actual.data(), scratch);
        std::size_t differing = 0;
        for (std::size_t i = 0; i < n; ++i) {
          for (unsigned c = 0; c < 2; ++c) {
            if (!(reference[i * 2 + c] == actual[c * n + i])) {
              ++differing;
            }
          }
        }
        REQUIRE(differing == 0);
      }
    }
  }

  SUBCASE("axis permutation is pure metadata") {
    // The same values laid out in file order (axis 0 slowest) must sample
    // identically once stride/num/min/deltaInv are permuted. If this ever fails,
    // the reversal has leaked into the hot path and a transposed velocity model
    // has somewhere to hide.
    auto f = [](const double* p, unsigned c) {
      return 1.0 + p[0] + 2.0 * p[1] * p[1] - 0.5 * p[2] + 0.25 * c;
    };
    const std::vector<unsigned> num = {9, 7, 11};
    const unsigned components = 2;
    auto natural = makeStore({-1.0, 0.5, 2.0}, {0.4, 0.25, 0.3}, num, components, f);

    std::vector<double> reversed(natural.f64.size());
    std::vector<std::size_t> reversedStride(3);
    std::size_t total = 1;
    for (int d = 2; d >= 0; --d) {
      reversedStride[d] = total;
      total *= num[d];
    }
    for (unsigned i = 0; i < num[0]; ++i) {
      for (unsigned j = 0; j < num[1]; ++j) {
        for (unsigned k = 0; k < num[2]; ++k) {
          const std::size_t from = i + num[0] * (j + static_cast<std::size_t>(num[1]) * k);
          const std::size_t to =
              i * reversedStride[0] + j * reversedStride[1] + k * reversedStride[2];
          for (unsigned c = 0; c < components; ++c) {
            reversed[to * components + c] = natural.f64[from * components + c];
          }
        }
      }
    }
    ArrayView view = natural.view;
    view.data = reversed.data();
    for (unsigned d = 0; d < 3; ++d) {
      view.stride[d] = reversedStride[d] * components;
    }

    std::mt19937_64 rng(31337);
    std::uniform_real_distribution<double> uniform(-2.0, 5.0);
    const std::size_t n = 128;
    std::vector<double> xs(n * 3);
    for (auto& e : xs) {
      e = uniform(rng);
    }
    for (auto type : {Interpolation::Nearest, Interpolation::Linear, Interpolation::Cubic}) {
      std::vector<double> a(n * components);
      std::vector<double> b(n * components);
      SampleScratch s1;
      SampleScratch s2;
      sampleBatchFlat(natural.view, type, xs.data(), n, a.data(), s1);
      sampleBatchFlat(view, type, xs.data(), n, b.data(), s2);
      for (std::size_t i = 0; i < n * components; ++i) {
        REQUIRE(a[i] == b[i]);
      }
    }
  }

  SUBCASE("f32 storage tracks f64 to float precision") {
    auto f = [](const double* p, unsigned) { return std::sin(1.1 * p[0]) + 0.3 * p[1]; };
    auto store = makeStore({0.0, 0.0}, {0.5, 0.5}, {12, 10}, 1, f);
    ArrayView view32 = store.view;
    view32.data = store.f32.data();
    view32.type = ElementType::Float;
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> uniform(0.0, 4.0);
    std::vector<std::vector<double>> points;
    for (int i = 0; i < 64; ++i) {
      points.push_back({uniform(rng), uniform(rng)});
    }
    const auto a = interp(store.view, Interpolation::Cubic, points);
    const auto b = interp(view32, Interpolation::Cubic, points);
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
      worst = std::max(worst, std::fabs(a[i] - b[i]));
    }
    REQUIRE(worst > 0.0);
    REQUIRE(worst < 1e-6);
  }
}

TEST_CASE("Datafield sampling honours the pointer-array contract") {
  using namespace datafield_test;

  // 3-D, three components, so the skip and the aliasing cases are not degenerate.
  const auto store =
      makeStore({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {6, 6, 6}, 3, [](const double* p, unsigned c) {
        return (c + 1) * (p[0] + 2.0 * p[1] + 3.0 * p[2]) + c;
      });

  constexpr std::size_t N = 7;
  std::vector<double> xs(3 * N);
  for (std::size_t i = 0; i < N; ++i) {
    xs[0 * N + i] = 0.3 + 0.6 * static_cast<double>(i);
    xs[1 * N + i] = 1.1 + 0.4 * static_cast<double>(i);
    xs[2 * N + i] = 2.2 + 0.3 * static_cast<double>(i);
  }
  const double* coords[3] = {&xs[0], &xs[N], &xs[2 * N]};

  for (const auto type : {Interpolation::Nearest, Interpolation::Linear, Interpolation::Cubic}) {
    // Reference: every component requested.
    std::vector<double> all(3 * N, 0.0);
    double* allPtr[3] = {&all[0], &all[N], &all[2 * N]};
    SampleScratch scratch;
    sampleBatch(store.view, type, coords, N, allPtr, scratch);

    SUBCASE("a skipped component does not change the ones that are kept") {
      // Skipping must be exactly that -- not a different code path that happens
      // to agree to within a tolerance. Bitwise, or the null-skip is a second
      // numerical path and Grid.h's cross-path argument no longer holds.
      std::vector<double> only(N, -1.0);
      double* onlyPtr[3] = {nullptr, &only[0], nullptr};
      SampleScratch s2;
      sampleBatch(store.view, type, coords, N, onlyPtr, s2);
      for (std::size_t i = 0; i < N; ++i) {
        REQUIRE(only[i] == all[N + i]);
      }
    }

    SUBCASE("out may alias x, because the slot allocator can make it so") {
      std::vector<double> aliased = xs;
      const double* aliasCoords[3] = {&aliased[0], &aliased[N], &aliased[2 * N]};
      // Component 1 lands on top of the y coordinate run it was computed from.
      double* aliasOut[3] = {nullptr, &aliased[N], nullptr};
      SampleScratch s3;
      sampleBatch(store.view, type, aliasCoords, N, aliasOut, s3);
      for (std::size_t i = 0; i < N; ++i) {
        REQUIRE(aliased[N + i] == all[N + i]);
      }
    }

    SUBCASE("the float overload differs from the f64 one only by the store") {
      std::vector<float> f(3 * N, 0.0F);
      float* fPtr[3] = {&f[0], &f[N], &f[2 * N]};
      SampleScratch s4;
      sampleBatch(store.view, type, coords, N, fPtr, s4);
      for (std::size_t k = 0; k < 3 * N; ++k) {
        REQUIRE(f[k] == static_cast<float>(all[k]));
      }
    }

    SUBCASE("out entries need not be contiguous with each other") {
      // The whole point of the pointer array: the destinations are separate
      // transient slots, which the allocator does not place adjacently.
      std::vector<double> a(N, 0.0);
      std::vector<double> b(N, 0.0);
      std::vector<double> c(N, 0.0);
      double* scattered[3] = {a.data(), b.data(), c.data()};
      SampleScratch s5;
      sampleBatch(store.view, type, coords, N, scattered, s5);
      for (std::size_t i = 0; i < N; ++i) {
        REQUIRE(a[i] == all[i]);
        REQUIRE(b[i] == all[N + i]);
        REQUIRE(c[i] == all[2 * N + i]);
      }
    }
  }
}

} // namespace seissol::unit_test
