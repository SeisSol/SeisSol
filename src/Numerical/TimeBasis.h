// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_NUMERICAL_TIMEBASIS_H_
#define SEISSOL_SRC_NUMERICAL_TIMEBASIS_H_

#include "Functions.h"
#include "Quadrature.h"

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>
namespace seissol::numerical {

/*
Abstracted time basis class. Provides coefficient vectors for integration, point evaluation and
derivatives.
*/
template <typename RealT>
class TimeBasis {
  public:
  explicit TimeBasis(std::size_t order) : order_(order) {}
  virtual ~TimeBasis() = default;

  [[nodiscard]] std::size_t order() const { return order_; }

  [[nodiscard]] virtual std::vector<RealT> derivative(double position, double timestep) const = 0;
  [[nodiscard]] virtual std::vector<RealT> point(double position, double timestep) const = 0;
  [[nodiscard]] virtual std::vector<RealT>
      integrate(double start, double end, double timestep) const = 0;

  /*
    Nodes and weights of a quadrature rule over one timestep, as (nodes,
    weights). Gauss-Legendre with as many nodes as the basis has functions,
    which integrates a polynomial of degree 2n-1 exactly.

    A solver whose flux is nonlinear in the state cannot integrate it in
    closed form and samples it instead. Where it samples is a question about
    the timestep rather than about the basis, so the rule lives here and not in
    the kernel that happens to need it first; the coefficients that read an
    expansion at its nodes are paired with it one level up.
  */
  [[nodiscard]] std::pair<std::vector<double>, std::vector<double>>
      quadrature(double timestep) const {
    std::vector<double> reference(order_);
    std::vector<double> referenceWeights(order_);
    seissol::quadrature::GaussJacobi(reference.data(), referenceWeights.data(), order_, 0, 0);

    // The rule comes back over [-1, 1] and in descending order. Time runs
    // forwards here: a predictor that carries an internal variable from one
    // node to the next needs them that way round.
    std::vector<double> nodes(order_);
    std::vector<double> weights(order_);
    for (std::size_t i = 0; i < order_; ++i) {
      const std::size_t source = order_ - 1 - i;
      nodes[i] = 0.5 * (reference[source] + 1.0) * timestep;
      weights[i] = 0.5 * referenceWeights[source] * timestep;
    }
    return {nodes, weights};
  }

  /*
    The same, with the ends of the timestep among the nodes: Gauss-Lobatto
    with one node more than the basis has functions.

    A rule whose nodes are all interior leaves the interval from the start of
    the step to its first node unsampled, which a predictor carrying an
    internal variable across the nodes cannot account for. Lobatto has the
    ends, so the chain of nodes covers the step without a gap and the first
    node is the state itself.

    n Lobatto nodes integrate a polynomial of degree 2n-3 exactly, so n =
    order + 1 matches what order Gauss-Legendre nodes do -- the endpoints cost
    one node, not accuracy.
  */
  [[nodiscard]] std::pair<std::vector<double>, std::vector<double>>
      quadratureWithEndpoints(double timestep) const {
    const std::size_t points = order_ + 1;
    std::vector<double> nodes(points);
    std::vector<double> weights(points);

    // The interior nodes are the roots of the derivative of the Legendre
    // polynomial of degree points-1, which are the Gauss-Jacobi nodes for
    // a = b = 1; the weights follow from that polynomial at the nodes.
    if (points > 2) {
      std::vector<double> interior(points - 2);
      std::vector<double> unused(points - 2);
      seissol::quadrature::GaussJacobi(interior.data(), unused.data(), points - 2, 1, 1);
      for (std::size_t i = 0; i < points - 2; ++i) {
        nodes[i + 1] = interior[points - 3 - i];
      }
    }
    nodes[0] = -1.0;
    nodes[points - 1] = 1.0;

    const double scale = 2.0 / static_cast<double>(points * (points - 1));
    for (std::size_t i = 0; i < points; ++i) {
      const double legendre = seissol::functions::JacobiP(points - 1, 0, 0, nodes[i]);
      weights[i] = scale / (legendre * legendre);
      // map [-1, 1] onto the timestep, keeping the order
      weights[i] *= 0.5 * timestep;
      nodes[i] = 0.5 * (nodes[i] + 1.0) * timestep;
    }
    return {nodes, weights};
  }

  private:
  std::size_t order_;
};

/*
A time basis in the form

f_i(x) = x^i / i! .

Used in the Cauchy-Kovalevskaya kernels.

*/
template <typename RealT>
class MonomialBasis : public TimeBasis<RealT> {
  public:
  ~MonomialBasis() override = default;
  explicit MonomialBasis(std::size_t order) : TimeBasis<RealT>(order), order_(order) {}

  [[nodiscard]] std::vector<RealT> derivative(double position, double /*timestep*/) const override {
    std::vector<RealT> coeffs(order_);
    coeffs[0] = 0;
    if (coeffs.size() > 1) {
      coeffs[1] = 1;
      double coeffCache = 1;
      for (std::size_t i = 1; i + 1 < order_; ++i) {
        coeffCache *= position / i;
        coeffs[i + 1] = coeffCache;
      }
    }
    return coeffs;
  }

  [[nodiscard]] std::vector<RealT> point(double position, double /*timestep*/) const override {
    std::vector<RealT> coeffs(order_);
    coeffs[0] = 1;
    double coeffCache = 1;
    for (std::size_t i = 1; i < order_; ++i) {
      coeffCache *= position / i;
      coeffs[i] = coeffCache;
    }
    return coeffs;
  }

  [[nodiscard]] std::vector<RealT>
      integrate(double start, double end, double /*timestep*/) const override {
    std::vector<RealT> coeffs(order_);
    double coeffStart = start;
    double coeffEnd = end;
    for (std::size_t i = 0; i < order_; ++i) {
      coeffs[i] = coeffEnd - coeffStart;
      coeffStart *= start / (i + 2);
      coeffEnd *= end / (i + 2);
    }
    return coeffs;
  }

  private:
  std::size_t order_;
};

/**
Legendre polynomial time basis.

Used in the Space-Time Predictor kernels.
*/
template <typename RealT>
class LegendreBasis : public TimeBasis<RealT> {
  public:
  ~LegendreBasis() override = default;
  explicit LegendreBasis(std::size_t order) : TimeBasis<RealT>(order), order_(order) {}

  [[nodiscard]] std::vector<RealT> derivative(double position, double timestep) const override {
    const double tau = position / timestep;
    std::vector<RealT> data(order_);
    for (std::size_t i = 0; i < order_; ++i) {
      data[i] = seissol::functions::shiftedLegendre(i, tau, 1) / timestep;
    }
    return data;
  }

  [[nodiscard]] std::vector<RealT> point(double position, double timestep) const override {
    const double tau = position / timestep;
    std::vector<RealT> data(order_);
    for (std::size_t i = 0; i < order_; ++i) {
      data[i] = seissol::functions::shiftedLegendre(i, tau, 0);
    }
    return data;
  }

  [[nodiscard]] std::vector<RealT>
      integrate(double start, double end, double timestep) const override {
    const double tauS = start / timestep;
    const double tauE = end / timestep;
    std::vector<RealT> data(order_);
    for (std::size_t i = 0; i < order_; ++i) {
      // apply integral transform with x |-> (x * timestep)

      const auto fE = seissol::functions::shiftedLegendre(i, tauE, -1);
      const auto fS = seissol::functions::shiftedLegendre(i, tauS, -1);
      data[i] = timestep * (fE - fS);
    }
    return data;
  }

  private:
  std::size_t order_;
};

/**
A basis for every expansion a solver keeps.

A solver that transports only its state keeps one, and the two sets of
coefficients are the same. One that transports more keeps the second in
another basis: the coefficients of the state fall out of a recursion, the rest
have to be won from samples, and the basis that is well conditioned for
winning them is not the one the recursion produces.

The pairing lives here rather than in the coefficients, so that a solver names
one basis and a time operation is one call. Where the two bases coincide the
work is done once.

`CoefficientsT` is whatever the caller wants filled; it needs a `state` and an
`extra` member that can be written through iterators. `QuadratureT` is what a
quadrature is returned in; it is built from the nodes, their weights, and one
`CoefficientsT` per node, in that order.
*/
template <typename CoefficientsT,
          typename QuadratureT,
          typename StateBasisT,
          typename ExtraBasisT = StateBasisT>
class CompoundTimeBasis {
  public:
  explicit CompoundTimeBasis(std::size_t order) : order_(order), state_(order), extra_(order) {}

  [[nodiscard]] std::size_t order() const { return order_; }

  [[nodiscard]] CoefficientsT derivative(double position, double timestep) const {
    return both([&](const auto& basis) { return basis.derivative(position, timestep); });
  }

  [[nodiscard]] CoefficientsT point(double position, double timestep) const {
    return both([&](const auto& basis) { return basis.point(position, timestep); });
  }

  [[nodiscard]] CoefficientsT integrate(double start, double end, double timestep) const {
    return both([&](const auto& basis) { return basis.integrate(start, end, timestep); });
  }

  /*
    A quadrature over one step: where to sample it, with what weight, and the
    coefficients that read every expansion at each of its nodes.

    A solver whose flux is nonlinear in the state cannot integrate that flux in
    closed form and samples it instead, so all three travel together -- and
    pairing a rule with the coefficients of both expansions is the same job
    this class does for the time operations. Doing it here is what keeps a
    quadrature one call, and keeps whoever samples from having to evaluate a
    basis at its own nodes.

    Which nodes a rule has is a question about the timestep rather than about a
    basis, so either member answers that part alike; the state's does.
  */
  [[nodiscard]] QuadratureT quadrature(double timestep) const {
    return collocate(state_.quadrature(timestep), timestep);
  }

  [[nodiscard]] QuadratureT quadratureWithEndpoints(double timestep) const {
    return collocate(state_.quadratureWithEndpoints(timestep), timestep);
  }

  private:
  /// A rule, with every expansion evaluated at each of its nodes.
  [[nodiscard]] QuadratureT collocate(std::pair<std::vector<double>, std::vector<double>> rule,
                                      double timestep) const {
    QuadratureT result{std::move(rule.first), std::move(rule.second), {}};
    result.coefficients.reserve(result.nodes.size());
    for (const auto& node : result.nodes) {
      result.coefficients.push_back(point(node, timestep));
    }
    return result;
  }

  template <typename Operation>
  [[nodiscard]] CoefficientsT both(Operation&& operation) const {
    CoefficientsT coefficients{};
    const auto state = operation(state_);
    std::copy(state.begin(), state.end(), coefficients.state.begin());
    if constexpr (std::is_same_v<StateBasisT, ExtraBasisT>) {
      std::copy(state.begin(), state.end(), coefficients.extra.begin());
    } else {
      const auto extra = operation(extra_);
      std::copy(extra.begin(), extra.end(), coefficients.extra.begin());
    }
    return coefficients;
  }

  std::size_t order_;
  StateBasisT state_;
  ExtraBasisT extra_;
};

} // namespace seissol::numerical
#endif // SEISSOL_SRC_NUMERICAL_TIMEBASIS_H_
