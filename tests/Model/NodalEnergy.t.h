// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_NODALENERGY_T_H_
#define SEISSOL_TESTS_MODEL_NODALENERGY_T_H_

// The volume energies of a material that varies inside a cell are integrated
// point by point, each point with the material there. A material that does not
// vary has to give what the cell moments give, and a density that varies has
// to be read at the points it is integrated at: a sample taken for the wrong
// point, or an interpolation that misses a parameter, shows up in either.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Equations/Datastructures.h"
#include "Equations/Energy.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Numerical/Quadrature.h"
#include "ResultWriter/EnergyQuadrature.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string_view>
#include <type_traits>

namespace seissol::unit_test {

namespace nodalenergy {

using Material = seissol::model::MaterialT;
using Energy = seissol::model::EnergyCompute<Material>;
using Quadrature = seissol::writer::EnergyQuadrature<Material>;

constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
/// What a material interpolated to a point may be off by, which follows the
/// precision the interpolation is stored in.
constexpr double MaterialTolerance = std::is_same_v<real, double> ? 1e-10 : 1e-5;

/// A point of a set given as a points x 3 matrix. A row that is zero, such as
/// a vertex at the origin, may be left out of the storage, so every entry is
/// read through the range of the view.
template <typename ViewT>
std::array<double, 3> pointOf(const ViewT& view, std::size_t row) {
  std::array<double, 3> point{};
  for (std::size_t j = 0; j < 3; ++j) {
    if (view.isInRange(row, j)) {
      point[j] = view(row, j);
    }
  }
  return point;
}

/// Random degrees of freedom, and the anelastic ones where the build has them.
struct Dofs {
  alignas(Alignment) std::array<real, tensor::Q::size()> dofs{};
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
  alignas(Alignment) std::array<real, tensor::Qane::size()> dofsAne{};
  [[nodiscard]] const real* anelastic() const { return dofsAne.data(); }
#else
  [[nodiscard]] const real* anelastic() const { return nullptr; }
#endif

  explicit Dofs(std::mt19937& rng) {
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (auto& value : dofs) {
      value = static_cast<real>(gauss(rng));
    }
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
    // the anelastic variables are strain rates, far smaller than the stresses
    for (auto& value : dofsAne) {
      value = static_cast<real>(1e-9 * gauss(rng));
    }
#endif
  }
};

/// What the output gives for a cell whose material is the same everywhere.
std::array<double, Quadrature::Count> fromMoments(const Material& material, const Dofs& dofs) {
  alignas(Alignment) std::array<real, tensor::momentQ::size()> linData{};
  kernel::momentQCompute linKrnl;
  linKrnl.bindGlobals(seissol::Pool::host());
  linKrnl.Q = dofs.dofs.data();
  linKrnl.momentQ = linData.data();
  linKrnl.execute();

  alignas(Alignment) std::array<real, tensor::momentQQ::size()> quadData{};
  kernel::momentQQCompute quadKrnl;
  quadKrnl.bindGlobals(seissol::Pool::host());
  quadKrnl.Q = dofs.dofs.data();
  quadKrnl.momentQQ = quadData.data();
  quadKrnl.execute();

  const auto moments =
      Energy::computeMoments(dofs.dofs.data(), dofs.anelastic(), seissol::Pool::host());
  const auto data = Energy::initEnergyData(material);

  auto lin = init::momentQ::view::create(linData.data());
  auto quad = init::momentQQ::view::create(quadData.data());
  std::array<double, Quadrature::Count> result{};
  for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
    const auto values = Energy::computeEnergies(material,
                                                data,
                                                multisim::simtensor(lin, sim),
                                                multisim::simtensor(quad, sim),
                                                moments,
                                                sim);
    for (std::size_t i = 0; i < values.size(); ++i) {
      result[values.size() * sim + i] = values[i];
    }
  }
  return result;
}

} // namespace nodalenergy

TEST_CASE("Nodal energies of a material that does not vary") {
  using namespace nodalenergy;

  std::mt19937 rng(20260928);
  const Quadrature quadrature;

  for (std::size_t attempt = 0; attempt < 4; ++attempt) {
    const auto material = coefficients::configuredMaterial<Material>(rng);
    std::array<Material, LTS::MaterialNodes> samples;
    samples.fill(material);
    const Dofs dofs(rng);

    const auto expected = fromMoments(material, dofs);
    std::array<double, Quadrature::Points> shearModulus{};
    const auto actual = quadrature.energies(
        samples, dofs.dofs.data(), dofs.anelastic(), seissol::Pool::host(), shearModulus);

    for (std::size_t i = 0; i < expected.size(); ++i) {
      CAPTURE(i);
      REQUIRE(actual[i] == doctest::Approx(expected[i]).epsilon(Tolerance));
    }
    for (const auto modulus : shearModulus) {
      REQUIRE(modulus == doctest::Approx(material.getMuBar()).epsilon(MaterialTolerance));
    }
  }
}

TEST_CASE("Nodal energies read the material where they integrate") {
  using namespace nodalenergy;

  // the density of a poroelastic material enters the kinetic energy together
  // with the fluid, which is beside the point here
  if constexpr (!std::is_base_of_v<seissol::model::PoroElasticMaterial, Material>) {
    std::mt19937 rng(20260929);
    const Quadrature quadrature;

    // a field of degree two, or one where the basis does not reach that,
    // which every sample set carries exactly
    constexpr double Quadratic = ConvergenceOrder > 2 ? 1.0 : 0.0;
    const auto field = [&](const double* point) {
      return 1.0 + 0.3 * point[0] - 0.2 * point[1] + 0.5 * point[2] +
             Quadratic * 0.4 * point[0] * point[2];
    };
    const bool hasShearModulus = Material::ParameterMap.count("mu") > 0;

    const auto material = coefficients::configuredMaterial<Material>(rng);
    std::array<Material, LTS::MaterialNodes> samples;
    const auto nodes = init::materialNodes::view::create(init::materialNodes::Values);
    for (std::size_t node = 0; node < LTS::MaterialNodes; ++node) {
      const auto point = pointOf(nodes, node);
      samples[node] = material;
      samples[node].rho = material.rho * field(point.data());
      if (hasShearModulus) {
        samples[node].*Material::ParameterMap.at("mu") = material.getMuBar() * field(point.data());
      }
    }
    const Dofs dofs(rng);

    std::array<double, Quadrature::Points> shearModulus{};
    const auto actual = quadrature.energies(
        samples, dofs.dofs.data(), dofs.anelastic(), seissol::Pool::host(), shearModulus);

    std::array<std::array<double, 3>, Quadrature::Points> points{};
    double weights[Quadrature::Points]{};
    seissol::quadrature::TetrahedronQuadrature(points.data(), weights, ConvergenceOrder + 1);

    alignas(Alignment) std::array<real, tensor::dofsQP::size()> atPointsData{};
    kernel::evalAtQP evalKrnl;
    evalKrnl.bindGlobals(seissol::Pool::host());
    evalKrnl.Q = dofs.dofs.data();
    evalKrnl.dofsQP = atPointsData.data();
    evalKrnl.execute();
    auto atPoints = init::dofsQP::view::create(atPointsData.data());

    // whichever of the two a material reports its kinetic energy as
    std::array<std::size_t, 2> kinetic{Energy::EnergyCount, Energy::EnergyCount};
    for (std::size_t i = 0; i < Energy::EnergyCount; ++i) {
      if (Energy::Energies[i].name == std::string_view("elastic_kinetic_energy")) {
        kinetic[0] = i;
      }
      if (Energy::Energies[i].name == std::string_view("acoustic_kinetic_energy")) {
        kinetic[1] = i;
      }
    }

    for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
      const auto values = multisim::simtensor(atPoints, sim);
      double expected = 0;
      for (std::size_t point = 0; point < Quadrature::Points; ++point) {
        double velocitySq = 0;
        for (std::size_t i = 0; i < 3; ++i) {
          const double velocity = values(point, Material::VelocityOffset + i);
          velocitySq += velocity * velocity;
        }
        expected += weights[point] * 0.5 * material.rho * field(points[point].data()) * velocitySq;
      }

      double reported = 0;
      for (const auto index : kinetic) {
        if (index < Energy::EnergyCount) {
          reported += actual[Energy::EnergyCount * sim + index];
        }
      }
      REQUIRE(reported == doctest::Approx(expected).epsilon(Tolerance));
    }

    if (hasShearModulus) {
      for (std::size_t point = 0; point < Quadrature::Points; ++point) {
        REQUIRE(shearModulus[point] ==
                doctest::Approx(material.getMuBar() * field(points[point].data()))
                    .epsilon(MaterialTolerance));
      }
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_NODALENERGY_T_H_
