// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Alignment.h"
#include "Common/Constants.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Solver/MultipleSimulations.h"
#include "TestHelper.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

namespace seissol::unit_test {

// The nodal boundary conditions (free surface with gravity, Dirichlet, analytical) evaluate the
// trace of I at the face nodes and take it back to the volume with project2nFaceTo3m. If the
// boundary values are just that trace, the result has to be the one of the modal local flux.
TEST_CASE("Nodal local flux matches the modal local flux" * doctest::test_suite("kernel")) {
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937 generator(1234);
  std::uniform_real_distribution<real> distribution(-1, 1);

  alignas(Alignment) real dofs[tensor::I::size()]{};
  auto dofsView = init::I::view::create(dofs);
  for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
    auto simulationDofs = multisim::simtensor(dofsView, s);
    for (std::size_t b = 0; b < tensor::I::Shape[multisim::BasisFunctionDimension]; ++b) {
      for (std::size_t q = 0; q < tensor::I::Shape[multisim::BasisFunctionDimension + 1]; ++q) {
        simulationDofs(b, q) = distribution(generator);
      }
    }
  }

  // The same flux solver for both paths; each is written through its own view.
  alignas(Alignment) real fluxPlus[tensor::AplusT::size()]{};
  alignas(Alignment) real fluxMinus[tensor::AminusT::size()]{};
  auto fluxPlusView = init::AplusT::view::create(fluxPlus);
  auto fluxMinusView = init::AminusT::view::create(fluxMinus);
  for (std::size_t i = 0; i < tensor::AplusT::Shape[0]; ++i) {
    for (std::size_t j = 0; j < tensor::AplusT::Shape[1]; ++j) {
      const auto value = distribution(generator);
      if (fluxPlusView.isInRange(i, j)) {
        fluxPlusView(i, j) = value;
      }
      if (fluxMinusView.isInRange(i, j)) {
        fluxMinusView(i, j) = value;
      }
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    alignas(Alignment) real modalQ[tensor::Q::size()]{};
    alignas(Alignment) real nodalQ[tensor::Q::size()]{};
    alignas(Alignment) real trace[tensor::INodal::size()]{};

    kernel::localFlux modalKrnl;
    for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
      modalKrnl.rDivM(i) = init::rDivM::Values[init::rDivM::index(i)];
      modalKrnl.fMrT(i) = init::fMrT::Values[init::fMrT::index(i)];
    }
    modalKrnl.AplusT = fluxPlus;
    modalKrnl.I = dofs;
    modalKrnl.Q = modalQ;
    modalKrnl._prefetch.I = dofs;
    modalKrnl._prefetch.Q = modalQ;
    modalKrnl.execute(face);

    kernel::projectToNodalBoundary projectKrnl;
    for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
      projectKrnl.V3mTo2nFace(i) =
          nodal::init::V3mTo2nFace::Values[nodal::init::V3mTo2nFace::index(i)];
    }
    projectKrnl.I = dofs;
    projectKrnl.INodal = trace;
    projectKrnl.execute(face);

    kernel::localFluxNodal nodalKrnl;
    for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
      nodalKrnl.project2nFaceTo3m(i) =
          init::project2nFaceTo3m::Values[init::project2nFaceTo3m::index(i)];
    }
    nodalKrnl.AminusT = fluxMinus;
    nodalKrnl.INodal = trace;
    nodalKrnl.Q = nodalQ;
    nodalKrnl._prefetch.I = dofs;
    nodalKrnl._prefetch.Q = nodalQ;
    nodalKrnl.execute(face);

    double scale = 0;
    for (const auto value : modalQ) {
      scale = std::max(scale, std::abs(static_cast<double>(value)));
    }
    REQUIRE(scale > 0);

    const double tolerance = 1e3 * std::numeric_limits<real>::epsilon() * scale;
    for (std::size_t j = 0; j < tensor::Q::size(); ++j) {
      REQUIRE(nodalQ[j] == AbsApprox(modalQ[j]).epsilon(tolerance).delta(0));
    }
  }
}

} // namespace seissol::unit_test
