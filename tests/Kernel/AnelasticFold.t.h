// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "TestConfigs.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

namespace seissol::unit_test {

namespace anelastic_fold {

template <typename RealT, std::size_t Size>
void fillRandom(RealT (&values)[Size], std::mt19937& generator) {
  std::uniform_real_distribution<RealT> distribution(-1, 1);
  for (auto& value : values) {
    value = distribution(generator);
  }
}

template <typename RealT, std::size_t Size>
void checkClose(const RealT (&actual)[Size], const RealT (&expected)[Size]) {
  const double tolerance = 1e3 * std::numeric_limits<RealT>::epsilon();
  for (std::size_t i = 0; i < Size; ++i) {
    REQUIRE(std::abs(actual[i] - expected[i]) <=
            tolerance * std::max(1.0, std::abs(static_cast<double>(expected[i]))));
  }
}

} // namespace anelastic_fold

// The anelastic solver adds the flux of an analytical boundary to the extended DOFs. On the host,
// it goes into them before the local integral folds them into the DOFs and the memory variables.
// On the device, the local integral has folded them already: the flux goes into them alone, and is
// folded on its own, with the kernel of the neighbor integral. Both give the same DOFs.
TEST_CASE_TEMPLATE_DEFINE("A nodal flux folded on its own adds what it adds in the local fold" *
                              doctest::test_suite("kernel"),
                          Cfg,
                          AnelasticNodalFluxFold) {
  using real = Real<Cfg>;
  using anelastic_fold::checkClose;
  using anelastic_fold::fillRandom;

  // NOLINTNEXTLINE(bugprone-random-generator-seed)
  std::mt19937 generator(20261104);

  alignas(Alignment) real q[tensor::Q<Cfg>::size()];
  alignas(Alignment) real qane[tensor::Qane<Cfg>::size()];
  alignas(Alignment) real qext[tensor::Qext<Cfg>::size()];
  alignas(Alignment) real iane[tensor::Iane<Cfg>::size()];
  alignas(Alignment) real iNodal[tensor::INodal<Cfg>::size()];
  alignas(Alignment) real aminusT[tensor::AminusT<Cfg>::size()];
  alignas(Alignment) real w[tensor::w<Cfg>::size()];
  alignas(Alignment) real bigW[tensor::W<Cfg>::size()];
  alignas(Alignment) real e[tensor::E<Cfg>::size()];
  fillRandom(q, generator);
  fillRandom(qane, generator);
  fillRandom(qext, generator);
  fillRandom(iane, generator);
  fillRandom(iNodal, generator);
  fillRandom(aminusT, generator);
  fillRandom(w, generator);
  fillRandom(bigW, generator);
  fillRandom(e, generator);

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    CAPTURE(face);

    const auto localFold = [&](real* dofs, real* dofsAne, real* dofsExt) {
      kernel::local<Cfg> krnl;
      krnl.Q = dofs;
      krnl.Qane = dofsAne;
      krnl.Qext = dofsExt;
      krnl.Iane = iane;
      krnl.W = bigW;
      krnl.w = w;
      krnl.E = e;
      krnl.execute();
    };
    const auto nodalFlux = [&](real* dofsExt) {
      kernel::localFluxNodal<Cfg> krnl;
      krnl.bindGlobals(seissol::Pool<Cfg>::host());
      krnl.Qext = dofsExt;
      krnl.INodal = iNodal;
      krnl.AminusT = aminusT;
      krnl.execute(face);
    };

    // the host: the flux, then the fold
    alignas(Alignment) real hostQ[tensor::Q<Cfg>::size()];
    alignas(Alignment) real hostQane[tensor::Qane<Cfg>::size()];
    alignas(Alignment) real hostQext[tensor::Qext<Cfg>::size()];
    std::copy_n(q, tensor::Q<Cfg>::size(), hostQ);
    std::copy_n(qane, tensor::Qane<Cfg>::size(), hostQane);
    std::copy_n(qext, tensor::Qext<Cfg>::size(), hostQext);
    nodalFlux(hostQext);
    localFold(hostQ, hostQane, hostQext);

    // the device: the fold, then the flux alone, folded on its own
    alignas(Alignment) real deviceQ[tensor::Q<Cfg>::size()];
    alignas(Alignment) real deviceQane[tensor::Qane<Cfg>::size()];
    alignas(Alignment) real deviceQext[tensor::Qext<Cfg>::size()];
    std::copy_n(q, tensor::Q<Cfg>::size(), deviceQ);
    std::copy_n(qane, tensor::Qane<Cfg>::size(), deviceQane);
    std::copy_n(qext, tensor::Qext<Cfg>::size(), deviceQext);
    localFold(deviceQ, deviceQane, deviceQext);

    std::fill_n(deviceQext, tensor::Qext<Cfg>::size(), real{0});
    nodalFlux(deviceQext);
    kernel::neighbor<Cfg> fold;
    fold.Q = deviceQ;
    fold.Qane = deviceQane;
    fold.Qext = deviceQext;
    fold.w = w;
    fold.execute();

    checkClose(deviceQ, hostQ);
    checkClose(deviceQane, hostQane);
  }
}

TEST_CASE_TEMPLATE_APPLY(AnelasticNodalFluxFold, ConfigsOfSolver<SolverType::LinearCKAnelastic>);

} // namespace seissol::unit_test
