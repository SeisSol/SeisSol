// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "Geometry/Refinement/MeshRefiner.h"
#include "Geometry/Refinement/RefinerUtils.h"
#include "Geometry/Refinement/VariableSubSampler.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "MockReader.h"
#include "Numerical/BasisFunction.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <random>

namespace seissol::unit_test {

TEST_CASE("Variable Subsampler") {
  constexpr double Epsilon = std::numeric_limits<real>::epsilon() * 1e1;

  // NOLINTNEXTLINE (-cert-dcl59-cpp)
  std::mt19937 rnggen(1234);

  SUBCASE("Divide by 4") {
    constexpr std::size_t Alignment = 4;
    constexpr std::size_t ConvergenceOrder = 3;
    constexpr std::size_t Quantities = 9;
    constexpr std::size_t AlignedQuantities =
        kernels::getNumberOfAlignedBasisFunctions(ConvergenceOrder, Alignment);
    constexpr std::size_t SubTriangles = 4;

    const seissol::refinement::DivideTetrahedronBy4<double> refineBy4;
    const seissol::refinement::VariableSubsampler<double> subsampler(
        1, refineBy4, ConvergenceOrder, Quantities, AlignedQuantities, false);

    std::array<real, Quantities * SubTriangles> expectedDOFs{};

    std::uniform_real_distribution<> rngdist(0.0, 1.0);

    // For order 3 there are 108 DOFs (taking the alignment of 4 into account)
    std::array<real, AlignedQuantities * Quantities> dofs{};
    for (std::size_t i = 0; i < dofs.size(); ++i) {
      dofs[i] = rngdist(rnggen);
    }

    // compute the expected result (and the subdivision) manually
    std::vector<basisFunction::SampledBasisFunctions<double>> sbf;
    sbf.emplace_back(ConvergenceOrder, 5. / 16, 1. / 16, 5. / 16);
    sbf.emplace_back(ConvergenceOrder, 1. / 16, 5. / 16, 5. / 16);
    sbf.emplace_back(ConvergenceOrder, 5. / 16, 5. / 16, 1. / 16);
    sbf.emplace_back(ConvergenceOrder, 5. / 16, 5. / 16, 5. / 16);

    for (std::size_t i = 0; i < Quantities; ++i) {
      for (std::size_t j = 0; j < SubTriangles; ++j) {
        expectedDOFs[i * SubTriangles + j] = sbf[j].evalWithCoeffs(&dofs[i * AlignedQuantities]);
      }
    }

    // now compare to the actual implementation

    std::array<unsigned int, 1> cellMap{};

    // A triangle is divided into four subtriangles; there are 9 quantities.
    std::array<real, Quantities * SubTriangles> outDofs{};

    for (std::size_t var = 0; var < Quantities; ++var) {
      subsampler.get(dofs.data(), cellMap.data(), var, &outDofs[var * SubTriangles]);
    }
    for (std::size_t i = 0; i < Quantities * SubTriangles; ++i) {
      REQUIRE(outDofs[i] == AbsApprox(expectedDOFs[i]).epsilon(Epsilon));
    }
  };
}

} // namespace seissol::unit_test
