// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <array>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

#include "Geometry/Refinement/MeshRefiner.h"
#include "Geometry/Refinement/RefinerUtils.h"
#include "Geometry/Refinement/VariableSubSampler.h"
#include "Kernels/Precision.h"
#include "MockReader.h"

namespace seissol::unit_test {

TEST_CASE("Variable Subsampler") {
  constexpr double Epsilon = std::numeric_limits<real>::epsilon() * 1e1;
  std::srand(1234);

  SUBCASE("Divide by 4") {
    const seissol::refinement::DivideTetrahedronBy4<double> refineBy4;
    const seissol::refinement::VariableSubsampler<double> subsampler(1, refineBy4, 3, 9, 12);

    const std::array<real, 36> expectedDOFs = {
        -0.95909429432054482678,  -0.24576668840548565598, -0.073841666364211855367,
        0.084639342883142787421,  -0.60318547945356015827, -0.29280576995908313975,
        -0.20578768130632243971,  0.23721562803389545371,  -0.14513201479740706068,
        0.18204270118367804621,   -0.66283340897768605604, 0.96121445469738475698,
        0.2514671156611722469,    0.14933892377708751775,  -0.16498441628616278276,
        0.74160144785464598982,   -0.55061439302327430667, -0.25879404893976487578,
        -0.79732770027133659241,  0.19534461286095072818,  0.24679981287805408119,
        0.39225418036479331452,   0.61942505362744282316,  1.3786782367932419735,
        -0.16608198189081485596,  0.1212460819861717054,   0.83965964835616058171,
        0.81367404419197641996,   0.066295701815709762172, 0.25114635060216367046,
        -0.091879653285143331187, 0.64540753085971003244,  -0.11038223602601260342,
        0.49194913342387541766,   0.74918009276514585526,  1.1485201151026944721};

    // For order 3 there are 108 DOFs (taking alignment into account)
    std::array<real, 108> dofs;
    for (int i = 0; i < 108; i++) {
      dofs[i] = static_cast<real>(std::rand()) / static_cast<real>(RAND_MAX);
    }
    unsigned int cellMap[1] = {0};
    // A triangle is divided into four subtriangles there are 9 quantities.
    real outDofs[36];
    std::fill(std::begin(outDofs), std::end(outDofs), 0);

    for (unsigned var = 0; var < 9; var++) {
      subsampler.get(dofs.data(), cellMap, var, &outDofs[var * 4]);
    }
    for (int i = 0; i < 36; i++) {
      REQUIRE(outDofs[i] == AbsApprox(expectedDOFs[i]).epsilon(Epsilon));
    }
  };
}

} // namespace seissol::unit_test
