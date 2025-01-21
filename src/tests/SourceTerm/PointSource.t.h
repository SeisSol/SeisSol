// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SourceTerm/PointSource.h"
#include "tests/TestHelper.h"

#include <memory>

namespace seissol::unit_test {

TEST_CASE("Transform moment tensor") {
  constexpr double Epsilon = 100 * std::numeric_limits<real>::epsilon();

  // strike = dip = rake = pi / 3
  real strike = M_PI / 3.0;
  real dip = M_PI / 3.0;
  real rake = M_PI / 3.0;

  // M_xy = M_yx = 1, others are zero
  real localMomentTensorXY[3][3] = {
      {0.0, 1.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
  };
  real localSolidVelocityComponent[3] = {0.0, 0.0, 0.0};
  const real localPressureComponent = 0.0;
  real localFluidVelocityComponent[3] = {0.0, 0.0, 0.0};

  auto momentTensor = seissol::memory::AlignedArray<real, sourceterm::PointSources::TensorSize>{};

  seissol::sourceterm::transformMomentTensor(localMomentTensorXY,
                                             localSolidVelocityComponent,
                                             localPressureComponent,
                                             localFluidVelocityComponent,
                                             strike,
                                             dip,
                                             rake,
                                             momentTensor);

  // Compare to hand-computed reference solution
  REQUIRE(momentTensor[0] == AbsApprox(-5.0 * std::sqrt(3.0) / 32.0).epsilon(Epsilon));
  REQUIRE(momentTensor[1] == AbsApprox(-7.0 * std::sqrt(3.0) / 32.0).epsilon(Epsilon));
  REQUIRE(momentTensor[2] == AbsApprox(3.0 * std::sqrt(3.0) / 8.0).epsilon(Epsilon));
  REQUIRE(momentTensor[3] == AbsApprox(19.0 / 32.0).epsilon(Epsilon));
  REQUIRE(momentTensor[4] == AbsApprox(-9.0 / 16.0).epsilon(Epsilon));
  REQUIRE(momentTensor[5] == AbsApprox(-std::sqrt(3.0) / 16.0).epsilon(Epsilon));
  REQUIRE(momentTensor[6] == 0);
  REQUIRE(momentTensor[7] == 0);
  REQUIRE(momentTensor[8] == 0);

  // strike = dip = rake = pi / 3
  strike = -1.349886940156521;
  dip = 3.034923466331855;
  rake = 0.725404224946106;

  // Random M
  real localMomentTensorXZ[3][3] = {
      {1.833885014595086, -0.970040810572334, 0.602398893453385},
      {-0.970040810572334, -1.307688296305273, 1.572402458710038},
      {0.602398893453385, 1.572402458710038, 2.769437029884877},
  };

  seissol::sourceterm::transformMomentTensor(localMomentTensorXZ,
                                             localSolidVelocityComponent,
                                             localPressureComponent,
                                             localFluidVelocityComponent,
                                             strike,
                                             dip,
                                             rake,
                                             momentTensor);

  // Compare to hand-computed reference solution
  REQUIRE(momentTensor[0] == AbsApprox(-0.415053502680640).epsilon(Epsilon));
  REQUIRE(momentTensor[1] == AbsApprox(0.648994284092410).epsilon(Epsilon));
  REQUIRE(momentTensor[2] == AbsApprox(3.061692966762920).epsilon(Epsilon));
  REQUIRE(momentTensor[3] == AbsApprox(1.909053142737053).epsilon(Epsilon));
  REQUIRE(momentTensor[4] == AbsApprox(0.677535767462651).epsilon(Epsilon));
  REQUIRE(momentTensor[5] == AbsApprox(-1.029826812214912).epsilon(Epsilon));
  REQUIRE(momentTensor[6] == 0.0);
  REQUIRE(momentTensor[7] == 0.0);
  REQUIRE(momentTensor[8] == 0.0);
}

} // namespace seissol::unit_test
