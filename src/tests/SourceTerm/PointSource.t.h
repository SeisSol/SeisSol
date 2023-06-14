#include <SourceTerm/PointSource.h>
#include <SourceTerm/PiecewiseLinearFunction1D.h>

#include <memory>

namespace seissol::unit_test {

TEST_CASE("Transform moment tensor") {
  constexpr double epsilon = 100 * std::numeric_limits<real>::epsilon();

  // strike = dip = rake = pi / 3
  real strike, dip, rake;
  strike = dip = rake = M_PI / 3.0;

  // M_xy = M_yx = 1, others are zero
  real l_localMomentTensorXY[3][3] = {
      {0.0, 1.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
  };
  real l_localSolidVelocityComponent[3] = {0.0, 0.0, 0.0};
  real l_localPressureComponent = 0.0;
  real l_localFluidVelocityComponent[3] = {0.0, 0.0, 0.0};

  auto l_momentTensor = sourceterm::AlignedArray<real, sourceterm::PointSources::TensorSize>{};

  seissol::sourceterm::transformMomentTensor(l_localMomentTensorXY,
                                             l_localSolidVelocityComponent,
                                             l_localPressureComponent,
                                             l_localFluidVelocityComponent,
                                             strike,
                                             dip,
                                             rake,
                                             l_momentTensor);

  // Compare to hand-computed reference solution
  REQUIRE(l_momentTensor[0] == AbsApprox(-5.0 * std::sqrt(3.0) / 32.0).epsilon(epsilon));
  REQUIRE(l_momentTensor[1] == AbsApprox(-7.0 * std::sqrt(3.0) / 32.0).epsilon(epsilon));
  REQUIRE(l_momentTensor[2] == AbsApprox(3.0 * std::sqrt(3.0) / 8.0).epsilon(epsilon));
  REQUIRE(l_momentTensor[3] == AbsApprox(19.0 / 32.0).epsilon(epsilon));
  REQUIRE(l_momentTensor[4] == AbsApprox(-9.0 / 16.0).epsilon(epsilon));
  REQUIRE(l_momentTensor[5] == AbsApprox(-std::sqrt(3.0) / 16.0).epsilon(epsilon));
  REQUIRE(l_momentTensor[6] == 0);
  REQUIRE(l_momentTensor[7] == 0);
  REQUIRE(l_momentTensor[8] == 0);

  // strike = dip = rake = pi / 3
  strike = -1.349886940156521;
  dip = 3.034923466331855;
  rake = 0.725404224946106;

  // Random M
  real l_localMomentTensorXZ[3][3] = {
      {1.833885014595086, -0.970040810572334, 0.602398893453385},
      {-0.970040810572334, -1.307688296305273, 1.572402458710038},
      {0.602398893453385, 1.572402458710038, 2.769437029884877},
  };

  seissol::sourceterm::transformMomentTensor(l_localMomentTensorXZ,
                                             l_localSolidVelocityComponent,
                                             l_localPressureComponent,
                                             l_localFluidVelocityComponent,
                                             strike,
                                             dip,
                                             rake,
                                             l_momentTensor);

  // Compare to hand-computed reference solution
  REQUIRE(l_momentTensor[0] == AbsApprox(-0.415053502680640).epsilon(epsilon));
  REQUIRE(l_momentTensor[1] == AbsApprox(0.648994284092410).epsilon(epsilon));
  REQUIRE(l_momentTensor[2] == AbsApprox(3.061692966762920).epsilon(epsilon));
  REQUIRE(l_momentTensor[3] == AbsApprox(1.909053142737053).epsilon(epsilon));
  REQUIRE(l_momentTensor[4] == AbsApprox(0.677535767462651).epsilon(epsilon));
  REQUIRE(l_momentTensor[5] == AbsApprox(-1.029826812214912).epsilon(epsilon));
  REQUIRE(l_momentTensor[6] == 0.0);
  REQUIRE(l_momentTensor[7] == 0.0);
  REQUIRE(l_momentTensor[8] == 0.0);
}

TEST_CASE("Samples to piecewise linear function 1D") {
  const real l_samples[] = {
      0.312858596637428, -0.864879917324456, -0.030051296196269, -0.164879019209038};
  const unsigned l_numberOfSamples = sizeof(l_samples) / sizeof(real);
  const real l_onsetTime = 1.2;
  const real l_samplingInterval = 0.05;
  auto alloc = std::allocator<real>();

  auto l_pwlf = sourceterm::PiecewiseLinearFunction1D(
      l_samples, l_numberOfSamples, l_onsetTime, l_samplingInterval, alloc);

  for (int i = 0; i < 3; ++i) {
    REQUIRE(l_pwlf.slopes[i] == (l_samples[i + 1] - l_samples[i]) / l_samplingInterval);
    REQUIRE(l_pwlf.intercepts[i] == l_samples[i] - (l_samples[i + 1] - l_samples[i]) /
                                                       l_samplingInterval *
                                                       (l_onsetTime + i * l_samplingInterval));
  }

  REQUIRE(l_pwlf.slopes.size() == 3);
  REQUIRE(l_pwlf.onsetTime == l_onsetTime);
  REQUIRE(l_pwlf.samplingInterval == l_samplingInterval);
}

TEST_CASE("Compute PwLFTimeIntegral") {
  /* Setup corresponds to
   *        |  40*t - 39       if t \in [1.00, 1.05),
   *        | -80*t + 87       if t \in [1.05, 1.10),
   * f(t) = |  60*t - 67       if t \in [1.10, 1.15),
   *        |  10*t - 9.5      if t \in [1.15, 1.20),
   *        |  0               else.
   */
  constexpr double epsilon = 200 * std::numeric_limits<real>::epsilon();
  const real l_samples[] = {1.0, 3.0, -1.0, 2.0, 2.5};
  const unsigned l_numberOfSamples = sizeof(l_samples) / sizeof(real);
  const real l_onsetTime = 1.0;
  const real l_samplingInterval = 0.05;
  auto alloc = std::allocator<real>();
  auto l_pwlf = sourceterm::PiecewiseLinearFunction1D(
      l_samples, l_numberOfSamples, l_onsetTime, l_samplingInterval, alloc);

  // integrate f(t) from -2 to 1.05 (only first term)
  REQUIRE(l_pwlf.timeIntegral(-2.0, 1.05) ==
          AbsApprox(0.5 * 40.0 * (1.05 * 1.05 - 1.0) - 39 * 0.05).epsilon(epsilon));

  // integrate f(t) from 1.04 to 1.06 (over boundary)
  REQUIRE(l_pwlf.timeIntegral(1.04, 1.06) ==
          AbsApprox(0.5 * 40.0 * (1.05 * 1.05 - 1.04 * 1.04) - 39 * 0.01 -
                    0.5 * 80 * (1.06 * 1.06 - 1.05 * 1.05) + 87 * 0.01)
              .epsilon(300 * epsilon));

  // integrate f(t) from 1.10 to 1.10 (on boundary)
  REQUIRE(l_pwlf.timeIntegral(1.1, 1.1) == AbsApprox(0.0).epsilon(epsilon));

  // integrate f(t) from 1.19 to 100 (only last term)
  REQUIRE(l_pwlf.timeIntegral(1.19, 100.0) ==
          AbsApprox(0.5 * 10.0 * (1.2 * 1.2 - 1.19 * 1.19) - 9.5 * 0.01).epsilon(epsilon));

  // integrate f(t) from -100 to 100 (integral over whole support)
  REQUIRE(l_pwlf.timeIntegral(-100.0, 100.0) ==
          AbsApprox(0.5 * 40.0 * (1.05 * 1.05 - 1.00 * 1.00) - 39.0 * 0.05 -
                    0.5 * 80.0 * (1.10 * 1.10 - 1.05 * 1.05) + 87.0 * 0.05 +
                    0.5 * 60.0 * (1.15 * 1.15 - 1.10 * 1.10) - 67.0 * 0.05 +
                    0.5 * 10.0 * (1.20 * 1.20 - 1.15 * 1.15) - 9.5 * 0.05)
              .epsilon(4 * epsilon));
}

} // namespace seissol::unit_test
