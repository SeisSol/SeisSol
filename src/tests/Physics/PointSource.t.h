/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 **/
#include <cxxtest/TestSuite.h>

#include <SourceTerm/PointSource.h>

#if defined(DOUBLE_PRECISION)
#define EPSILON 1e-16
#elif defined(SINGLE_PRECISION)
#define EPSILON 1e-8
#endif

namespace seissol {
  namespace unit_test {
    class PointSourceTestSuite;
  }
}

class seissol::unit_test::PointSourceTestSuite : public CxxTest::TestSuite
{
public:
	void testTransformMomentTensor()
	{
    // strike = dip = rake = pi / 3 
    real strike, dip, rake;
    strike = dip = rake = M_PI / 3.0;

    // M_xy = M_yx = 1, others are zero
    real l_localMomentTensorXY[3][3] = {
      { 0.0, 1.0, 0.0 },
      { 1.0, 0.0, 0.0 },
      { 0.0, 0.0, 0.0 },
    };
    real l_localVelocityComponent[3] = {0.0, 0.0, 0.0};

    real l_momentTensor[NUMBER_OF_QUANTITIES];
    
    seissol::sourceterm::transformMomentTensor(l_localMomentTensorXY, l_localVelocityComponent, strike, dip, rake, l_momentTensor);
    
    // Compare to hand-computed reference solution
    TS_ASSERT_DELTA(l_momentTensor[0], -5.0*sqrt(3.0)/32.0, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[1], -7.0*sqrt(3.0)/32.0, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[2], 3.0*sqrt(3.0)/8.0, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[3], 19.0/32.0, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[4], -9.0/16.0, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[5], -sqrt(3.0)/16.0, 100 * EPSILON);
    TS_ASSERT_EQUALS(l_momentTensor[6], 0.0);
    TS_ASSERT_EQUALS(l_momentTensor[7], 0.0);
    TS_ASSERT_EQUALS(l_momentTensor[8], 0.0);
    
    // strike = dip = rake = pi / 3 
    strike = -1.349886940156521;
    dip = 3.034923466331855;
    rake = 0.725404224946106;

    // Random M
    real l_localMomentTensorXZ[3][3] = {
      {  1.833885014595086,  -0.970040810572334,   0.602398893453385 },
      { -0.970040810572334,  -1.307688296305273,   1.572402458710038 },
      {  0.602398893453385,   1.572402458710038,   2.769437029884877 },
    };
    
    seissol::sourceterm::transformMomentTensor(l_localMomentTensorXZ, l_localVelocityComponent, strike, dip, rake, l_momentTensor);
    
    // Compare to hand-computed reference solution
    TS_ASSERT_DELTA(l_momentTensor[0], -0.415053502680640, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[1], 0.648994284092410, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[2], 3.061692966762920, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[3], 1.909053142737053, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[4], 0.677535767462651, 100 * EPSILON);
    TS_ASSERT_DELTA(l_momentTensor[5], -1.029826812214912, 100 * EPSILON);
    TS_ASSERT_EQUALS(l_momentTensor[6], 0.0);
    TS_ASSERT_EQUALS(l_momentTensor[7], 0.0);
    TS_ASSERT_EQUALS(l_momentTensor[8], 0.0);
	}
  
  void testSamplesToPiecewiseLinearFunction1D()
  {
    real l_samples[] = { 0.312858596637428,  -0.864879917324456,  -0.030051296196269,  -0.164879019209038 };
    unsigned l_numberOfSamples = sizeof(l_samples) / sizeof(real);
    real l_onsetTime = 1.2;
    real l_samplingInterval = 0.05;
    PiecewiseLinearFunction1D l_pwlf;
    
    seissol::sourceterm::samplesToPiecewiseLinearFunction1D(l_samples, l_numberOfSamples, l_onsetTime, l_samplingInterval, &l_pwlf);
    
    for (int i = 0; i < 3; ++i) {
      TS_ASSERT_EQUALS(l_pwlf.slopes[i], (l_samples[i+1] - l_samples[i]) / l_samplingInterval);
      TS_ASSERT_EQUALS(l_pwlf.intercepts[i], l_samples[i] - (l_samples[i+1] - l_samples[i]) / l_samplingInterval * (l_onsetTime + i * l_samplingInterval));
    }
    
    TS_ASSERT_EQUALS(l_pwlf.numberOfPieces, 3);
    TS_ASSERT_EQUALS(l_pwlf.onsetTime, l_onsetTime);
    TS_ASSERT_EQUALS(l_pwlf.samplingInterval, l_samplingInterval);
  }
  
  void testComputePwLFTimeIntegral()
  {
    /* Setup corresponds to
     *        |  40*t - 39       if t \in [1.00, 1.05),
     *        | -80*t + 87       if t \in [1.05, 1.10),
     * f(t) = |  60*t - 67       if t \in [1.10, 1.15),
     *        |  10*t - 9.5      if t \in [1.15, 1.20),
     *        |  0               else.
     */
    real l_samples[] = { 1.0,  3.0,  -1.0,  2.0,  2.5 };
    unsigned l_numberOfSamples = sizeof(l_samples) / sizeof(real);
    real l_onsetTime = 1.0;
    real l_samplingInterval = 0.05;
    PiecewiseLinearFunction1D l_pwlf;
    seissol::sourceterm::samplesToPiecewiseLinearFunction1D(l_samples, l_numberOfSamples, l_onsetTime, l_samplingInterval, &l_pwlf);
    
    // integrate f(t) from -2 to 1.05 (only first term)
    TS_ASSERT_DELTA(seissol::sourceterm::computePwLFTimeIntegral(&l_pwlf, -2.0, 1.05), 0.5*40.0*(1.05*1.05 - 1.0) - 39*0.05, 200 * EPSILON);
    
    // integrate f(t) from 1.04 to 1.06 (over boundary)
    TS_ASSERT_DELTA(seissol::sourceterm::computePwLFTimeIntegral(&l_pwlf, 1.04, 1.06), 0.5*40.0*(1.05*1.05 - 1.04*1.04) - 39*0.01 - 0.5*80*(1.06*1.06 - 1.05*1.05) + 87*0.01, 600 * EPSILON);
    
    // integrate f(t) from 1.10 to 1.10 (on boundary)
    TS_ASSERT_DELTA(seissol::sourceterm::computePwLFTimeIntegral(&l_pwlf, 1.1, 1.1), 0.0, 200 * EPSILON);
    
    // integrate f(t) from 1.19 to 100 (only last term)
    TS_ASSERT_DELTA(seissol::sourceterm::computePwLFTimeIntegral(&l_pwlf, 1.19, 100.0), 0.5*10.0*(1.2*1.2 - 1.19*1.19) - 9.5*0.01, 200 * EPSILON);
    
    // integrate f(t) from -100 to 100 (integral over whole support)
    TS_ASSERT_DELTA(seissol::sourceterm::computePwLFTimeIntegral(&l_pwlf, -100.0, 100.0),
        0.5*40.0*(1.05*1.05 - 1.00*1.00) - 39.0*0.05
      - 0.5*80.0*(1.10*1.10 - 1.05*1.05) + 87.0*0.05
      + 0.5*60.0*(1.15*1.15 - 1.10*1.10) - 67.0*0.05
      + 0.5*10.0*(1.20*1.20 - 1.15*1.15) -  9.5*0.05,
      800 * EPSILON);
  }
  
  void addPointSourceToDOFs()
  {
    /// \todo Write a test if the function's implementation gets non-trivial.
  }
};
