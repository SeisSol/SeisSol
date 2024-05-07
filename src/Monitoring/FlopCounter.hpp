/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 *
 * @section DESCRIPTION
 * Counts the floating point operations in SeisSol.
 **/

#ifndef FLOPCOUNTER_HPP
#define FLOPCOUNTER_HPP

#include <fstream>

// Floating point operations performed in the matrix kernels.
// Remark: These variables are updated by the matrix kernels (subroutine.cpp) only in debug builds.
extern long long libxsmm_num_total_flops;
extern long long pspamm_num_total_flops;

namespace seissol::monitoring {
struct FlopCounter {
  public:
  void init(std::string outputFileNamePrefix);
  void printPerformanceUpdate(double wallTime);
  void printPerformanceSummary(double wallTime);
  void incrementNonZeroFlopsLocal(long long update);
  void incrementHardwareFlopsLocal(long long update);
  void incrementNonZeroFlopsNeighbor(long long update);
  void incrementHardwareFlopsNeighbor(long long update);
  void incrementNonZeroFlopsOther(long long update);
  void incrementHardwareFlopsOther(long long update);
  void incrementNonZeroFlopsDynamicRupture(long long update);
  void incrementHardwareFlopsDynamicRupture(long long update);
  void incrementNonZeroFlopsPlasticity(long long update);
  void incrementHardwareFlopsPlasticity(long long update);

  private:
  std::ofstream out;
  long long previousTotalFlops = 0;
  double previousWallTime = 0;
  // global variables for summing-up SeisSol internal counters
  long long nonZeroFlopsLocal = 0;
  long long hardwareFlopsLocal = 0;
  long long nonZeroFlopsNeighbor = 0;
  long long hardwareFlopsNeighbor = 0;
  long long nonZeroFlopsOther = 0;
  long long hardwareFlopsOther = 0;
  long long nonZeroFlopsDynamicRupture = 0;
  long long hardwareFlopsDynamicRupture = 0;
  long long nonZeroFlopsPlasticity = 0;
  long long hardwareFlopsPlasticity = 0;
};
} // namespace seissol::monitoring

#endif
