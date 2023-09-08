/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2015, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#include "SeisSol.h"
#include "Initializer/InitProcedure/Init.hpp"
#include "Common/filesystem.h"

#ifdef ACL_DEVICE
#include "device.h"
#endif

int main(int argc, char* argv[]) {
#ifdef ACL_DEVICE
#ifdef USE_MPI
  seissol::MPI::mpi.bindAcceleratorDevice();
#endif // USE_MPI
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->initialize();
  device.api->allocateStackMem();
#endif // ACL_DEVICE

  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("SeisSol");
    LIKWID_MARKER_REGISTER("computeDynamicRuptureFrictionLaw");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePostHook");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePostcomputeImposedState");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePreHook");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePrecomputeStress");
    LIKWID_MARKER_REGISTER("computeDynamicRuptureSpaceTimeInterpolation");
    LIKWID_MARKER_REGISTER("computeDynamicRuptureUpdateFrictionAndSlip");
  }

#pragma omp parallel
  { LIKWID_MARKER_START("SeisSol"); }

  EPIK_TRACER("SeisSol");
  SCOREP_USER_REGION("SeisSol", SCOREP_USER_REGION_TYPE_FUNCTION);

  // Initialize SeisSol
  bool runSeisSol = seissol::SeisSol::main.init(argc, argv);

  const auto stamp = utils::TimeUtils::timeAsString("%F_%T", time(0L));
  seissol::SeisSol::main.setBackupTimeStamp(stamp);

  // Run SeisSol
  if (runSeisSol) {
    seissol::initializer::initprocedure::seissolMain();
  }

#pragma omp parallel
  { LIKWID_MARKER_STOP("SeisSol"); }

  LIKWID_MARKER_CLOSE;
  // Finalize SeisSol
  seissol::SeisSol::main.finalize();

#ifdef ACL_DEVICE
  device.api->finalize();
#endif
  return 0;
}
