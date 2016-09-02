/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 */

#include "utils/logger.h"

#include "Backend.h"
#include "posix/Fault.h"
#include "posix/Wavefield.h"
#include "h5/Wavefield.h"
#include "h5/Fault.h"
#include "mpio/Wavefield.h"
#include "mpio/WavefieldAsync.h"
#include "mpio/Fault.h"
#include "mpio/FaultAsync.h"
#ifdef USE_SIONLIB
#include "sionlib/Fault.h"
#include "sionlib/Wavefield.h"
#endif // USE_SIONLIB

void seissol::checkpoint::createBackend(Backend backend, Wavefield* &waveField, Fault* &fault)
{
	switch (backend) {
	case POSIX:
		waveField = new posix::Wavefield();
		fault = new posix::Fault();
		break;
	case HDF5:
		waveField = new h5::Wavefield();
		fault = new h5::Fault();
		break;
	case MPIO:
		waveField = new mpio::Wavefield();
		fault = new mpio::Fault();
		break;
	case MPIO_ASYNC:
		waveField = new mpio::WavefieldAsync();
		fault = new mpio::FaultAsync();
		break;
	case SIONLIB:
#ifdef USE_SIONLIB
		waveField = new sionlib::Wavefield();
		fault = new sionlib::Fault();
		break;
#else //USE_SIONLIB
		logError() << "SIONlib checkpoint backend unsupported";
		break;
#endif //USE_SIONLIB
	default:
		logError() << "Unsupported checkpoint backend";
	}
}