/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2017, SeisSol Group
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

#ifndef CHECKPOINT_MPIO_WAVEFIELD_ASYNC_H
#define CHECKPOINT_MPIO_WAVEFIELD_ASYNC_H

#ifndef USE_MPI
#include "Checkpoint/WavefieldDummy.h"
#else // USE_MPI

#include "Wavefield.h"

#endif // USE_MPI

namespace seissol
{

namespace checkpoint
{

namespace mpio
{

#ifndef USE_MPI
typedef WavefieldDummy WavefieldAsync;
#else // USE_MPI

class WavefieldAsync : public Wavefield
{
private:
	/** Buffer for storing the DOFs copy */
	real* m_dofsCopy;

	/** True if a checkpoint was started */
	bool m_started;

public:
	WavefieldAsync()
		: seissol::checkpoint::CheckPoint(IDENTIFIER),
		seissol::checkpoint::Wavefield(IDENTIFIER),
		m_dofsCopy(0L), m_started(false)
	{
	}

	bool init(unsigned long numDofs, unsigned int groupSize = 1);

	void writePrepare(const void* header, size_t headerSize);

	void write(const void* header, size_t headerSize);

	void close();
};

#endif // USE_MPI

}

}

}

#endif // CHECKPOINT_MPIO_WAVEFIELD_ASYNC_H

