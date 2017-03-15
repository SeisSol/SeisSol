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

#ifndef CHECKPOINT_MPIO_WAVEFIELD_H
#define CHECKPOINT_MPIO_WAVEFIELD_H

#ifndef USE_MPI
#include "Checkpoint/WavefieldDummy.h"
#else // USE_MPI

#include <mpi.h>

#include "CheckPoint.h"
#include "Checkpoint/Wavefield.h"
#include "Checkpoint/DynStruct.h"

#endif // USE_MPI

namespace seissol
{

namespace checkpoint
{

namespace mpio
{

#ifndef USE_MPI
typedef WavefieldDummy Wavefield;
#else // USE_MPI

class Wavefield : public CheckPoint, virtual public seissol::checkpoint::Wavefield
{
private:
	/** The partition component in the header */
	DynStruct::Component<int> m_partitionComp;

	/** MPI-IO supports buffer larger > 2 GB */
	bool m_useLargeBuffer;

public:
	Wavefield()
		: seissol::checkpoint::CheckPoint(IDENTIFIER),
		seissol::checkpoint::Wavefield(IDENTIFIER),
		CheckPoint(IDENTIFIER),
		m_useLargeBuffer(true)
	{
	}

	void setHeader(WavefieldHeader &header);

	bool init(size_t headerSize, unsigned long numDofs, unsigned int groupSize = 1);

	void load(real* dofs);

	void initHeader(WavefieldHeader &header);

	void write(const void* header, size_t headerSize);

protected:
	bool validate(MPI_File file);

	/**
	 * Write the header information to the file
	 *
	 * @param time
	 * @param timestepWaveField
	 */
	void writeHeader(const void* header, size_t headerSize);

protected:
	static const unsigned long IDENTIFIER = 0x7A3B4;
};

#endif // USE_MPI

}

}

}

#endif // CHECKPOINT_MPIO_WAVEFIELD_H

