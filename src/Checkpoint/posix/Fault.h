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

#ifndef CHECKPOINT_POSIX_FAULT_H
#define CHECKPOINT_POSIX_FAULT_H

#include "CheckPoint.h"
#include "Checkpoint/Fault.h"

namespace seissol
{

namespace checkpoint
{

namespace posix
{

class Fault : public CheckPoint, virtual public seissol::checkpoint::Fault
{
public:
	Fault()
		: seissol::checkpoint::CheckPoint(IDENTIFIER),
		seissol::checkpoint::Fault(IDENTIFIER),
		CheckPoint(IDENTIFIER, sizeof(int))
	{}

	bool init(unsigned int numSides, unsigned int numBndGP,
		unsigned int groupSize = 1);

	/**
	 * @param[out] timestepFault Time step of the fault writer in the checkpoint
	 *  (if the fault writer was active)
	 */
	void load(int &timestepFault, double* mu, double* slipRate1, double* slipRate2,
		double* slip, double* slip1, double* slip2, double* state, double* strength);

	void write(int timestepFault);

	void close()
	{
		if (numSides() == 0)
			return;

		CheckPoint::close();
	}

private:
	static const unsigned long IDENTIFIER = 0x7A127;
};

}

}

}

#endif // CHECKPOINT_POSIX_FAULT_H
