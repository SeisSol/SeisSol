/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Gilbert Brietzke (gilbert.brietzke AT lrz.de, http://www.lrz.de)
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

#ifndef CHECKPOINT_SIONLIB_WAVEFIELD_H
#define CHECKPOINT_SIONLIB_WAVEFIELD_H

#include "CheckPoint.h"
#include "Checkpoint/Wavefield.h"

namespace seissol
{

namespace checkpoint
{

namespace sionlib
{

class Wavefield : public CheckPoint, virtual public seissol::checkpoint::Wavefield
{
public:
	Wavefield()
		: seissol::checkpoint::CheckPoint(IDENTIFIER),
		seissol::checkpoint::Wavefield(IDENTIFIER),
		CheckPoint(IDENTIFIER)
	{}

	bool init(size_t headerSize, unsigned long numDofs, unsigned int groupSize = 1);

	void load(real* dofs);

	void write(const void* header, size_t headerSize);

private:
	static const unsigned long IDENTIFIER = 0x7A389;
};

}

}

}

#endif // CHECKPOINT_SIONLIB_WAVEFIELD_H

