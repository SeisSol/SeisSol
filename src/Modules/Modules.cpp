/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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

#include <cassert>

#include "Modules.h"

void seissol::Modules::_registerHook(Module &module, Hook hook, int priority)
{
	assert(hook < MAX_HOOKS);

	if (m_nextHook >= MAX_INIT_HOOKS)
		logError() << "Trying to register for a hook after initialization phase";
	if (hook < m_nextHook)
		logError() << "Trying to register for hook" << strHook(hook)
			<< "but SeisSol was already processing" << strHook(static_cast<Hook>(m_nextHook-1));

	m_hooks[hook].insert(std::pair<int, Module*>(priority, &module));
}

const char* seissol::Modules::strHook(Hook hook)
{
	switch (hook) {
	case PRE_MPI:
		return "PRE_MPI";
	case POST_MPI_INIT:
		return "POST_MPI_INIT";
	case POST_MESH:
		return "POST_MESH";
	case PRE_MODEL:
		return "PRE_MODEL";
	case POST_MODEL:
		return "POST_MODEL";
	case SIMULATION_START:
		return "SIMULATION_START";
	case SYNCHRONIZATION_POINT:
		return "SYNCHRONIZATION_POINT";
	default:
		return "unknown hook";
	}
}
