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

#ifndef MODULES_H
#define MODULES_H

#include <map>
#include <utility>

#include "utils/logger.h"

#include "Module.h"

namespace seissol
{

/**
 * Possible hooks for modules
 *
 * To add a new hook, add an additional enum entry here, add the corresponding
 * function to {@link Module}, add a specialization of {@link Modules::call(Module*)}
 * and update the function {@link Modules::strHook}.
 *
 * @warning The order of the hooks has to be the same they are called in SeisSol.
 */
enum Hook
{
	PRE_MPI,
	POST_MPI_INIT,
	/** @warning This will be triggered only with C++ mesh readers (not when Gambit3D-mixed is used) */
	POST_MESH,
	PRE_MODEL,
	POST_MODEL,
	FIRST_HOOK = PRE_MPI,
	MAX_INIT_HOOKS = POST_MODEL + 1,
	MAX_HOOKS = POST_MODEL + 1
};

/**
 * Manages SeisSol modules
 */
class Modules
{
private:
	std::multimap<int, Module*> m_hooks[MAX_HOOKS];

	// The hook that should be called next
	Hook m_nextHook;

private:
	Modules()
		: m_nextHook(FIRST_HOOK)
	{
	}

	/**
	 * Do the real work for registering for a hook
	 */
	void _registerHook(Module &module, Hook hook, int priority);

	/**
	 * Do the real work for handling a hook
	 */
	template<Hook hook>
	void _callHook()
	{
		for (std::multimap<int, Module*>::iterator it = m_hooks[hook].begin();
				it != m_hooks[hook].end(); it++) {
			call<hook>(it->second);
		}
	}

private:
	template<Hook hook>
	static void call(Module* module);

	static const char* strHook(Hook hook);

	/**
	 * The only instance of this class
	 *
	 * We need to use a static function to prevent the "static initialization order fiasco".
	 */
	static Modules& instance()
	{
		static Modules _instance;
		return _instance;
	}

// Public interface
public:
	static void registerHook(Module &module, Hook hook, int priority = Module::DEFAULT)
	{
		instance()._registerHook(module, hook, priority);
	}

	template<Hook hook>
	static void callHook()
	{
		instance()._callHook<hook>();
	}

};

// Create all template instances for call
#define MODULES_CALL_INSTANCE(enum, func)                          \
	template<> inline                                          \
	void seissol::Modules::call<seissol::enum>(Module* module) \
	{                                                          \
		module->func();                                    \
	}

MODULES_CALL_INSTANCE(PRE_MPI, preMPI)
MODULES_CALL_INSTANCE(POST_MPI_INIT, postMPIInit)
MODULES_CALL_INSTANCE(POST_MESH, postMesh)
MODULES_CALL_INSTANCE(PRE_MODEL, preModel)
MODULES_CALL_INSTANCE(POST_MODEL, postModel)

}

#endif // MODULES_H
