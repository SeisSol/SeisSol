/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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

#ifndef FREESURFACEWRITER_H
#define FREESURFACEWRITER_H

#include <Geometry/MeshReader.h>
#include <Parallel/MPI.h>
#include <utils/logger.h>
#include <async/Module.h>
#include <Modules/Module.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include "FreeSurfaceWriterExecutor.h"

namespace seissol
{
namespace writer
{

class FreeSurfaceWriter : private async::Module<FreeSurfaceWriterExecutor, FreeSurfaceInitParam, FreeSurfaceParam>, public seissol::Module
{
private:
	/** Is enabled? */
	bool m_enabled;

	/** The asynchronous executor */
	FreeSurfaceWriterExecutor m_executor;
  
  /** free surface integration module. */
  seissol::solver::FreeSurfaceIntegrator* m_freeSurfaceIntegrator;
  
  void constructSurfaceMesh(  MeshReader const& meshReader,
                              unsigned*&        cells,
                              double*&          vertices,
                              unsigned&         nCells,
                              unsigned&         nVertices );

public:
	FreeSurfaceWriter() : m_enabled(false), m_freeSurfaceIntegrator(NULL) {}

	/**
	 * Called by ASYNC on all ranks
	 */
	void setUp()
	{
		setExecutor(m_executor);
	}

	void init(  MeshReader const&                       meshReader,
              seissol::solver::FreeSurfaceIntegrator* freeSurfaceIntegrator,
              char const*                             outputPrefix,
              double                                  interval );

	void write(double time)
	{
		SCOREP_USER_REGION("FreeSurfaceWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION)

		if (!m_enabled) {
      logError() << "Trying to write free surface output, but it is disabled.";
    }

		int const rank = seissol::MPI::mpi.rank();

		wait();

		logInfo(rank) << "Writing free surface at time" << utils::nospace << time << ".";

		FreeSurfaceParam param;
		param.time = time;

		for (unsigned i = 0; i < 2*FREESURFACE_NUMBER_OF_COMPONENTS; ++i) {
			sendBuffer(FreeSurfaceWriterExecutor::VARIABLES0 + i);
    }

		call(param);

		logInfo(rank) << "Writing free surface at time" << utils::nospace << time << ". Done.";
	}

	void close()
	{
		if (!m_enabled) {
			return;
    }

		wait();
	}

	void tearDown()
	{
		if (!m_enabled) {
			return;
    }

		m_executor.finalize();
	}

	//
	// Hooks
	//
	void simulationStart();

	void syncPoint(double currentTime);
};

}

}

#endif // FREESURFACEWRITER_H
