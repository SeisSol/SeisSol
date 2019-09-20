/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <Geometry/MeshReader.h>
#include <utils/logger.h>
#include <async/Module.h>
#include <Modules/Module.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include "Checkpoint/DynStruct.h"
#include "Monitoring/Stopwatch.h"
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

	/** Timestep component in the checkpoint header */
	DynStruct::Component<int> m_timestepComp;

	/** The asynchronous executor */
	FreeSurfaceWriterExecutor m_executor;

	/** Frontend stopwatch */
	Stopwatch m_stopwatch;

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
		if (isAffinityNecessary()) {
		  const auto freeCpus = parallel::getFreeCPUsMask();
		  logInfo(seissol::MPI::mpi.rank()) << "Free surface writer thread affinity:" << parallel::maskToString(parallel::getFreeCPUsMask());
		  if (parallel::freeCPUsMaskEmpty(freeCpus)) {
		    logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
		  }
		  setAffinityIfNecessary(freeCpus);
		}
	}

	void enable();

	void init(  MeshReader const&                       meshReader,
              seissol::solver::FreeSurfaceIntegrator* freeSurfaceIntegrator,
              char const*                             outputPrefix,
              double                                  interval,
              xdmfwriter::BackendType                 backend );

	void write(double time);

	void close()
	{
		if (m_enabled)
			wait();

		finalize();

		if (!m_enabled)
			return;

		m_stopwatch.printTime("Time free surface writer frontend:");
	}

	void tearDown()
	{
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
