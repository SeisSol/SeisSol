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

#ifndef WAVE_FIELD_WRITER_H
#define WAVE_FIELD_WRITER_H

#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#include "utils/logger.h"

#include "async/Module.h"

#include "Checkpoint/DynStruct.h"
#include "Geometry/refinement/VariableSubSampler.h"
#include "Monitoring/Stopwatch.h"
#include "WaveFieldWriterExecutor.h"
#include <Modules/Module.h>

namespace seissol
{
namespace refinement {
  template<typename T>
  class MeshRefiner;
}

namespace writer
{

class WaveFieldWriter : private async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>, public seissol::Module
{
	/** True if wave field output is enabled */
	bool m_enabled;

	/** The timestep component in the checkpoint header */
	DynStruct::Component<int> m_timestepComp;

	/** False if entire region is to be written */
	bool m_extractRegion;

	/** The asynchronous executor */
	WaveFieldWriterExecutor m_executor;

	/** Variable buffer ids (high and low order variables) */
	int m_variableBufferIds[2];

	/** The output prefix for the filename */
	std::string m_outputPrefix;

	/** The variable subsampler for the refined mesh */
	refinement::VariableSubsampler<double>* m_variableSubsampler;

	/** Number of variables */
	unsigned int m_numVariables;

	/** Number of integrated variables */
	unsigned int m_numIntegratedVariables;

	/** Flag indicated which variables should be written */
	bool* m_outputFlags;

	/** Flag indicated which low variables should be written */
	bool* m_lowOutputFlags;

	/** Refined number of cells */
	unsigned int m_numCells;

	/** Unrefined (low order) number of cells */
	unsigned int m_numLowCells;

	/** Pointer to the degrees of freedom */
	const real* m_dofs;

	/** Pointer to the plastic strain */
	const real* m_pstrain;

	/** Pointer to the integrals */
	const real* m_integrals;

	/** Mapping from the cell order to dofs order */
	unsigned int* m_map;

	/** The stopwatch for the frontend */
	Stopwatch m_stopwatch;

	/** Checks if a vertex given by the vertexCoords lies inside the boxBounds */
	/*   The boxBounds is in the format: xMin, xMax, yMin, yMax, zMin, zMax */
	bool vertexInBox(const double * const boxBounds, const double * const vertexCoords) {
		if (vertexCoords[0] <= boxBounds[1] && vertexCoords[0] >= boxBounds[0] &&
				vertexCoords[1] <= boxBounds[3] && vertexCoords[1] >= boxBounds[2] &&
				vertexCoords[2] <= boxBounds[5] && vertexCoords[2] >= boxBounds[4]) {
			return true;
		} else {
			return false;
		}
	}
  
  refinement::TetrahedronRefiner<double>* createRefiner(int refinement);
  
  unsigned const* adjustOffsets(refinement::MeshRefiner<double>* meshRefiner);

public:
	WaveFieldWriter()
		: m_enabled(false),
		  m_extractRegion(false),
		  m_variableSubsampler(0L),
		  m_numVariables(0),
		  m_outputFlags(0L),
		  m_lowOutputFlags(0L),
		  m_numCells(0), m_numLowCells(0),
		  m_dofs(0L), m_pstrain(0L), m_integrals(0L),
		  m_map(0L)
	{
	}

	/**
	 * Activate the wave field output
	 */
	void enable();

	/**
	 * @return True if wave field output is enabled, false otherwise
	 */
	bool isEnabled() const
	{
		return m_enabled;
	}

	/**
	 * Set the output prefix for the filename
	 */
	void setFilename(const char* outputPrefix)
	{
		m_outputPrefix = outputPrefix;
	}

	/**
	 * Called by ASYNC on all ranks
	 */
	void setUp();

  void setWaveFieldInterval(double interval) {
    setSyncInterval(interval);
  }

	/**
	 * Initialize the wave field ouput
	 *
	 * @param map The mapping from the cell order to dofs order
	 * @param timeTolerance The tolerance in the time for ignoring duplicate time steps
	 */
	void init(unsigned int numVars, int order, int numAlignedDOF,
			const MeshReader &meshReader,  const std::vector<unsigned> &LtsClusteringData,
			const real* dofs,  const real* pstrain, const real* integrals,
			unsigned int* map,
			int refinement, int* outputMask, double* outputRegionBounds,
      xdmfwriter::BackendType backend);

	/**
	 * Write a time step
	 */
	void write(double time);

	/**
	 * Close wave field writer and free resources
	 */
	void close()
	{
		// Cleanup the executor
		if (m_enabled)
			wait();

		finalize();

		if (!m_enabled)
			return;

		m_stopwatch.printTime("Time wave field writer frontend:");

		delete m_variableSubsampler;
		m_variableSubsampler = 0L;
		delete [] m_outputFlags;
		m_outputFlags = 0L;
		delete [] m_lowOutputFlags;
		m_lowOutputFlags = 0L;
		if (m_extractRegion) {
			delete [] m_map;
			m_map = 0L;
		}
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

#endif // WAVE_FIELD_WRITER_H
