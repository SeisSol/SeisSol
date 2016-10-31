/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2016, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Main C++ SeisSol file
 */

#ifndef SEISSOL_H
#define SEISSOL_H

#include <string>

#include "utils/logger.h"

#ifdef GENERATEDKERNELS
#include "Solver/time_stepping/TimeManager.h"
#include "Solver/Simulator.h"
#include "Initializer/time_stepping/LtsLayout.h"
#include "Checkpoint/Manager.h"
#include "SourceTerm/Manager.h"
#endif // GENERATEDKERNELS

#include "ResultWriter/AsyncIO.h"
#include "ResultWriter/WaveFieldWriter.h"
#include "ResultWriter/PostProcessor.h"

#define SEISSOL_VERSION_STRING "SVN Mainline"

class MeshReader;

namespace seissol
{

/**
 * @todo Initialize rank
 */
class SeisSol
{
private:
	/** The name of the parameter file */
	std::string m_parameterFile;

	/** Async I/O handler (needs to be initialize before other I/O modules) */
	io::AsyncIO m_asyncIO;

	MeshReader* m_meshReader;

#ifdef GENERATEDKERNELS
	/*
	 * initializers
	 */
	initializers::time_stepping::LtsLayout m_ltsLayout;

	initializers::MemoryManager m_memoryManager;

	//! time manager
	time_stepping::TimeManager  m_timeManager;

	//! simulator
	Simulator m_simulator;

	/** Check pointing module */
	checkpoint::Manager m_checkPointManager;

	/** Source term module */
	sourceterm::Manager m_sourceTermManager;
#endif // GENERATEDKERNELS

	/** Wavefield output module */
	writer::WaveFieldWriter m_waveFieldWriter;

	/** PostProcessor module **/
	writer::PostProcessor m_postProcessor;

private:
	/**
	 * Only one instance of this class should exist (private constructor).
	 */
	SeisSol()
		: m_meshReader(0L)
	{}

public:
	/**
	 * Cleanup data structures
	 */
	virtual ~SeisSol()
	{
		delete m_meshReader;
	}

	/**
	 * Initialize C++ part of the program
	 */
	bool init(int argc, char* argv[]);

	/**
	 * Finalize SeisSol
	 */
	void finalize();

	const char* parameterFile() const
	{
		return m_parameterFile.c_str();
	}

#ifdef GENERATEDKERNELS
	initializers::time_stepping::LtsLayout& getLtsLayout(){ return m_ltsLayout; }

	initializers::MemoryManager& getMemoryManager() {
    return m_memoryManager;
  }

	time_stepping::TimeManager& timeManager()
	{
		return m_timeManager;
	}

	Simulator& simulator() { return m_simulator; }

	checkpoint::Manager& checkPointManager()
	{
		return m_checkPointManager;
	}

	sourceterm::Manager& sourceTermManager()
	{
		return m_sourceTermManager;
	}
#endif // GENERATEDKERNELS

	io::AsyncIO& asyncIO()
	{
		return m_asyncIO;
	}

	/**
	 * Get the wave field writer module
	 */
	writer::WaveFieldWriter& waveFieldWriter()
	{
		return m_waveFieldWriter;
	}

	/** Get the post processor module
	 */
	 writer::PostProcessor& postProcessor()
	 {
		 return m_postProcessor;
	 }

	/**
	 * Set the mesh reader
	 */
	void setMeshReader(MeshReader* meshReader)
	{
		if (m_meshReader != 0L)
			logError() << "Mesh reader already initialized";

		m_meshReader = meshReader;
	}

	/**
	 * Delete the mesh reader to free memory resources.
	 *
	 * Should be called after initialization
	 */
	void freeMeshReader()
	{
		delete m_meshReader;
		m_meshReader = 0L;
	}

	/**
	 * Get the mesh reader
	 */
	const MeshReader& meshReader() const
	{
		return *m_meshReader;
	}

	/**
	 * Get the mesh reader
	 */
	MeshReader& meshReader()
	{
		return *m_meshReader;
	}

public:
	/** The only instance of this class; the main C++ functionality */
	static SeisSol main;
};

}

#endif // SEISSOL_H
