// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_
#define SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_

#include "xdmfwriter/XdmfWriter.h"
#include "async/ExecInfo.h"

#include "Monitoring/Stopwatch.h"

namespace seissol
{
namespace writer
{
struct FreeSurfaceInitParam
{
	int timestep;
	xdmfwriter::BackendType backend;
	std::string backupTimeStamp;
};

struct FreeSurfaceParam
{
	double time;
};

class FreeSurfaceWriterExecutor
{
public:
	enum BufferIds {
		OUTPUT_PREFIX = 0,
		CELLS = 1,
		VERTICES = 2,
		LOCATIONFLAGS = 3,
		VARIABLES0 = 4,
	};

private:
#ifdef USE_MPI
	/** The MPI communicator for the writer */
	MPI_Comm m_comm;
#endif // USE_MPI

	xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>* m_xdmfWriter;
  unsigned m_numVariables;

	/** Backend stopwatch */
	Stopwatch m_stopwatch;

public:
	FreeSurfaceWriterExecutor()
		:
#ifdef USE_MPI
		m_comm(MPI_COMM_NULL),
#endif // USE_MPI
		m_xdmfWriter(0L),
		m_numVariables(0) {}

	/**
	 * Initialize the XDMF writer
	 */
	void execInit(const async::ExecInfo &info, const FreeSurfaceInitParam &param);

	void exec(const async::ExecInfo &info, const FreeSurfaceParam &param)
	{
		if (!m_xdmfWriter) {
			return;
		}

		m_stopwatch.start();

		m_xdmfWriter->addTimeStep(param.time);

		for (unsigned int i = 0; i < m_numVariables; i++) {
			m_xdmfWriter->writeCellData(i, static_cast<const real*>(info.buffer(VARIABLES0 + i)));
    }

		m_xdmfWriter->flush();

		m_stopwatch.pause();
	}

	void setLocationFlagData(const unsigned int *locationFlags) {
		m_xdmfWriter->writeExtraIntCellData(locationFlags);
	}

	void finalize()
	{
		if (m_xdmfWriter) {
			m_stopwatch.printTime("Time free surface writer backend:"
#ifdef USE_MPI
				, m_comm
#endif // USE_MPI
			);
		}

#ifdef USE_MPI
		if (m_comm != MPI_COMM_NULL) {
			MPI_Comm_free(&m_comm);
			m_comm = MPI_COMM_NULL;
		}
#endif // USE_MPI

		delete m_xdmfWriter;
		m_xdmfWriter = 0L;
	}

private:
	/** Variable names in the output */
	static char const * const LABELS[];
};

}

}


#endif // SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_

