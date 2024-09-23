// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_
#define SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "xdmfwriter/XdmfWriter.h"
#include "async/ExecInfo.h"
#include "Monitoring/Stopwatch.h"
#include "Kernels/Precision.h"

namespace seissol
{

namespace writer
{

struct FaultInitParam
{
	static const unsigned int OUTPUT_MASK_SIZE = 20;

	bool outputMask[OUTPUT_MASK_SIZE];
	int timestep;
	xdmfwriter::BackendType backend;
	std::string backupTimeStamp;
};

struct FaultParam
{
	double time;
};

class FaultWriterExecutor
{
public:
	enum BufferIds {
		OUTPUT_PREFIX = 0,
		CELLS = 1,
		VERTICES = 2,
		FAULTTAGS = 3,
		VARIABLES0 = 4
	};

private:
	xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>* m_xdmfWriter;

#ifdef USE_MPI
	/** The MPI communicator for the writer */
	MPI_Comm m_comm;
#endif // USE_MPI

	/** The number of variables that should be written */
	unsigned int m_numVariables;

	/** Backend stopwatch */
	Stopwatch m_stopwatch;

public:
	FaultWriterExecutor()
		: m_xdmfWriter(0L),
#ifdef USE_MPI
		m_comm(MPI_COMM_NULL),
#endif // USE_MPI
		m_numVariables(0)
	{
	}

	/**
	 * Initialize the XDMF writer
	 */
	void execInit(const async::ExecInfo &info, const FaultInitParam &param);

	void exec(const async::ExecInfo &info, const FaultParam &param)
	{
		if (!m_xdmfWriter)
			return;

		m_stopwatch.start();

		m_xdmfWriter->addTimeStep(param.time);

		for (unsigned int i = 0; i < m_numVariables; i++)
			m_xdmfWriter->writeCellData(i, static_cast<const real *>(info.buffer(VARIABLES0 + i)));

		m_xdmfWriter->flush();

		m_stopwatch.pause();
	}

	void setFaultTagsData(const unsigned int *faultTags) {
		m_xdmfWriter->writeExtraIntCellData(faultTags);
	}

	void finalize()
	{
		if (m_xdmfWriter) {
			m_stopwatch.printTime("Time fault writer backend:"
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

    static std::string getLabelName(size_t index) {
	  return LABELS[index];
	}

private:
	/** Variable names in the output */
	static char const * const LABELS[];
};

}

}


#endif // SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_

