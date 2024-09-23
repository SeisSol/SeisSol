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

#include "Parallel/MPI.h"

#include <string>
#include <vector>

#include "Solver/FreeSurfaceIntegrator.h"
#include "utils/logger.h"
#include "FreeSurfaceWriterExecutor.h"

/**
 * Initialize the XDMF writers
 */
void seissol::writer::FreeSurfaceWriterExecutor::execInit(const async::ExecInfo &info, const seissol::writer::FreeSurfaceInitParam &param)
{
	if (m_xdmfWriter) {
		logError() << "Free surface writer already initialized.";
  }

	unsigned int nCells = info.bufferSize(CELLS) / (3 * sizeof(int));
	unsigned int nVertices = info.bufferSize(VERTICES) / (3 * sizeof(double));

#ifdef USE_MPI
	MPI_Comm_split(seissol::MPI::mpi.comm(), (nCells > 0 ? 0 : MPI_UNDEFINED), 0, &m_comm);
#endif // USE_MPI

	if (nCells > 0) {
		int rank = 0;
#ifdef USE_MPI
		MPI_Comm_rank(m_comm, &rank);
#endif // USE_MPI

		std::string outputName(static_cast<const char*>(info.buffer(OUTPUT_PREFIX)));
		outputName += "-surface";

    m_numVariables = 2*FREESURFACE_NUMBER_OF_COMPONENTS;
		std::vector<const char*> variables;
		for (unsigned int i = 0; i < m_numVariables; i++) {
			variables.push_back(LABELS[i]);
		}

		// TODO get the timestep from the checkpoint
		m_xdmfWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>(param.backend,
		                                                                              outputName.c_str(),
		                                                                              param.timestep);

#ifdef USE_MPI
		m_xdmfWriter->setComm(m_comm);
#endif // USE_MPI
		m_xdmfWriter->setBackupTimeStamp(param.backupTimeStamp);
		std::string extraIntVarName = "locationFlag";

		m_xdmfWriter->init(variables, std::vector<const char*>(), extraIntVarName.c_str());
		m_xdmfWriter->setMesh(nCells,
		                      static_cast<const unsigned int*>(info.buffer(CELLS)),
		                      nVertices,
		                      static_cast<const double*>(info.buffer(VERTICES)),
		                      param.timestep != 0);
		setLocationFlagData(static_cast<const unsigned int*>(info.buffer(LOCATIONFLAGS)));

		logInfo(rank) << "Initializing free surface output. Done.";
	}
}

char const * const seissol::writer::FreeSurfaceWriterExecutor::LABELS[] = {
	"v1", "v2", "v3", "u1", "u2", "u3"
};

