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

#include "Parallel/MPI.h"

#include <string>
#include <vector>

#include "utils/logger.h"

#include "FaultWriterExecutor.h"

/**
 * Initialize the XDMF writers
 */
void seissol::writer::FaultWriterExecutor::execInit(const async::ExecInfo &info, const seissol::writer::FaultInitParam &param)
{
	if (m_xdmfWriter)
		logError() << "Wave field writer already initialized";

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
		outputName += "-fault";

		std::vector<const char*> variables;
		for (unsigned int i = 0; i < FaultInitParam::OUTPUT_MASK_SIZE; i++) {
			if (param.outputMask[i])
				variables.push_back(LABELS[i]);
		}
		m_numVariables = variables.size();

		// TODO get the timestep from the checkpoint
		m_xdmfWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>(param.backend,
			outputName.c_str(), param.timestep);

#ifdef USE_MPI
		m_xdmfWriter->setComm(m_comm);
#endif // USE_MPI
		m_xdmfWriter->setBackupTimeStamp(param.backupTimeStamp);

		m_xdmfWriter->init(variables, std::vector<const char*>(), "fault-tag", true, true);
		m_xdmfWriter->setMesh(nCells, static_cast<const unsigned int*>(info.buffer(CELLS)),
			nVertices, static_cast<const double*>(info.buffer(VERTICES)),
			param.timestep != 0);
		setFaultTagsData(static_cast<const unsigned int*>(info.buffer(FAULTTAGS)));


		logInfo(rank) << "Initializing XDMF fault output. Done.";
	}
}

char const * const seissol::writer::FaultWriterExecutor::LABELS[] = {
	"SRs", "SRd", "T_s", "T_d", "P_n", "u_n", "Mud", "StV", "Ts0", "Td0", "Pn0", "Sls", "Sld", "Vr", "ASl","PSR", "RT", "DS", "P_f", "Tmp"
};

