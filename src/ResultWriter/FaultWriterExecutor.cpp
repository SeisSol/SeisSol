/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

		m_xdmfWriter->init(variables, std::vector<const char*>());
		m_xdmfWriter->setMesh(nCells, static_cast<const unsigned int*>(info.buffer(CELLS)),
			nVertices, static_cast<const double*>(info.buffer(VERTICES)),
			param.timestep != 0);

		logInfo(rank) << "Initializing XDMF fault output. Done.";
	}
}

char const * const seissol::writer::FaultWriterExecutor::LABELS[] = {
	"SRs", "SRd", "T_s", "T_d", "P_n", "u_n", "Mud", "StV", "Ts0", "Td0", "Pn0", "Sls", "Sld", "Vr", "ASl","PSR", "RT", "DS", "P_f", "Tmp"
};
