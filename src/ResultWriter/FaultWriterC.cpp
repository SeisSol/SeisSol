/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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

#include "xdmfwriter/XdmfWriter.h"
#include "FaultWriterC.h"

static xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE>* xdmfWriter = 0L;

static char const * const labels[] = {
  "SRs", "SRd", "T_s", "T_d", "P_n", "u_n", "Mud", "StV", "Ts0", "Td0", "Pn0", "Sls", "Sld", "Vr", "ASl","PSR", "RT", "DS"
};

#ifdef USE_MPI
static MPI_Comm comm = MPI_COMM_NULL;
#endif // USE_MPI

extern "C"
{

void fault_create_comm(int dr)
{
#ifdef USE_MPI
	MPI_Comm_split(seissol::MPI::mpi.comm(), (dr ? 0 : MPI_UNDEFINED), 0, &comm);
#endif // USE_MPI
}

void fault_hdf_init(const int* cells, const double* vertices,
		int nCells, int nVertices,
		int* outputMask, const char* outputPrefix)
{
	if (nCells > 0) {
		int rank = 0;
#ifdef USE_MPI
		MPI_Comm_rank(comm, &rank);
#endif // USE_MPI

		logInfo(rank) << "Initializing fault output.";

		std::string outputName(outputPrefix);
		outputName += "-fault";

		std::vector<const char*> variables;
		if (outputMask[0]) {
			variables.push_back(labels[0]);
			variables.push_back(labels[1]);
		}
		if (outputMask[1]) {
			variables.push_back(labels[2]);
			variables.push_back(labels[3]);
			variables.push_back(labels[4]);
		}
		if (outputMask[2])
			variables.push_back(labels[5]);
		if (outputMask[3]) {
			variables.push_back(labels[6]);
			variables.push_back(labels[7]);
		}
		if (outputMask[4]) {
			variables.push_back(labels[8]);
			variables.push_back(labels[9]);
			variables.push_back(labels[10]);
		}
		if (outputMask[5]) {
			variables.push_back(labels[11]);
			variables.push_back(labels[12]);
		}
		if (outputMask[6])
			variables.push_back(labels[13]);
		if (outputMask[7])
			variables.push_back(labels[14]);
		if (outputMask[8])
			variables.push_back(labels[15]);
		if (outputMask[9])
			variables.push_back(labels[16]);
		if (outputMask[10])
			variables.push_back(labels[17]);

		// TODO get the timestep from the checkpoint
		xdmfWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE>(rank,
				outputName.c_str(), variables, 0);
#ifdef USE_MPI
		xdmfWriter->setComm(comm);
#endif // USE_MPI
		xdmfWriter->init(nCells, reinterpret_cast<const unsigned int*>(cells), nVertices, vertices, true);
	}
}

void fault_hdf_close()
{
	if (xdmfWriter) {
		xdmfWriter->close();
		delete xdmfWriter;
	}
}

void fault_hdf_write(double time)
{
	if (xdmfWriter) {
		xdmfWriter->addTimeStep(time);
	}
}

void fault_hdf_write_data(int id, const double* data)
{
	if (xdmfWriter) {
		xdmfWriter->writeData(id, data);
	}
}

void fault_hdf_flush()
{
	if (xdmfWriter) {
		xdmfWriter->flush();
	}
}

}
