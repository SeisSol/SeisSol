/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#include <Kernels/Precision.h>
#include <async/ExecInfo.h>
#include <mpi.h>
#include <string>
#include <utils/env.h>
#include <vector>

#include "FreeSurfaceWriterExecutor.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "utils/logger.h"

/**
 * Initialize the XDMF writers
 */
void seissol::writer::FreeSurfaceWriterExecutor::execInit(
    const async::ExecInfo& info, const seissol::writer::FreeSurfaceInitParam& param) {
  if (m_xdmfWriter != nullptr) {
    logError() << "Free surface writer already initialized.";
  }

  const unsigned int nCells = info.bufferSize(Cells) / (3 * sizeof(int));
  const unsigned int nVertices = info.bufferSize(Vertices) / (3 * sizeof(double));

#ifdef USE_MPI
  MPI_Comm_split(seissol::MPI::mpi.comm(), (nCells > 0 ? 0 : MPI_UNDEFINED), 0, &m_comm);
#endif // USE_MPI

  if (nCells > 0) {
    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(m_comm, &rank);
#endif // USE_MPI

    std::string outputName(static_cast<const char*>(info.buffer(OutputPrefix)));
    outputName += "-surface";

    m_numVariables = 2 * FREESURFACE_NUMBER_OF_COMPONENTS;
    std::vector<const char*> variables;
    variables.reserve(m_numVariables);
    for (unsigned int i = 0; i < m_numVariables; i++) {
      variables.push_back(Labels[i]);
    }

    // TODO get the timestep from the checkpoint
    m_xdmfWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>(
        param.backend, outputName.c_str(), param.timestep);

#ifdef USE_MPI
    m_xdmfWriter->setComm(m_comm);
#endif // USE_MPI
    m_xdmfWriter->setBackupTimeStamp(param.backupTimeStamp);
    const std::string extraIntVarName = "locationFlag";
    const auto vertexFilter = utils::Env::get<bool>("SEISSOL_VERTEXFILTER", true);
    m_xdmfWriter->init(
        variables, std::vector<const char*>(), extraIntVarName.c_str(), vertexFilter);
    m_xdmfWriter->setMesh(nCells,
                          static_cast<const unsigned int*>(info.buffer(Cells)),
                          nVertices,
                          static_cast<const double*>(info.buffer(Vertices)),
                          param.timestep != 0);
    setLocationFlagData(static_cast<const unsigned int*>(info.buffer(LocationFlags)));

    logInfo() << "Initializing free surface output. Done.";
  }
}

const char* const seissol::writer::FreeSurfaceWriterExecutor::Labels[] = {
    "v1", "v2", "v3", "u1", "u2", "u3"};
