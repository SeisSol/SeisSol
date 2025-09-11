// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

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

  MPI_Comm_split(seissol::MPI::mpi.comm(), (nCells > 0 ? 0 : MPI_UNDEFINED), 0, &m_comm);

  m_enabled = true;

  if (nCells > 0) {
    int rank = 0;
    MPI_Comm_rank(m_comm, &rank);

    std::string outputName(static_cast<const char*>(info.buffer(OutputPrefix)));
    outputName += "-surface";

    m_numVariables = 2 * seissol::solver::FreeSurfaceIntegrator::NumComponents;
    std::vector<const char*> variables;
    variables.reserve(m_numVariables);
    for (unsigned int i = 0; i < m_numVariables; i++) {
      variables.push_back(Labels[i]);
    }

    // TODO get the timestep from the checkpoint
    m_xdmfWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>(
        param.backend, outputName.c_str(), param.timestep);

    m_xdmfWriter->setComm(m_comm);
    m_xdmfWriter->setBackupTimeStamp(param.backupTimeStamp);
    const std::string extraIntVarName = "locationFlag";
    const auto vertexFilter = utils::Env("").get<bool>("SEISSOL_VERTEXFILTER", true);
    m_xdmfWriter->init(
        variables, std::vector<const char*>(), {extraIntVarName, "global-id"}, vertexFilter);
    m_xdmfWriter->setMesh(nCells,
                          static_cast<const unsigned int*>(info.buffer(Cells)),
                          nVertices,
                          static_cast<const double*>(info.buffer(Vertices)),
                          param.timestep != 0);
    setLocationFlagData(static_cast<const unsigned int*>(info.buffer(LocationFlags)));
    m_xdmfWriter->writeExtraIntCellData(1,
                                        static_cast<const unsigned int*>(info.buffer(GlobalIds)));

    logInfo() << "Initializing free surface output. Done.";
  }
}

const char* const seissol::writer::FreeSurfaceWriterExecutor::Labels[] = {
    "v1", "v2", "v3", "u1", "u2", "u3"};
