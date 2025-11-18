// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "FreeSurfaceWriterExecutor.h"

#include "Kernels/Precision.h"
#include "Parallel/MPI.h"
#include "Solver/FreeSurfaceIntegrator.h"

#include <async/ExecInfo.h>
#include <mpi.h>
#include <string>
#include <utils/env.h>
#include <utils/logger.h>
#include <vector>

/**
 * Initialize the XDMF writers
 */
void seissol::writer::FreeSurfaceWriterExecutor::execInit(
    const async::ExecInfo& info, const seissol::writer::FreeSurfaceInitParam& param) {
  if (m_xdmfWriter != nullptr) {
    logError() << "Free surface writer already initialized.";
  }

  logInfo() << "X1";

  const unsigned int nCells = info.bufferSize(Cells) / (3 * sizeof(int));
  const unsigned int nVertices = info.bufferSize(Vertices) / (3 * sizeof(double));

  MPI_Comm_split(seissol::Mpi::mpi.comm(), (nCells > 0 ? 0 : MPI_UNDEFINED), 0, &m_comm);

  logInfo() << "X2";
  m_enabled = true;

  if (nCells > 0) {
    int rank = 0;
    MPI_Comm_rank(m_comm, &rank);
    logInfo() << "X3";

    std::string outputName(static_cast<const char*>(info.buffer(OutputPrefix)));
    outputName += "-surface";

    m_numVariables = 2 * seissol::solver::FreeSurfaceIntegrator::NumComponents;
    std::vector<const char*> variables;
    variables.reserve(m_numVariables);
    for (unsigned int i = 0; i < m_numVariables; i++) {
      variables.push_back(Labels[i]);
    }

    logInfo() << "X4";

    // TODO get the timestep from the checkpoint
    m_xdmfWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>(
        param.backend, outputName.c_str(), param.timestep);

    logInfo() << "X5";

    m_xdmfWriter->setComm(m_comm);
    logInfo() << "X6";
    m_xdmfWriter->setBackupTimeStamp(param.backupTimeStamp);
    logInfo() << "X7";
    const std::string extraIntVarName = "locationFlag";
    const auto vertexFilter = utils::Env("").get<bool>("SEISSOL_VERTEXFILTER", true);
    m_xdmfWriter->init(
        variables, std::vector<const char*>(), {extraIntVarName, "global-id"}, vertexFilter);
    logInfo() << "X8";
    m_xdmfWriter->setMesh(nCells,
                          static_cast<const unsigned int*>(info.buffer(Cells)),
                          nVertices,
                          static_cast<const double*>(info.buffer(Vertices)),
                          param.timestep != 0);
    logInfo() << "X9";
    setLocationFlagData(static_cast<const unsigned int*>(info.buffer(LocationFlags)));
    logInfo() << "X10";
    m_xdmfWriter->writeExtraIntCellData(1,
                                        static_cast<const unsigned int*>(info.buffer(GlobalIds)));
    logInfo() << "X11";

    logInfo() << "Initializing free surface output. Done.";
  }
  logInfo() << "X12";
}

const char* const seissol::writer::FreeSurfaceWriterExecutor::Labels[] = {
    "v1", "v2", "v3", "u1", "u2", "u3"};
