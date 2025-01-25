// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "Parallel/MPI.h"

#include <Kernels/Precision.h>
#include <async/ExecInfo.h>
#include <mpi.h>
#include <string>
#include <utils/env.h>
#include <vector>

#include "utils/logger.h"

#include "FaultWriterExecutor.h"

/**
 * Initialize the XDMF writers
 */
void seissol::writer::FaultWriterExecutor::execInit(const async::ExecInfo& info,
                                                    const seissol::writer::FaultInitParam& param) {
  if (m_xdmfWriter != nullptr) {
    logError() << "Wave field writer already initialized";
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
    outputName += "-fault";

    std::vector<const char*> variables;
    for (unsigned int i = 0; i < FaultInitParam::OutputMaskSize; i++) {
      if (param.outputMask[i]) {
        variables.push_back(Labels[i]);
      }
    }
    m_numVariables = variables.size();

    // TODO get the timestep from the checkpoint
    m_xdmfWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>(
        param.backend, outputName.c_str(), param.timestep);

#ifdef USE_MPI
    m_xdmfWriter->setComm(m_comm);
#endif // USE_MPI
    m_xdmfWriter->setBackupTimeStamp(param.backupTimeStamp);
    const auto vertexFilter = utils::Env::get<bool>("SEISSOL_VERTEXFILTER", true);
    m_xdmfWriter->init(variables, std::vector<const char*>(), "fault-tag", vertexFilter, true);
    m_xdmfWriter->setMesh(nCells,
                          static_cast<const unsigned int*>(info.buffer(Cells)),
                          nVertices,
                          static_cast<const double*>(info.buffer(Vertices)),
                          param.timestep != 0);
    setFaultTagsData(static_cast<const unsigned int*>(info.buffer(FaultTags)));

    logInfo() << "Initializing XDMF fault output. Done.";
  }
}

const char* const seissol::writer::FaultWriterExecutor::Labels[] = {
    "SRs", "SRd", "T_s", "T_d", "P_n", "u_n", "Mud", "StV", "Ts0", "Td0",
    "Pn0", "Sls", "Sld", "Vr",  "ASl", "PSR", "RT",  "DS",  "P_f", "Tmp"};
