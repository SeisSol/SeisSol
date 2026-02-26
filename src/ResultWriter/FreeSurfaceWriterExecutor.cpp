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
  if (xdmfWriter_ != nullptr) {
    logError() << "Free surface writer already initialized.";
  }

  const unsigned int nCells = info.bufferSize(Cells) / (3 * sizeof(int));
  const unsigned int nVertices = info.bufferSize(Vertices) / (3 * sizeof(double));

  MPI_Comm_split(seissol::Mpi::mpi.comm(), (nCells > 0 ? 0 : MPI_UNDEFINED), 0, &comm_);

  enabled_ = true;

  if (nCells > 0) {
    int rank = 0;
    MPI_Comm_rank(comm_, &rank);

    std::string outputName(static_cast<const char*>(info.buffer(OutputPrefix)));
    outputName += "-surface";

    numVariables_ = 2 * seissol::solver::FreeSurfaceIntegrator::NumComponents;
    std::vector<const char*> variables;
    variables.reserve(numVariables_);
    for (unsigned int i = 0; i < numVariables_; i++) {
      variables.push_back(Labels[i]);
    }

    // TODO get the timestep from the checkpoint
    xdmfWriter_ = new xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>(
        param.backend, outputName.c_str(), param.timestep);

    xdmfWriter_->setComm(comm_);
    xdmfWriter_->setBackupTimeStamp(param.backupTimeStamp);
    const std::string extraIntVarName = "locationFlag";
    const auto vertexFilter = utils::Env("").get<bool>("SEISSOL_VERTEXFILTER", true);
    xdmfWriter_->init(
        variables, std::vector<const char*>(), {extraIntVarName, "global-id"}, vertexFilter);
    xdmfWriter_->setMesh(nCells,
                         static_cast<const unsigned int*>(info.buffer(Cells)),
                         nVertices,
                         static_cast<const double*>(info.buffer(Vertices)),
                         param.timestep != 0);
    setLocationFlagData(static_cast<const unsigned int*>(info.buffer(LocationFlags)));
    xdmfWriter_->writeExtraIntCellData(1, static_cast<const unsigned int*>(info.buffer(GlobalIds)));

    logInfo() << "Initializing free surface output. Done.";
  }
}

const char* const seissol::writer::FreeSurfaceWriterExecutor::Labels[] = {
    "v1", "v2", "v3", "u1", "u2", "u3"};
