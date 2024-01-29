/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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

#include "MPI.h"
#include "utils/stringutils.h"
#include <unistd.h>
#include <cstdlib>
#include <cctype>

#ifdef ACL_DEVICE
#include "Parallel/AcceleratorDevice.h"
#endif

void seissol::MPI::init(int& argc, char**& argv) {
  // Note: Strictly speaking, we only require MPI_THREAD_MULTIPLE if using
  // a communication thread and/or async I/O.
  // The safer (and more sane) option is to enable it by default.
  int required = MPI_THREAD_MULTIPLE;
  int provided;
  MPI_Init_thread(&argc, &argv, required, &provided);

  setComm(MPI_COMM_WORLD);

  std::string hostName(256, ' ');
  if (gethostname(const_cast<char*>(hostName.c_str()), 256) != 0) {
    hostName = "unknown-host";
  } else {
    utils::StringUtils::rtrim(hostName);
    hostName.pop_back();
  }
  hostNames = collectContainer(hostName);

  // Test this after setComm() to get the correct m_rank
  if (provided < required) {
    logError() << utils::nospace << "Provided MPI thread support (" << provided
               << ") is smaller than required thread support (" << required << ").";
  }
}

void seissol::MPI::setComm(MPI_Comm comm) {
  m_comm = comm;

  MPI_Comm_rank(comm, &m_rank);
  MPI_Comm_size(comm, &m_size);

  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &m_sharedMemComm);
  MPI_Comm_rank(m_sharedMemComm, &m_sharedMemMpiRank);
  MPI_Comm_size(m_sharedMemComm, &m_sharedMemMpiSize);
}

#ifdef ACL_DEVICE
void seissol::MPI::bindAcceleratorDevice() {
  auto& instance = seissol::AcceleratorDevice::getInstance();
  instance.bindAcceleratorDevice(0);
}
#endif

void seissol::MPI::setDataTransferModeFromEnv() {
  // TODO (Ravil, David): switch to reading this option from the parameter-file.
  // Waiting for David to finish his `no-fortran` PR
  if (const char* envVariable = std::getenv("SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE")) {
    std::string option{envVariable};
    std::transform(option.begin(), option.end(), option.begin(), [](unsigned char c) {
      return std::tolower(c);
    });

    if (option == "direct") {
      preferredDataTransferMode = DataTransferMode::Direct;
    } else if (option == "host") {
      preferredDataTransferMode = DataTransferMode::CopyInCopyOutHost;
    } else {
      logWarning(m_rank) << "Ignoring `SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE`."
                         << "Expected values: direct, host.";
      option = "direct";
    }
#ifndef ACL_DEVICE
    if (preferredDataTransferMode != DataTransferMode::Direct) {
      logWarning(m_rank) << "The CPU version of SeisSol supports"
                         << "only the `direct` MPI transfer mode.";
      option = "direct";
      preferredDataTransferMode = DataTransferMode::Direct;
    }
#endif
      logInfo(m_rank) << "Selected" << option << "MPI data transfer mode as the preferred one";
  }
}

seissol::MPI seissol::MPI::mpi;
