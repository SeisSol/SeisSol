// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_HELPER_H_
#define SEISSOL_SRC_PARALLEL_HELPER_H_

#include "Common/Marker.h"
#include "Parallel/MPI.h"

#include <utils/env.h>
#include <utils/logger.h>

namespace seissol {

// TODO: refactor so that we always "carry" the env object; and not need to re-create it

void printCommThreadInfo(const Mpi& mpiBasic, utils::Env& env);

bool useCommThread(const Mpi& mpiBasic, utils::Env& env);

bool usePersistentMpi(utils::Env& env);

void printPersistentMpiInfo(utils::Env& env);

bool useUSM(SEISSOL_GPU_PARAM utils::Env& env);

bool useUSM();

void printUSMInfo(utils::Env& env);

bool useMPIUSM(SEISSOL_GPU_PARAM utils::Env& env);

bool useMPIUSM();

void printMPIUSMInfo(utils::Env& env);

bool useDeviceL2Compress(utils::Env& env);

bool useDeviceL2Compress();

void printDeviceL2Compress(utils::Env& env);

enum class DataTransferMode { Direct, CopyInCopyOutHost };

DataTransferMode getDataTransferMode(utils::Env& env);

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_HELPER_H_
