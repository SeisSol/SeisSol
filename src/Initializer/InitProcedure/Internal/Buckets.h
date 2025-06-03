// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BUCKETS_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BUCKETS_H_

#include <Common/Real.h>
#include <Initializer/InitProcedure/Internal/MeshLayout.h>
#include <Memory/MemoryContainer.h>
#include <vector>

namespace seissol::initializer::internal {
struct CommunicationInfo {
  void* buffer;
  std::size_t offset;
  RealType precision;
  std::size_t count;
  int tag;
  int rank;
};

void bucketsAndCommunication(memory::MemoryContainer& container,
                             const std::vector<ClusterLayout>& meshLayout);
} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BUCKETS_H_
