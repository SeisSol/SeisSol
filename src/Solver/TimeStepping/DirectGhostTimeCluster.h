// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_DIRECTGHOSTTIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_DIRECTGHOSTTIMECLUSTER_H_

#include "Initializer/Typedefs.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <list>

namespace seissol::time_stepping {
class DirectGhostTimeCluster : public AbstractGhostTimeCluster {
  protected:
  void sendCopyLayer() override;
  void receiveGhostLayer() override;
  bool testForGhostLayerReceives() override;

  public:
  DirectGhostTimeCluster(double maxTimeStepSize,
                         std::uint64_t timeStepRate,
                         std::size_t globalTimeClusterId,
                         std::size_t otherGlobalTimeClusterId,
                         const std::string& displayName,
                         const std::string& otherDisplayName,
                         const seissol::solver::HaloCommunication& meshStructure,
                         bool persistent);
  void finalize() override;

  private:
  bool persistent_;
};
} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_DIRECTGHOSTTIMECLUSTER_H_
