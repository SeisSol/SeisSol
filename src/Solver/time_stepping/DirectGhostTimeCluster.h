// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIME_STEPPING_DIRECTGHOSTTIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIME_STEPPING_DIRECTGHOSTTIMECLUSTER_H_

#include <list>
#include "Initializer/Typedefs.h"
#include "Solver/time_stepping/AbstractGhostTimeCluster.h"


namespace seissol::time_stepping {
class DirectGhostTimeCluster : public AbstractGhostTimeCluster {
protected:
  virtual void sendCopyLayer();
  virtual void receiveGhostLayer();
  virtual bool testForGhostLayerReceives();

public:
    DirectGhostTimeCluster(double maxTimeStepSize,
                           int timeStepRate,
                           int globalTimeClusterId,
                           int otherGlobalTimeClusterId,
                           const MeshStructure* meshStructure,
                           bool persistent);
    void finalize() override;
private:
  bool persistent;
};
} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIME_STEPPING_DIRECTGHOSTTIMECLUSTER_H_

