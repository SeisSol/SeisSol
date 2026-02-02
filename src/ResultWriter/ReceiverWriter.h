// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_

#include "Geometry/MeshReader.h"
#include "Kernels/Receiver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"

#include <Eigen/Dense>
#include <string_view>
#include <vector>

namespace seissol {
struct LocalIntegrationData;
struct GlobalData;
class SeisSol;
namespace initializer::parameters {
struct ReceiverOutputParameters;
} // namespace initializer::parameters
} // namespace seissol

namespace seissol::writer {
Eigen::Vector3d parseReceiverLine(const std::string& line);
std::vector<Eigen::Vector3d> parseReceiverFile(const std::string& receiverFileName);

class ReceiverWriter : public seissol::Module {
  private:
  seissol::SeisSol& seissolInstance_;

  public:
  explicit ReceiverWriter(seissol::SeisSol& seissolInstance) : seissolInstance_(seissolInstance) {}

  void init(const std::string& fileNamePrefix,
            double endTime,
            const seissol::initializer::parameters::ReceiverOutputParameters& parameters);

  void addPoints(const seissol::geometry::MeshReader& mesh,
                 const LTS::Backmap& backmap,
                 const CompoundGlobalData& global);

  kernels::ReceiverCluster* receiverCluster(std::size_t id);
  //
  // Hooks
  //
  void syncPoint(double /*currentTime*/) override;
  void simulationStart(std::optional<double> checkpointTime) override;
  void shutdown() override;

  private:
  [[nodiscard]] std::string fileName(unsigned pointId) const;
  void writeHeader(unsigned pointId, const Eigen::Vector3d& point);

  std::string receiverFileName_;
  std::string fileNamePrefix_;
  double samplingInterval_{};
  std::vector<std::shared_ptr<kernels::DerivedReceiverQuantity>> derivedQuantities_;
  std::vector<std::shared_ptr<kernels::ReceiverCluster>> receiverClusters_;
  Stopwatch stopwatch_;
};
} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
