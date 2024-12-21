// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_

#include <string_view>
#include <vector>

#include "Geometry/MeshReader.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/Lut.h"
#include "Kernels/Receiver.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"
#include <Eigen/Dense>

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
  seissol::SeisSol& seissolInstance;

  public:
  ReceiverWriter(seissol::SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

  void init(const std::string& fileNamePrefix,
            double endTime,
            const seissol::initializer::parameters::ReceiverOutputParameters& parameters);

  void addPoints(const seissol::geometry::MeshReader& mesh,
                 const seissol::initializer::Lut& ltsLut,
                 const seissol::initializer::LTS& lts,
                 const GlobalData* global);

  kernels::ReceiverCluster* receiverCluster(unsigned clusterId, LayerType layer);
  //
  // Hooks
  //
  void syncPoint(double /*currentTime*/) override;
  void simulationStart() override;
  void shutdown() override;

  private:
  [[nodiscard]] std::string fileName(unsigned pointId) const;
  void writeHeader(unsigned pointId, const Eigen::Vector3d& point);

  std::string m_receiverFileName;
  std::string m_fileNamePrefix;
  double m_samplingInterval;
  std::vector<std::shared_ptr<kernels::DerivedReceiverQuantity>> derivedQuantities;
  // Map needed because LayerType enum casts weirdly to int.
  std::unordered_map<LayerType, std::vector<kernels::ReceiverCluster>> m_receiverClusters;
  Stopwatch m_stopwatch;
};
} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_RECEIVERWRITER_H_
