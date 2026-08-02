// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_
#define SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_

#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/MemoryAllocator.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <fstream>
#include <iostream>
#include <string>

namespace seissol {
class SeisSol;
namespace writer {

class EnergiesStorage {
  public:
  void setSimcount(size_t count);

  size_t addEnergy(const std::string& name);

  double& energy(size_t handle, size_t sim);

  [[nodiscard]] double energy(size_t handle, size_t sim) const;

  double& energy(const std::string& name, size_t sim);

  [[nodiscard]] double energy(const std::string& name, size_t sim) const;

  [[nodiscard]] const std::map<std::string, size_t>& handles() const;

  std::vector<double>& values();

  void reset();

  private:
  std::size_t simcount_{1};
  std::vector<double> values_;
  std::map<std::string, size_t> handles_;
};

class EnergyOutput : public Module {
  public:
  void init(GlobalData* newGlobal,
            const DynamicRupture::Storage& newDynRuptTree,
            const seissol::geometry::MeshReader& newMeshReader,
            const LTS::Storage& newStorage,
            bool newIsPlasticityEnabled,
            const std::string& outputFileNamePrefix,
            const seissol::initializer::parameters::EnergyOutputParameters& parameters);

  void syncPoint(double time) override;

  void simulationStart(std::optional<double> checkpointTime) override;

  explicit EnergyOutput(seissol::SeisSol& seissolInstance) : seissolInstance_(seissolInstance) {}

  ~EnergyOutput() override;

  auto operator=(const EnergyOutput&) = delete;
  auto operator=(EnergyOutput&&) = delete;
  EnergyOutput(const EnergyOutput&) = delete;
  EnergyOutput(EnergyOutput&&) = delete;

  private:
  void computeDynamicRuptureEnergies();

  void computeVolumeEnergies();

  void computeEnergies();

  void reduceEnergies();

  void reduceMinTimeSinceSlipRateBelowThreshold();

  void printEnergies();

  void checkAbortCriterion(const std::array<double, multisim::NumSimulations>& timeSinceThreshold,
                           const std::string& prefixMessage);

  void writeHeader();

  void writeEnergies(double time);

  seissol::SeisSol& seissolInstance_;

  bool shouldComputeVolumeEnergies() const;

  bool isEnabled_ = false;
  bool isTerminalOutputEnabled_ = false;
  bool isFileOutputEnabled_ = false;
  bool isPlasticityEnabled_ = false;
  bool isCheckAbortCriteraSlipRateEnabled_ = false;
  bool isCheckAbortCriteraMomentRateEnabled_ = false;
  int computeVolumeEnergiesEveryOutput_ = 1;
  int outputId_ = 0;

  std::string outputFileName_;
  std::ofstream out_;

  const GlobalData* global_ = nullptr;
  const DynamicRupture::Storage* drStorage_ = nullptr;
  const seissol::geometry::MeshReader* meshReader_ = nullptr;
  const LTS::Storage* ltsStorage_ = nullptr;

  EnergiesStorage energiesStorage_{};
  std::array<double, multisim::NumSimulations> minTimeSinceSlipRateBelowThreshold_{};
  std::array<double, multisim::NumSimulations> minTimeSinceMomentRateBelowThreshold_{};
  double terminatorMaxTimePostRupture_{};
  double energyOutputInterval_{};
  double terminatorMomentRateThreshold_{};
  std::array<double, multisim::NumSimulations> seismicMomentPrevious_{};
};

} // namespace writer
} // namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_
