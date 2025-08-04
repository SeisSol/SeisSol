// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_
#define SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_

#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include "Geometry/MeshReader.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"

#include "Initializer/Parameters/SeisSolParameters.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"
#include <Solver/MultipleSimulations.h>

namespace seissol {
class SeisSol;
namespace writer {

struct EnergiesStorage {
  static constexpr size_t NumberOfEnergies = 13;
  std::array<double, multisim::NumSimulations * NumberOfEnergies> energies{};

  double& gravitationalEnergy(size_t sim);

  double& acousticEnergy(size_t sim);

  double& acousticKineticEnergy(size_t sim);

  double& elasticEnergy(size_t sim);

  double& elasticKineticEnergy(size_t sim);

  double& totalFrictionalWork(size_t sim);

  double& staticFrictionalWork(size_t sim);

  double& plasticMoment(size_t sim);

  double& seismicMoment(size_t sim);

  double& potency(size_t sim);

  double& totalMomentumX(size_t sim);
  double& totalMomentumY(size_t sim);
  double& totalMomentumZ(size_t sim);
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

  explicit EnergyOutput(seissol::SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

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

  seissol::SeisSol& seissolInstance;

  bool shouldComputeVolumeEnergies() const;

  bool isEnabled = false;
  bool isTerminalOutputEnabled = false;
  bool isFileOutputEnabled = false;
  bool isPlasticityEnabled = false;
  bool isCheckAbortCriteraSlipRateEnabled = false;
  bool isCheckAbortCriteraMomentRateEnabled = false;
  int computeVolumeEnergiesEveryOutput = 1;
  int outputId = 0;

  std::string outputFileName;
  std::ofstream out;

#ifdef ACL_DEVICE
  void* timeDerivativePlusHost = nullptr;
  void* timeDerivativeMinusHost = nullptr;
  void* timeDerivativePlusHostMapped = nullptr;
  void* timeDerivativeMinusHostMapped = nullptr;
#endif

  const GlobalData* global = nullptr;
  const DynamicRupture::Storage* drStorage = nullptr;
  const seissol::geometry::MeshReader* meshReader = nullptr;
  const LTS::Storage* ltsStorage = nullptr;

  EnergiesStorage energiesStorage{};
  std::array<double, multisim::NumSimulations> minTimeSinceSlipRateBelowThreshold;
  std::array<double, multisim::NumSimulations> minTimeSinceMomentRateBelowThreshold;
  double terminatorMaxTimePostRupture{};
  double energyOutputInterval{};
  double terminatorMomentRateThreshold{};
  std::array<double, multisim::NumSimulations> seismicMomentPrevious;

  std::size_t approxElements;
};

} // namespace writer
} // namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_
