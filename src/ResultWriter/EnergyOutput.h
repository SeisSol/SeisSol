// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
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
#include <memory>
#include <string>

#include "Geometry/MeshReader.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Lut.h"
#include "Initializer/Typedefs.h"

#include "Initializer/Parameters/SeisSolParameters.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"

namespace seissol {
class SeisSol;
namespace writer {

struct EnergiesStorage {
  static constexpr size_t NumberOfEnergies = 10;
  std::array<double, MULTIPLE_SIMULATIONS * NumberOfEnergies> energies{};

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
            std::array<std::shared_ptr<seissol::initializer::DynamicRupture>, MULTIPLE_SIMULATIONS>&
                newDynRup,
            std::array<seissol::initializer::LTSTree*, MULTIPLE_SIMULATIONS>& newDynRuptTree,
            seissol::geometry::MeshReader* newMeshReader,
            seissol::initializer::LTSTree* newLtsTree,
            seissol::initializer::LTS* newLts,
            seissol::initializer::Lut* newLtsLut,
            bool newIsPlasticityEnabled,
            const std::string& outputFileNamePrefix,
            const seissol::initializer::parameters::EnergyOutputParameters& parameters);

  void syncPoint(double time) override;

  void simulationStart() override;

  EnergyOutput(seissol::SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

  ~EnergyOutput() override;

  auto operator=(const EnergyOutput&) = delete;
  auto operator=(EnergyOutput&&) = delete;
  EnergyOutput(const EnergyOutput&) = delete;
  EnergyOutput(EnergyOutput&&) = delete;

  private:
  std::array<real, MULTIPLE_SIMULATIONS>
      computeStaticWork(const real* degreesOfFreedomPlus,
                        const real* degreesOfFreedomMinus,
                        const DRFaceInformation& faceInfo,
                        const DRGodunovData& godunovData,
                        const real slip[seissol::tensor::slipInterpolated::size()]);

  void computeDynamicRuptureEnergies();

  void computeVolumeEnergies();

  void computeEnergies();

  void reduceEnergies();

  void reduceMinTimeSinceSlipRateBelowThreshold();

  void printEnergies();

  void checkAbortCriterion(
      const real (&timeSinceThreshold)[MULTIPLE_SIMULATIONS],
      const std::string& prefix_message);

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
  real* timeDerivativePlusHost = nullptr;
  real* timeDerivativeMinusHost = nullptr;
  real* timeDerivativePlusHostMapped = nullptr;
  real* timeDerivativeMinusHostMapped = nullptr;
#endif

  const GlobalData* global = nullptr;
  // seissol::initializer::DynamicRupture* dynRup = nullptr;
  // seissol::initializer::LTSTree* dynRupTree = nullptr;
  std::array<std::shared_ptr<seissol::initializer::DynamicRupture>, MULTIPLE_SIMULATIONS> dynRup;
  std::array<seissol::initializer::LTSTree*, MULTIPLE_SIMULATIONS> dynRupTree;
  seissol::geometry::MeshReader* meshReader = nullptr;
  seissol::initializer::LTSTree* ltsTree = nullptr;
  seissol::initializer::LTS* lts = nullptr;
  seissol::initializer::Lut* ltsLut = nullptr;

  EnergiesStorage energiesStorage{};
  real minTimeSinceSlipRateBelowThreshold[MULTIPLE_SIMULATIONS] = {0.0};
  real minTimeSinceMomentRateBelowThreshold[MULTIPLE_SIMULATIONS] = {0.0};
  double terminatorMaxTimePostRupture;
  double energyOutputInterval;
  double terminatorMomentRateThreshold;
  double seismicMomentPrevious[MULTIPLE_SIMULATIONS] = {0.0};
};

} // namespace writer
} // namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_
