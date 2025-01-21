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
#include <string>

#include "Geometry/MeshReader.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Lut.h"

#include "Initializer/Parameters/SeisSolParameters.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"

namespace seissol {
class SeisSol;
namespace writer {

struct EnergiesStorage {
  std::array<double, 13> energies{};

  double& gravitationalEnergy();

  double& acousticEnergy();

  double& acousticKineticEnergy();

  double& elasticEnergy();

  double& elasticKineticEnergy();

  double& totalFrictionalWork();

  double& staticFrictionalWork();

  double& plasticMoment();

  double& seismicMoment();

  double& potency();

  double& totalMomentumX();
  double& totalMomentumY();
  double& totalMomentumZ();
};

class EnergyOutput : public Module {
  public:
  void init(GlobalData* newGlobal,
            seissol::initializer::DynamicRupture* newDynRup,
            seissol::initializer::LTSTree* newDynRuptTree,
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
  real computeStaticWork(const real* degreesOfFreedomPlus,
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

  void checkAbortCriterion(real timeSinceThreshold, const std::string& prefixMessage);

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
  seissol::initializer::DynamicRupture* dynRup = nullptr;
  seissol::initializer::LTSTree* dynRupTree = nullptr;
  seissol::geometry::MeshReader* meshReader = nullptr;
  seissol::initializer::LTSTree* ltsTree = nullptr;
  seissol::initializer::LTS* lts = nullptr;
  seissol::initializer::Lut* ltsLut = nullptr;

  EnergiesStorage energiesStorage{};
  real minTimeSinceSlipRateBelowThreshold{};
  real minTimeSinceMomentRateBelowThreshold = 0.0;
  double terminatorMaxTimePostRupture{};
  double energyOutputInterval{};
  double terminatorMomentRateThreshold{};
  double seismicMomentPrevious = 0.0;
};

} // namespace writer
} // namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_ENERGYOUTPUT_H_
