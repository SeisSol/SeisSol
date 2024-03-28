#ifndef ENERGYOUTPUT_H
#define ENERGYOUTPUT_H

#include <array>
#include <string>
#include <fstream>
#include <iostream>

#include <Initializer/typedefs.hpp>
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/Lut.hpp>

#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "Initializer/Parameters/SeisSolParameters.h"

namespace seissol {
class SeisSol;
namespace writer {

struct EnergiesStorage {
  std::array<double, 10> energies{};

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

  private:
  real computeStaticWork(const real* degreesOfFreedomPlus,
                         const real* degreesOfFreedomMinus,
                         DRFaceInformation const& faceInfo,
                         DRGodunovData const& godunovData,
                         const real slip[seissol::tensor::slipInterpolated::size()]);

  void computeDynamicRuptureEnergies();

  void computeVolumeEnergies();

  void computeEnergies();

  void reduceEnergies();

  void reduceMinTimeSinceSlipRateBelowThreshold();

  void printEnergies();

  void checkAbortCriterion(real timeSinceThreshold, const std::string& prefix_message);

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

  const GlobalData* global = nullptr;
  seissol::initializer::DynamicRupture* dynRup = nullptr;
  seissol::initializer::LTSTree* dynRupTree = nullptr;
  seissol::geometry::MeshReader* meshReader = nullptr;
  seissol::initializer::LTSTree* ltsTree = nullptr;
  seissol::initializer::LTS* lts = nullptr;
  seissol::initializer::Lut* ltsLut = nullptr;

  EnergiesStorage energiesStorage{};
  real minTimeSinceSlipRateBelowThreshold;
  real minTimeSinceMomentRateBelowThreshold = 0.0;
  double terminatorMaxTimePostRupture;
  double energyOutputInterval;
  double terminatorMomentRateThreshold;
  double seismicMomentPrevious = 0.0;
};

} // namespace writer
} // namespace seissol

#endif // ENERGYOUTPUT_H
