#ifndef ENERGYOUTPUT_H
#define ENERGYOUTPUT_H

#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "Geometry/MeshReader.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/LTS.h"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/typedefs.hpp"
#include <Geometry/MeshReader.h>
#include <Initializer/DynamicRupture.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>
#include <Initializer/typedefs.hpp>
#include <Solver/MultipleSimulations.h>

#include "Initializer/Parameters/SeisSolParameters.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"

namespace seissol {
class SeisSol;
namespace writer {

struct EnergiesStorage {
  static constexpr size_t NumberOfEnergies = 10;
  std::array<double, multipleSimulations::numberOfSimulations * NumberOfEnergies> energies{};

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
            std::array<std::shared_ptr<seissol::initializer::DynamicRupture>, MULTIPLE_SIMULATIONS>& newDynRup,
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

  private:
  std::array<real, multipleSimulations::numberOfSimulations>
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
      const real (&timeSinceThreshold)[multipleSimulations::numberOfSimulations],
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
  real minTimeSinceSlipRateBelowThreshold[multipleSimulations::numberOfSimulations] = {0.0};
  real minTimeSinceMomentRateBelowThreshold[multipleSimulations::numberOfSimulations] = {0.0};
  double terminatorMaxTimePostRupture;
  double energyOutputInterval;
  double terminatorMomentRateThreshold;
  double seismicMomentPrevious[multipleSimulations::numberOfSimulations] = {0.0};
};

} // namespace writer
} // namespace seissol

#endif // ENERGYOUTPUT_H
