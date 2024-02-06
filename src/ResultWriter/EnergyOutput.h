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
#include "Initializer/InputParameters.hpp"

namespace seissol::writer {

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
            seissol::initializers::DynamicRupture* newDynRup,
            seissol::initializers::LTSTree* newDynRuptTree,
            seissol::geometry::MeshReader* newMeshReader,
            seissol::initializers::LTSTree* newLtsTree,
            seissol::initializers::LTS* newLts,
            seissol::initializers::Lut* newLtsLut,
            bool newIsPlasticityEnabled,
            const std::string& outputFileNamePrefix,
            const seissol::initializer::parameters::EnergyOutputParameters& parameters);

  void syncPoint(double time) override;

  void simulationStart() override;

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

  void printEnergies();

  void writeHeader();

  void writeEnergies(double time);

  bool shouldComputeVolumeEnergies() const;

  bool isEnabled = false;
  bool isTerminalOutputEnabled = false;
  bool isFileOutputEnabled = false;
  bool isPlasticityEnabled = false;
  int computeVolumeEnergiesEveryOutput = 1;
  int outputId = 0;

  std::string outputFileName;
  std::ofstream out;

  const GlobalData* global = nullptr;
  seissol::initializers::DynamicRupture* dynRup = nullptr;
  seissol::initializers::LTSTree* dynRupTree = nullptr;
  seissol::geometry::MeshReader* meshReader = nullptr;
  seissol::initializers::LTSTree* ltsTree = nullptr;
  seissol::initializers::LTS* lts = nullptr;
  seissol::initializers::Lut* ltsLut = nullptr;

  EnergiesStorage energiesStorage{};
};

} // namespace seissol::writer

#endif // ENERGYOUTPUT_H
