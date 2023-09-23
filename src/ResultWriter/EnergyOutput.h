#ifndef ENERGYOUTPUT_H
#define ENERGYOUTPUT_H

#include <Initializer/tree/LTSForest.hpp>
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
  std::array<double, 9> energies{};

  double& gravitationalEnergy();

  double& acousticEnergy();

  double& acousticKineticEnergy();

  double& elasticEnergy();

  double& elasticKineticEnergy();

  double& totalFrictionalWork();

  double& staticFrictionalWork();

  double& plasticMoment();

  double& seismicMoment();
};

class EnergyOutput : public Module {
  public:
  void init(seissol::initializers::GlobalDataStorage* newGlobal,
            seissol::initializers::DynRupLTSForest* dynrupForest,
            seissol::geometry::MeshReader* newMeshReader,
            seissol::initializers::ClusterLTSForest* clusterForest,
            seissol::initializers::ClusterBackmap* backmap,
            bool newIsPlasticityEnabled,
            const std::string& outputFileNamePrefix,
            const seissol::initializer::parameters::EnergyOutputParameters& parameters);

  void syncPoint(double time) override;

  void simulationStart() override;

  private:
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

  seissol::initializers::GlobalDataStorage* global = nullptr;
  seissol::initializers::DynRupLTSForest* dynrupForest = nullptr;
  seissol::geometry::MeshReader* meshReader = nullptr;
  seissol::initializers::ClusterLTSForest* clusterForest = nullptr;
  seissol::initializers::ClusterBackmap* backmap = nullptr;

  EnergiesStorage energiesStorage{};
};

} // namespace seissol::writer

#endif // ENERGYOUTPUT_H
