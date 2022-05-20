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
  void init(GlobalData* newGlobal,
            seissol::initializers::DynamicRupture* newDynRup,
            seissol::initializers::LTSTree* newDynRuptTree,
            MeshReader* newMeshReader,
            seissol::initializers::LTSTree* newLtsTree,
            seissol::initializers::LTS* newLts,
            seissol::initializers::Lut* newLtsLut,
            bool newIsPlasticityEnabled,
            bool newIsTerminalOutputEnabled,
            const std::string& outputFileNamePrefix,
            double newSyncPointInterval);

  void syncPoint(double time) override;

  void simulationStart() override;

  private:
  real computeStaticWork(const real* degreesOfFreedomPlus,
                         const real* degreesOfFreedomMinus,
                         DRFaceInformation const& faceInfo,
                         DRGodunovData const& godunovData,
                         const real slip[seissol::tensor::slipInterpolated::size()]);

  void computeDynamicRuptureEnergies();

  void computeEnergies();

  void reduceEnergies();

  void printEnergies();

  void writeHeader();

  void writeEnergies(double time);

  bool isEnabled = false;
  bool isTerminalOutputEnabled = false;
  bool isFileOutputEnabled = false;
  bool isPlasticityEnabled = false;

  std::string outputFileName;
  std::ofstream out;

  const GlobalData* global = nullptr;
  seissol::initializers::DynamicRupture* dynRup = nullptr;
  seissol::initializers::LTSTree* dynRupTree = nullptr;
  MeshReader* meshReader = nullptr;
  seissol::initializers::LTSTree* ltsTree = nullptr;
  seissol::initializers::LTS* lts = nullptr;
  seissol::initializers::Lut* ltsLut = nullptr;

  EnergiesStorage energiesStorage{};
};

} // namespace seissol::writer

#endif // ENERGYOUTPUT_H
