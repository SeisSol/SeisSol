#ifndef ENERGYOUTPUT_H
#define ENERGYOUTPUT_H

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
  std::array<double, 8> energies{};

  double& gravitationalEnergy() {
    return energies[0];
  }

  double& acousticEnergy() {
    return energies[1];
  }

  double& acousticKineticEnergy() {
    return energies[2];
  }

  double& elasticEnergy() {
    return energies[3];
  }

  double& elasticKineticEnergy() {
    return energies[4];
  }

  double& totalFrictionalWork() {
    return energies[5];
  }

  double& staticFrictionalWork() {
    return energies[6];
  }

  double& plasticMoment() {
    return energies[7];
  }


};

class EnergyOutput : public Module {
public:
  void init(
      GlobalData* newGlobal,
      seissol::initializers::DynamicRupture* newDynRup,
      seissol::initializers::LTSTree* newDynRuptTree,
      MeshReader* newMeshReader,
      seissol::initializers::LTSTree* newLtsTree,
      seissol::initializers::LTS* newLts,
      seissol::initializers::Lut* newLtsLut,
      bool newUsePlasticity,
      double newSyncPointInterval) {

    global = newGlobal;
    dynRup = newDynRup;
    dynRupTree = newDynRuptTree;
    meshReader = newMeshReader;
    ltsTree = newLtsTree;
    lts = newLts;
    ltsLut = newLtsLut;
    usePlasticity = newUsePlasticity;

    Modules::registerHook(*this, SIMULATION_START);
    Modules::registerHook(*this, SYNCHRONIZATION_POINT);
    setSyncInterval(newSyncPointInterval);
  }

  void syncPoint(double time) override {
    logInfo() << "Energies at time" << time;
    computeEnergies();
    reduceEnergies();
    printEnergies();
  }

private:
  real computePlasticMoment();

  real computeStaticWork(
                         real* degreesOfFreedomPlus,
                         real* degreesOfFreedomMinus,
                         DRFaceInformation const& faceInfo,
                         DRGodunovData const& godunovData,
                         real slip[seissol::tensor::slipInterpolated::size()]);

  void computeDynamicRuptureEnergies();

  void computeEnergies();

  void reduceEnergies();

  void printEnergies();

  const GlobalData* global;
  seissol::initializers::DynamicRupture* dynRup;
  seissol::initializers::LTSTree* dynRupTree;
  MeshReader* meshReader;
  seissol::initializers::LTSTree* ltsTree;
  seissol::initializers::LTS* lts;
  seissol::initializers::Lut* ltsLut;
  bool usePlasticity;

  EnergiesStorage energiesStorage{};
};

} // namespace seissol::writer

#endif // ENERGYOUTPUT_H
