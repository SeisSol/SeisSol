#include "InstantaneousTimeMirrorManager.h"
#include "Modules/Modules.h"
#include "Initializer/CellLocalMatrices.h"


namespace seissol {

void InstantaneousTimeMirrorManager::init(
    double velocityScalingFactor,
    double triggerTime,
    MeshReader* meshReader,
    initializers::LTSTree* ltsTree,
    initializers::LTS* lts,
    initializers::Lut* ltsLut) {
  isEnabled = true;
  this->velocityScalingFactor = velocityScalingFactor;
  this->triggerTime = triggerTime;
  this->meshReader = meshReader;
  this->ltsTree = ltsTree;
  this->lts = lts;
  this->ltsLut = ltsLut;
  setSyncInterval(triggerTime);
  Modules::registerHook(*this, SYNCHRONIZATION_POINT);
}

void InstantaneousTimeMirrorManager::syncPoint(double currentTime) {
  Module::syncPoint(currentTime);


  const auto rank = MPI::mpi.rank();
  logInfo(rank) << "InstantaneousTimeMirrorManager: Factor " << velocityScalingFactor;
  if (!isEnabled) {
    logInfo(rank) << "InstantaneousTimeMirrorManager: Skipping syncing at " << currentTime << "as it is disabled";
    return;
  }

  logInfo(rank) << "InstantaneousTimeMirrorManager Syncing at " << currentTime;

  logInfo(rank) << "Scaling velocitites by factor of " << velocityScalingFactor;
  updateVelocities();

  logInfo(rank) << "Updating CellLocalMatrices";
  initializers::initializeCellLocalMatrices(*meshReader, ltsTree, lts, ltsLut);

  logInfo(rank) << "Finished flipping.";
  isEnabled = false;
}

void InstantaneousTimeMirrorManager::updateVelocities() {
  for (auto it = ltsTree->beginLeaf(initializers::LayerMask(Ghost)); it != ltsTree->endLeaf(); ++it) {
    CellMaterialData* materials = it->var(lts->material);
    for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
      auto& material = materials[cell];
      material.local.rho *= this->velocityScalingFactor * this->velocityScalingFactor;
    }
  }

}

void initializeTimeMirrorManagers(
    double scalingFactor, double triggerTime,
    MeshReader* meshReader,
    initializers::LTSTree* ltsTree,
    initializers::LTS* lts,
    initializers::Lut* ltsLut,
    InstantaneousTimeMirrorManager& increaseManager,
    InstantaneousTimeMirrorManager& decreaseManager
) {
  increaseManager.init(scalingFactor, triggerTime, meshReader, ltsTree, lts, ltsLut);
  const double eps = 1e-7; // TODO(Lukas) Use CFL condition for this
  decreaseManager.init(1 / scalingFactor, triggerTime + eps, meshReader, ltsTree, lts, ltsLut);
};

} // namespace seissol