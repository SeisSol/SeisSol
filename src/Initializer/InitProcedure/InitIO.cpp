#include "InitIO.hpp"
#include "Common/filesystem.h"
#include "DynamicRupture/Misc.h"
#include "Init.hpp"
#include "Initializer/BasicTypedefs.hpp"
#include "SeisSol.h"
#include <cstring>
#include <vector>

#include "Parallel/MPI.h"

namespace {

// static void setupCheckpointing(seissol::SeisSol& seissolInstance) {
//   const auto& seissolParams = seissolInstance.getSeisSolParameters();
//   auto& memoryManager = seissolInstance.getMemoryManager();

//   auto* lts = memoryManager.getLts();
//   auto* ltsTree = memoryManager.getLtsTree();
//   auto dynRup = memoryManager.getDynamicRupture();
//   auto dynRupTree = memoryManager.getDynamicRuptureTree();

//   // Initialize checkpointing
//   int faultTimeStep;

//   // Only R&S friction explicitly stores the state variable, otherwise use the accumulated slip
//   // magnitude
//   real* stateVariable{nullptr};
//   if (dynamic_cast<seissol::initializer::LTSRateAndState*>(dynRup)) {
//     stateVariable = reinterpret_cast<real*>(dynRupTree->var(
//         dynamic_cast<seissol::initializer::LTSRateAndState*>(dynRup)->stateVariable));
//   } else {
//     stateVariable = reinterpret_cast<real*>(dynRupTree->var(dynRup->accumulatedSlipMagnitude));
//   }
//   // Only with prakash-clifton regularization, we store the fault strength, otherwise use the
//   // friction coefficient
//   real* strength{nullptr};
//   if (dynamic_cast<seissol::initializer::LTSLinearSlipWeakeningBimaterial*>(dynRup)) {
//     stateVariable = reinterpret_cast<real*>(dynRupTree->var(
//         dynamic_cast<seissol::initializer::LTSLinearSlipWeakeningBimaterial*>(dynRup)
//             ->regularisedStrength));
//   } else {
//     stateVariable = reinterpret_cast<real*>(dynRupTree->var(dynRup->mu));
//   }

//   size_t numSides = seissolInstance.meshReader().getFault().size();
//   unsigned int numBndGP = seissol::dr::misc::numberOfBoundaryGaussPoints;

//   bool hasCheckpoint = seissolInstance.checkPointManager().init(
//       reinterpret_cast<real*>(ltsTree->var(lts->dofs)),
//       ltsTree->getNumberOfCells(lts->dofs.mask) * tensor::Q::size(),
//       reinterpret_cast<real*>(dynRupTree->var(dynRup->mu)),
//       reinterpret_cast<real*>(dynRupTree->var(dynRup->slipRate1)),
//       reinterpret_cast<real*>(dynRupTree->var(dynRup->slipRate2)),
//       reinterpret_cast<real*>(dynRupTree->var(dynRup->accumulatedSlipMagnitude)),
//       reinterpret_cast<real*>(dynRupTree->var(dynRup->slip1)),
//       reinterpret_cast<real*>(dynRupTree->var(dynRup->slip2)),
//       stateVariable,
//       strength,
//       numSides,
//       numBndGP,
//       faultTimeStep);
//   if (hasCheckpoint) {
//     seissolInstance.simulator().setCurrentTime(seissolInstance.checkPointManager().header().time());
//     seissolInstance.faultWriter().setTimestep(faultTimeStep);
//   }
// }

static void setupOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
  auto* lts = memoryManager.getLts();
  auto* ltsTree = memoryManager.getLtsTree();
  auto* ltsLut = memoryManager.getLtsLut();
  auto dynRup = memoryManager.getDynamicRupture();
  auto dynRupTree = memoryManager.getDynamicRuptureTree();
  auto* globalData = memoryManager.getGlobalDataOnHost();
  const auto& backupTimeStamp = seissolInstance.getBackupTimeStamp();

  constexpr auto numberOfQuantities =
      tensor::Q::Shape[sizeof(tensor::Q::Shape) / sizeof(tensor::Q::Shape[0]) - 1];
  // TODO(David): handle attenuation properly here. We'll probably not want it to be contained in
  // numberOfQuantities. But the compile-time parameter
  // seissol::model::Material_t::NumberOfQuantities contains it nonetheless.

  if (seissolParams.output.waveFieldParameters.enabled) {
    // record the clustering info i.e., distribution of elements within an LTS tree
    const std::vector<Element>& meshElements = seissolInstance.meshReader().getElements();
    std::vector<unsigned> ltsClusteringData(meshElements.size());
    auto& ltsLayout = seissolInstance.getLtsLayout();
    for (const auto& element : meshElements) {
      ltsClusteringData[element.localId] = ltsLayout.getGlobalClusterId(element.localId);
    }
    // Initialize wave field output
    seissolInstance.waveFieldWriter().init(
        numberOfQuantities,
        ConvergenceOrder,
        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
        seissolInstance.meshReader(),
        ltsClusteringData,
        reinterpret_cast<const real*>(ltsTree->var(lts->dofs)),
        reinterpret_cast<const real*>(ltsTree->var(lts->pstrain)),
        seissolInstance.postProcessor().getIntegrals(ltsTree),
        ltsLut->getMeshToLtsLut(lts->dofs.mask)[0],
        seissolParams.output.waveFieldParameters,
        seissolParams.output.xdmfWriterBackend,
        backupTimeStamp);
  }

  if (seissolParams.output.freeSurfaceParameters.enabled) {
    // Initialize free surface output
    seissolInstance.freeSurfaceWriter().init(seissolInstance.meshReader(),
                                             &seissolInstance.freeSurfaceIntegrator(),
                                             seissolParams.output.prefix.c_str(),
                                             seissolParams.output.freeSurfaceParameters.interval,
                                             seissolParams.output.xdmfWriterBackend,
                                             backupTimeStamp);
  }

  if (seissolParams.output.receiverParameters.enabled) {
    auto& receiverWriter = seissolInstance.receiverWriter();
    // Initialize receiver output
    receiverWriter.init(seissolParams.output.prefix,
                        seissolParams.timeStepping.endTime,
                        seissolParams.output.receiverParameters);
    receiverWriter.addPoints(seissolInstance.meshReader(), *ltsLut, *lts, globalData);
    seissolInstance.timeManager().setReceiverClusters(receiverWriter);
  }

  if (seissolParams.output.energyParameters.enabled) {
    auto& energyOutput = seissolInstance.energyOutput();

    energyOutput.init(globalData,
                      dynRup,
                      dynRupTree,
                      &seissolInstance.meshReader(),
                      ltsTree,
                      lts,
                      ltsLut,
                      seissolParams.model.plasticity,
                      seissolParams.output.prefix,
                      seissolParams.output.energyParameters);
  }

  seissolInstance.flopCounter().init(seissolParams.output.prefix.c_str());

  seissolInstance.analysisWriter().init(&seissolInstance.meshReader(),
                                        seissolParams.output.prefix.c_str());
}

static void enableCheckpointing(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  if (seissolParams.output.checkpointParameters.enabled) {
    seissolInstance.simulator().setCheckPointInterval(
        seissolParams.output.checkpointParameters.interval);
    seissolInstance.checkPointManager().setBackend(
        seissolParams.output.checkpointParameters.backend);
    seissolInstance.checkPointManager().setFilename(
        seissolParams.output.checkpointParameters.fileName.c_str());
  }
}

static void initFaultOutputManager(seissol::SeisSol& seissolInstance) {
  const auto& backupTimeStamp = seissolInstance.getBackupTimeStamp();
  seissolInstance.getMemoryManager().initFaultOutputManager(backupTimeStamp);
  auto faultOutputManager = seissolInstance.getMemoryManager().getFaultOutputManager();
  seissolInstance.timeManager().setFaultOutputManager(faultOutputManager);
  for (unsigned int i = 0; i < MULTIPLE_SIMULATIONS; i++) {
    seissolInstance.getMemoryManager().getFaultOutputManager()[i]->initFaceToLtsMap();
  }
}

static void enableWaveFieldOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  if (seissolParams.output.waveFieldParameters.enabled) {
    seissolInstance.waveFieldWriter().enable();
    seissolInstance.waveFieldWriter().setFilename(seissolParams.output.prefix.c_str());
    seissolInstance.waveFieldWriter().setWaveFieldInterval(
        seissolParams.output.waveFieldParameters.interval);
  }
}

static void enableFreeSurfaceOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
  if (seissolParams.output.freeSurfaceParameters.enabled) {
    seissolInstance.freeSurfaceWriter().enable();

    seissolInstance.freeSurfaceIntegrator().initialize(
        seissolParams.output.freeSurfaceParameters.refinement,
        memoryManager.getGlobalDataOnHost(),
        memoryManager.getLts(),
        memoryManager.getLtsTree(),
        memoryManager.getLtsLut());
  }
}

static void setIntegralMask(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  seissolInstance.postProcessor().setIntegrationMask(
      seissolParams.output.waveFieldParameters.integrationMask);
}

} // namespace

void seissol::initializer::initprocedure::initIO(seissol::SeisSol& seissolInstance) {
  const auto rank = MPI::mpi.rank();
  logInfo(rank) << "Begin init output.";

  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const filesystem::path outputPath(seissolParams.output.prefix);
  const auto outputDir = filesystem::directory_entry(outputPath.parent_path());
  if (!filesystem::exists(outputDir)) {
    logWarning(rank) << "Output directory does not exist yet. We therefore create it now.";
    if (rank == 0) {
      filesystem::create_directory(outputDir);
    }
  }
  MPI::mpi.barrier(MPI::mpi.comm());

  // always enable checkpointing first
  // enableCheckpointing(seissolInstance);
  enableWaveFieldOutput(seissolInstance);
  setIntegralMask(seissolInstance);
  enableFreeSurfaceOutput(seissolInstance);
  initFaultOutputManager(seissolInstance);
  // setupCheckpointing(seissolInstance);
  setupOutput(seissolInstance);
  logInfo(rank) << "End init output.";
}
