#include "Init.hpp"
#include "InitIO.hpp"
#include "Initializer/BasicTypedefs.hpp"
#include <IO/Instance/CheckpointWriter.hpp>
#include <IO/Writer/Writer.hpp>
#include <SeisSol.h>
#include <cstring>
#include <vector>
#include "DynamicRupture/Misc.h"
#include "Common/filesystem.h"

#include "Parallel/MPI.h"

namespace {

static void setupCheckpointing(seissol::SeisSol& seissolInstance) {
  seissolInstance.getOutputManager().loadCheckpoint(seissolInstance.getSeisSolParameters().output.checkpointParameters.fileName);

  if (seissolInstance.getSeisSolParameters().output.checkpointParameters.enabled) {
    // TODO: for now, allow only _one_ checkpoint interval which checkpoints everything existent
    seissolInstance.getOutputManager().setupCheckpoint(seissolInstance.getSeisSolParameters().output.checkpointParameters.fileName, seissolInstance.getSeisSolParameters().output.checkpointParameters.interval);
  }

  /*if (hasCheckpoint) {
    seissolInstance.simulator().setCurrentTime(seissolInstanc7e.checkPointManager().header().time());
    seissolInstance.faultWriter().setTimestep(faultTimeStep);
  }*/
}

static void setupOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
  auto* lts = memoryManager.getLts();
  auto* ltsTree = memoryManager.getLtsTree();
  auto* ltsLut = memoryManager.getLtsLut();
  auto* dynRup = memoryManager.getDynamicRupture();
  auto* dynRupTree = memoryManager.getDynamicRuptureTree();
  auto* globalData = memoryManager.getGlobalDataOnHost();
  const auto& backupTimeStamp = seissolInstance.getBackupTimeStamp();

  constexpr auto numberOfQuantities =
      tensor::Q::Shape[sizeof(tensor::Q::Shape) / sizeof(tensor::Q::Shape[0]) - 1];
  // TODO(David): handle attenuation properly here. We'll probably not want it to be contained in
  // numberOfQuantities. But the compile-time parameter NUMBER_OF_QUANTITIES contains it
  // nonetheless.

  if (seissolParams.output.waveFieldParameters.enabled) {
    // seissolInstance.getOutputManager().addOutput(WaveField);
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
        CONVERGENCE_ORDER,
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

static void initFaultOutputManager(seissol::SeisSol& seissolInstance) {
  const auto& backupTimeStamp = seissolInstance.getBackupTimeStamp();
  seissolInstance.getMemoryManager().initFaultOutputManager(backupTimeStamp);
  auto* faultOutputManager = seissolInstance.getMemoryManager().getFaultOutputManager();
  seissolInstance.timeManager().setFaultOutputManager(faultOutputManager);

  seissolInstance.getMemoryManager().getFaultOutputManager()->initFaceToLtsMap();
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
    MPI::mpi.barrier(MPI::mpi.comm());
  }

  enableWaveFieldOutput(seissolInstance);
  setIntegralMask(seissolInstance);
  enableFreeSurfaceOutput(seissolInstance);
  initFaultOutputManager(seissolInstance);
  setupCheckpointing(seissolInstance);
  setupOutput(seissolInstance);
  logInfo(rank) << "End init output.";
}
