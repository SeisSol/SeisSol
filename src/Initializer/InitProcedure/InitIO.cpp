#include "Init.hpp"
#include "InitIO.hpp"
#include "Initializer/BasicTypedefs.hpp"
#include <SeisSol.h>
#include <cstring>
#include <vector>
#include "DynamicRupture/Misc.h"

#include "Parallel/MPI.h"

static void setupCheckpointing() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();

  auto* lts = memoryManager.getLts();
  auto* ltsTree = memoryManager.getLtsTree();
  auto* dynRup = memoryManager.getDynamicRupture();
  auto* dynRupTree = memoryManager.getDynamicRuptureTree();

  // Initialize checkpointing
  int faultTimeStep;

  // Only R&S friction explicitly stores the state variable, otherwise use the accumulated slip
  // magnitude
  real* stateVariable{nullptr};
  if (dynamic_cast<seissol::initializers::LTSRateAndState*>(dynRup)) {
    stateVariable = reinterpret_cast<real*>(dynRupTree->var(
        dynamic_cast<seissol::initializers::LTSRateAndState*>(dynRup)->stateVariable));
  } else {
    stateVariable = reinterpret_cast<real*>(dynRupTree->var(dynRup->accumulatedSlipMagnitude));
  }
  // Only with prakash-clifton regularization, we store the fault strength, otherwise use the
  // friction coefficient
  real* strength{nullptr};
  if (dynamic_cast<seissol::initializers::LTSLinearSlipWeakeningBimaterial*>(dynRup)) {
    stateVariable = reinterpret_cast<real*>(dynRupTree->var(
        dynamic_cast<seissol::initializers::LTSLinearSlipWeakeningBimaterial*>(dynRup)
            ->regularisedStrength));
  } else {
    stateVariable = reinterpret_cast<real*>(dynRupTree->var(dynRup->mu));
  }

  size_t numSides = seissol::SeisSol::main.meshReader().getFault().size();
  unsigned int numBndGP = seissol::dr::misc::numberOfBoundaryGaussPoints;

  bool hasCheckpoint = seissol::SeisSol::main.checkPointManager().init(
      reinterpret_cast<real*>(ltsTree->var(lts->dofs)),
      ltsTree->getNumberOfCells(lts->dofs.mask) * tensor::Q::size(),
      reinterpret_cast<real*>(dynRupTree->var(dynRup->mu)),
      reinterpret_cast<real*>(dynRupTree->var(dynRup->slipRate1)),
      reinterpret_cast<real*>(dynRupTree->var(dynRup->slipRate2)),
      reinterpret_cast<real*>(dynRupTree->var(dynRup->accumulatedSlipMagnitude)),
      reinterpret_cast<real*>(dynRupTree->var(dynRup->slip1)),
      reinterpret_cast<real*>(dynRupTree->var(dynRup->slip2)),
      stateVariable,
      strength,
      numSides,
      numBndGP,
      faultTimeStep);
  if (hasCheckpoint) {
    seissol::SeisSol::main.simulator().setCurrentTime(
        seissol::SeisSol::main.checkPointManager().header().time());
    seissol::SeisSol::main.faultWriter().setTimestep(faultTimeStep);
  }
}

static void setupOutput() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();
  auto* lts = memoryManager.getLts();
  auto* ltsTree = memoryManager.getLtsTree();
  auto* ltsLut = memoryManager.getLtsLut();
  auto* dynRup = memoryManager.getDynamicRupture();
  auto* dynRupTree = memoryManager.getDynamicRuptureTree();
  auto* globalData = memoryManager.getGlobalDataOnHost();

  constexpr auto numberOfQuantities =
      tensor::Q::Shape[sizeof(tensor::Q::Shape) / sizeof(tensor::Q::Shape[0]) - 1];
  // TODO(David): handle attenuation properly here. We'll probably not want it to be contained in
  // numberOfQuantities. But the compile-time parameter NUMBER_OF_QUANTITIES contains it
  // nonetheless.

  if (seissolParams.output.waveFieldParameters.enabled) {
    // record the clustering info i.e., distribution of elements within an LTS tree
    const std::vector<Element>& meshElements = seissol::SeisSol::main.meshReader().getElements();
    std::vector<unsigned> ltsClusteringData(meshElements.size());
    auto& ltsLayout = seissol::SeisSol::main.getLtsLayout();
    for (const auto& element : meshElements) {
      ltsClusteringData[element.localId] = ltsLayout.getGlobalClusterId(element.localId);
    }
    // Initialize wave field output
    seissol::SeisSol::main.waveFieldWriter().init(
        numberOfQuantities,
        seissol::ConvergenceOrder,
        seissol::kernels::NumberOfAlignedBasisFunctions(),
        seissol::SeisSol::main.meshReader(),
        ltsClusteringData,
        reinterpret_cast<const real*>(ltsTree->var(lts->dofs)),
        reinterpret_cast<const real*>(ltsTree->var(lts->pstrain)),
        seissol::SeisSol::main.postProcessor().getIntegrals(ltsTree),
        ltsLut->getMeshToLtsLut(lts->dofs.mask)[0],
        seissolParams.output.waveFieldParameters,
        seissolParams.output.xdmfWriterBackend);
  }

  if (seissolParams.output.freeSurfaceParameters.enabled) {
    // Initialize free surface output
    seissol::SeisSol::main.freeSurfaceWriter().init(
        seissol::SeisSol::main.meshReader(),
        &seissol::SeisSol::main.freeSurfaceIntegrator(),
        seissolParams.output.prefix.c_str(),
        seissolParams.output.freeSurfaceParameters.interval,
        seissolParams.output.xdmfWriterBackend);
  }

  if (seissolParams.output.receiverParameters.enabled) {
    auto& receiverWriter = seissol::SeisSol::main.receiverWriter();
    // Initialize receiver output
    receiverWriter.init(seissolParams.output.prefix, seissolParams.output.receiverParameters);
    receiverWriter.addPoints(seissol::SeisSol::main.meshReader(), *ltsLut, *lts, globalData);
    seissol::SeisSol::main.timeManager().setReceiverClusters(receiverWriter);
  }

  if (seissolParams.output.energyParameters.enabled) {
    auto& energyOutput = seissol::SeisSol::main.energyOutput();

    energyOutput.init(globalData,
                      dynRup,
                      dynRupTree,
                      &seissol::SeisSol::main.meshReader(),
                      ltsTree,
                      lts,
                      ltsLut,
                      seissolParams.model.plasticity,
                      seissolParams.output.prefix,
                      seissolParams.output.energyParameters);
  }

  seissol::SeisSol::main.flopCounter().init(seissolParams.output.prefix.c_str());

  seissol::SeisSol::main.analysisWriter().init(&seissol::SeisSol::main.meshReader(),
                                               seissolParams.output.prefix.c_str());
}

static void enableCheckpointing() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  if (seissolParams.output.checkpointParameters.enabled) {
    seissol::SeisSol::main.simulator().setCheckPointInterval(
        seissolParams.output.checkpointParameters.interval);
    seissol::SeisSol::main.checkPointManager().setBackend(
        seissolParams.output.checkpointParameters.backend);
    seissol::SeisSol::main.checkPointManager().setFilename(
        seissolParams.output.checkpointParameters.fileName.c_str());
  }
}

static void initFaultOutputManager() {
  seissol::SeisSol::main.getMemoryManager().initFaultOutputManager();

  auto* faultOutputManager = seissol::SeisSol::main.getMemoryManager().getFaultOutputManager();
  seissol::SeisSol::main.timeManager().setFaultOutputManager(faultOutputManager);

  seissol::SeisSol::main.getMemoryManager().getFaultOutputManager()->initFaceToLtsMap();
}

static void enableWaveFieldOutput() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  if (seissolParams.output.waveFieldParameters.enabled) {
    seissol::SeisSol::main.waveFieldWriter().enable();
    seissol::SeisSol::main.waveFieldWriter().setFilename(seissolParams.output.prefix.c_str());
    seissol::SeisSol::main.waveFieldWriter().setWaveFieldInterval(
        seissolParams.output.waveFieldParameters.interval);
  }
}

static void enableFreeSurfaceOutput() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();
  if (seissolParams.output.freeSurfaceParameters.enabled) {
    seissol::SeisSol::main.freeSurfaceWriter().enable();

    seissol::SeisSol::main.freeSurfaceIntegrator().initialize(
        seissolParams.output.freeSurfaceParameters.refinement,
        memoryManager.getGlobalDataOnHost(),
        memoryManager.getLts(),
        memoryManager.getLtsTree(),
        memoryManager.getLtsLut());
  }
}

static void setIntegralMask() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  seissol::SeisSol::main.postProcessor().setIntegrationMask(
      seissolParams.output.waveFieldParameters.integrationMask);
}

void seissol::initializer::initprocedure::initIO() {
  logInfo(seissol::MPI::mpi.rank()) << "Begin init output.";
  // always enable checkpointing first
  enableCheckpointing();
  enableWaveFieldOutput();
  setIntegralMask();
  enableFreeSurfaceOutput();
  initFaultOutputManager();
  setupCheckpointing();
  setupOutput();
  logInfo(seissol::MPI::mpi.rank()) << "End init output.";
}
