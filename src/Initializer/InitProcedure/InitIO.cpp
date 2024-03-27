#include "Init.hpp"
#include "InitIO.hpp"
#include "Initializer/BasicTypedefs.hpp"
#include <IO/Instance/Mesh/VtkHdf.hpp>
#include <IO/Writer/Writer.hpp>
#include <Numerical_aux/Transformation.h>
#include <SeisSol.h>
#include <cstring>
#include <kernel.h>
#include <string>
#include <tensor.h>
#include <vector>
#include "DynamicRupture/Misc.h"
#include "Common/filesystem.h"

#include "Parallel/MPI.h"

namespace {

static void setupCheckpointing(seissol::SeisSol& seissolInstance) {
  auto& checkpoint = seissolInstance.getOutputManager().getCheckpointManager();

  {
    auto* tree = seissolInstance.getMemoryManager().getLtsTree();
    checkpoint.registerTree<void>("lts", tree, nullptr);
    seissolInstance.getMemoryManager().getLts()->registerCheckpointVariables(checkpoint, tree);
  }

  {
    auto* tree = seissolInstance.getMemoryManager().getDynamicRuptureTree();
    checkpoint.registerTree<void>("dynrup", tree, nullptr);
    seissolInstance.getMemoryManager().getDynamicRupture()->registerCheckpointVariables(checkpoint,
                                                                                        tree);
  }

  {
    auto* tree = seissolInstance.getMemoryManager().getBoundaryTree();
    checkpoint.registerTree<void>("boundary", tree, nullptr);
    seissolInstance.getMemoryManager().getBoundary()->registerCheckpointVariables(checkpoint, tree);
  }

  // seissolInstance.getOutputManager().loadCheckpoint(
  //     seissolInstance.getSeisSolParameters().output.checkpointParameters.fileName);

  if (seissolInstance.getSeisSolParameters().output.checkpointParameters.enabled) {
    // TODO: for now, allow only _one_ checkpoint interval which checkpoints everything existent
    seissolInstance.getOutputManager().setupCheckpoint(
        seissolInstance.getSeisSolParameters().output.checkpointParameters.fileName,
        seissolInstance.getSeisSolParameters().output.checkpointParameters.interval);
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

  if (seissolParams.output.waveFieldParameters.enabled &&
      seissolParams.output.waveFieldParameters.vtkorder < 0) {
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

  if (seissolParams.output.waveFieldParameters.enabled &&
      seissolParams.output.waveFieldParameters.vtkorder >= 0) {
    auto order = seissolParams.output.waveFieldParameters.vtkorder;
    auto& meshReader = seissolInstance.meshReader();
    io::writer::ScheduledWriter schedWriter;
    schedWriter.name = "wavefield";
    schedWriter.interval = seissolParams.output.waveFieldParameters.interval;
    auto writer = io::instance::mesh::VtkHdfWriter(
        "wavefield", seissolInstance.meshReader().getElements().size(), 3, order);
    writer.addPointProjector([=](double* target, std::size_t index) {
      const auto& element = meshReader.getElements()[index];
      const auto& vertexArray = meshReader.getVertices();

      // for the very time being, circumvent the bounding box mechanism of Yateto as follows.
      const double zero[3] = {0, 0, 0};
      seissol::transformations::tetrahedronReferenceToGlobal(
          vertexArray[element.vertices[0]].coords,
          vertexArray[element.vertices[1]].coords,
          vertexArray[element.vertices[2]].coords,
          vertexArray[element.vertices[3]].coords,
          zero,
          &target[0]);
      for (std::size_t i = 1; i < tensor::vtk3d::Shape[order][1]; ++i) {
        seissol::transformations::tetrahedronReferenceToGlobal(
            vertexArray[element.vertices[0]].coords,
            vertexArray[element.vertices[1]].coords,
            vertexArray[element.vertices[2]].coords,
            vertexArray[element.vertices[3]].coords,
            &init::vtk3d::Values[order][i * 3 - 3],
            &target[i * 3]);
      }
    });
    std::vector<std::string> quantityLabels = {
        "sigma_xx",
        "sigma_yy",
        "sigma_zz",
        "sigma_xy",
        "sigma_yz",
        "sigma_xz",
        "u",
        "v",
        "w",
#ifdef USE_POROELASTIC
        "p",
        "u_f",
        "v_f",
        "w_f",
#endif
    };
    std::vector<std::string> plasticityLabels = {
        "ep_xx", "ep_yy", "ep_zz", "ep_xy", "ep_yz", "ep_xz", "eta"};
    for (std::size_t quantity = 0; quantity < NUMBER_OF_QUANTITIES; ++quantity) {
      if (seissolParams.output.waveFieldParameters.outputMask[quantity]) {
        writer.addPointData<real>(
            quantityLabels[quantity], {}, [=](real* target, std::size_t index) {
              const auto* dofsAllQuantities = ltsLut->lookup(lts->dofs, index);
              const auto* dofsSingleQuantity =
                  dofsAllQuantities + NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * quantity;
              kernel::projectBasisToVtkVolume vtkproj;
              vtkproj.qb = dofsSingleQuantity;
              vtkproj.xv(order) = target;
              vtkproj.collvv(CONVERGENCE_ORDER, order) =
                  init::collvv::Values[CONVERGENCE_ORDER + (CONVERGENCE_ORDER + 1) * order];
              vtkproj.execute(order);
            });
      }
    }
    if (seissolParams.model.plasticity) {
      for (std::size_t quantity = 0; quantity < 7; ++quantity) {
        if (seissolParams.output.waveFieldParameters.plasticityMask[quantity]) {
          writer.addPointData<real>(
              plasticityLabels[quantity], {}, [=](real* target, std::size_t index) {
                const auto* dofsAllQuantities = ltsLut->lookup(lts->pstrain, index);
                const auto* dofsSingleQuantity =
                    dofsAllQuantities + NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * quantity;
                kernel::projectBasisToVtkVolume vtkproj;
                vtkproj.qb = dofsSingleQuantity;
                vtkproj.xv(order) = target;
                vtkproj.collvv(CONVERGENCE_ORDER, order) =
                    init::collvv::Values[CONVERGENCE_ORDER + (CONVERGENCE_ORDER + 1) * order];
                vtkproj.execute(order);
              });
        }
      }
    }
    schedWriter.planWrite = writer.makeWriter();
    seissolInstance.getOutputManager().addOutput(schedWriter);
  }

  if (seissolParams.output.freeSurfaceParameters.enabled &&
      seissolParams.output.freeSurfaceParameters.vtkorder < 0) {
    // Initialize free surface output
    seissolInstance.freeSurfaceWriter().init(seissolInstance.meshReader(),
                                             &seissolInstance.freeSurfaceIntegrator(),
                                             seissolParams.output.prefix.c_str(),
                                             seissolParams.output.freeSurfaceParameters.interval,
                                             seissolParams.output.xdmfWriterBackend,
                                             backupTimeStamp);
  }

  if (seissolParams.output.freeSurfaceParameters.enabled &&
      seissolParams.output.freeSurfaceParameters.vtkorder >= 0) {
    logWarning() << "High-order free surface output has not yet been implemented. The output will "
                    "be disabled.";
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
  if (seissolParams.output.waveFieldParameters.enabled &&
      seissolParams.output.waveFieldParameters.vtkorder < 0) {
    seissolInstance.waveFieldWriter().enable();
    seissolInstance.waveFieldWriter().setFilename(seissolParams.output.prefix.c_str());
    seissolInstance.waveFieldWriter().setWaveFieldInterval(
        seissolParams.output.waveFieldParameters.interval);
  }
}

static void enableFreeSurfaceOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
  if (seissolParams.output.freeSurfaceParameters.enabled &&
      seissolParams.output.freeSurfaceParameters.vtkorder < 0) {
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
