// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitIO.h"
#include "Common/Filesystem.h"
#include "Equations/Datastructures.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Writer/Writer.h"
#include "Init.h"
#include "Numerical/Transformation.h"
#include "SeisSol.h"
#include <Common/Constants.h>
#include <Geometry/MeshDefinition.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Tree/Layer.h>
#include <Model/Plasticity.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include <algorithm>
#include <cstring>
#include <init.h>
#include <kernel.h>
#include <string>
#include <tensor.h>
#include <utils/logger.h>
#include <vector>

#include "Parallel/MPI.h"

namespace {

void setupCheckpointing(seissol::SeisSol& seissolInstance) {
  auto& checkpoint = seissolInstance.getOutputManager().getCheckpointManager();

  {
    auto* tree = seissolInstance.getMemoryManager().getLtsTree();
    std::vector<std::size_t> globalIds(
        tree->getNumberOfCells(seissol::initializer::LayerMask(Ghost)));
    const auto* ltsToMesh = seissolInstance.getMemoryManager().getLtsLut()->getLtsToMeshLut(
        seissol::initializer::LayerMask(Ghost));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < globalIds.size(); ++i) {
      globalIds[i] = seissolInstance.meshReader().getElements()[ltsToMesh[i]].globalId;
    }
    checkpoint.registerTree("lts", tree, globalIds);
    seissolInstance.getMemoryManager().getLts()->registerCheckpointVariables(checkpoint, tree);
  }

  {
    auto* tree = seissolInstance.getMemoryManager().getDynamicRuptureTree();
    auto* dynrup = seissolInstance.getMemoryManager().getDynamicRupture();
    std::vector<std::size_t> faceIdentifiers(
        tree->getNumberOfCells(seissol::initializer::LayerMask(Ghost)));
    const auto* drFaceInformation = tree->var(dynrup->faceInformation);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
      auto faultFace = drFaceInformation[i].meshFace;
      const auto& fault = seissolInstance.meshReader().getFault()[faultFace];
      // take the positive cell and side as fault face identifier
      // (should result in roughly twice as large numbers as when indexing all faces; cf. handshake
      // theorem)
      faceIdentifiers[i] = fault.globalId * 4 + fault.side;
    }
    checkpoint.registerTree("dynrup", tree, faceIdentifiers);
    dynrup->registerCheckpointVariables(checkpoint, tree);
  }

  if (seissolInstance.getCheckpointLoadFile().has_value()) {
    const double time = seissolInstance.getOutputManager().loadCheckpoint(
        seissolInstance.getCheckpointLoadFile().value());
    seissolInstance.simulator().setCurrentTime(time);
  }

  if (seissolInstance.getSeisSolParameters().output.checkpointParameters.enabled) {
    // FIXME: for now, we allow only _one_ checkpoint interval which checkpoints everything existent
    seissolInstance.getOutputManager().setupCheckpoint(
        seissolInstance.getSeisSolParameters().output.checkpointParameters.interval);
  }
}

void setupOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
  auto* lts = memoryManager.getLts();
  auto* ltsTree = memoryManager.getLtsTree();
  auto* ltsLut = memoryManager.getLtsLut();
  auto* dynRup = memoryManager.getDynamicRupture();
  auto* dynRupTree = memoryManager.getDynamicRuptureTree();
  auto* globalData = memoryManager.getGlobalDataOnHost();
  const auto& backupTimeStamp = seissolInstance.getBackupTimeStamp();

  constexpr auto NumQuantities =
      tensor::Q::Shape[sizeof(tensor::Q::Shape) / sizeof(tensor::Q::Shape[0]) - 1];
  // TODO(David): handle attenuation properly here. We'll probably not want it to be contained in
  // numberOfQuantities. But the compile-time parameter
  // seissol::model::MaterialT::NumQuantities contains it nonetheless.

  if (seissolParams.output.waveFieldParameters.enabled &&
      seissolParams.output.waveFieldParameters.vtkorder < 0) {
    // record the clustering info i.e., distribution of elements within an LTS tree
    const std::vector<Element>& meshElements = seissolInstance.meshReader().getElements();
    std::vector<unsigned> ltsClusteringData(meshElements.size());
    for (const auto& element : meshElements) {
      ltsClusteringData[element.localId] = element.clusterId;
    }
    // Initialize wave field output
    seissolInstance.waveFieldWriter().init(
        NumQuantities,
        ConvergenceOrder,
        NumAlignedBasisFunctions,
        seissolInstance.meshReader(),
        ltsClusteringData,
        reinterpret_cast<const real*>(ltsTree->var(lts->dofs)),
        reinterpret_cast<const real*>(ltsTree->var(lts->pstrain)),
        seissolInstance.postProcessor().getIntegrals(ltsTree),
        ltsLut->getMeshToLtsLut(lts->dofs.mask)[0].data(),
        seissolParams.output.waveFieldParameters,
        seissolParams.output.xdmfWriterBackend,
        backupTimeStamp);
  }

  // TODO(David): change Yateto/TensorForge interface to make padded sizes more accessible
  constexpr auto QDofSizePadded = tensor::Q::Size / tensor::Q::Shape[1];
  constexpr auto FaceDisplacementPadded =
      tensor::faceDisplacement::Size / tensor::faceDisplacement::Shape[1];

  if (seissolParams.output.waveFieldParameters.enabled &&
      seissolParams.output.waveFieldParameters.vtkorder >= 0) {

    // Effectively temporary code for now. To be refactored.
    if (seissolParams.output.waveFieldParameters.vtkorder == 0) {
      logError() << "VTK order 0 is currently not supported for the wavefield output.";
    }

    auto order = seissolParams.output.waveFieldParameters.vtkorder;
    auto& meshReader = seissolInstance.meshReader();

    // TODO: store somewhere
    std::vector<std::size_t> celllist;
    celllist.reserve(meshReader.getElements().size());
    if (seissolParams.output.waveFieldParameters.bounds.enabled ||
        !seissolParams.output.waveFieldParameters.groups.empty()) {
      const auto& vertexArray = meshReader.getVertices();
      for (std::size_t i = 0; i < meshReader.getElements().size(); ++i) {
        const auto& element = meshReader.getElements()[i];
        const auto& vertex0 = vertexArray[element.vertices[0]].coords;
        const auto& vertex1 = vertexArray[element.vertices[1]].coords;
        const auto& vertex2 = vertexArray[element.vertices[2]].coords;
        const auto& vertex3 = vertexArray[element.vertices[3]].coords;
        const bool inGroup = seissolParams.output.waveFieldParameters.groups.empty() ||
                             seissolParams.output.waveFieldParameters.groups.find(element.group) !=
                                 seissolParams.output.waveFieldParameters.groups.end();
        const bool inRegion = !seissolParams.output.waveFieldParameters.bounds.enabled ||
                              (seissolParams.output.waveFieldParameters.bounds.contains(
                                   vertex0[0], vertex0[1], vertex0[2]) ||
                               seissolParams.output.waveFieldParameters.bounds.contains(
                                   vertex1[0], vertex1[1], vertex1[2]) ||
                               seissolParams.output.waveFieldParameters.bounds.contains(
                                   vertex2[0], vertex2[1], vertex2[2]) ||
                               seissolParams.output.waveFieldParameters.bounds.contains(
                                   vertex3[0], vertex3[1], vertex3[2]));
        if (inGroup && inRegion) {
          celllist.push_back(i);
        }
      }
    } else {
      for (std::size_t i = 0; i < meshReader.getElements().size(); ++i) {
        celllist.push_back(i);
      }
    }
    auto* cellIndices = new std::size_t[celllist.size()];
    std::copy(celllist.begin(), celllist.end(), cellIndices);

    io::writer::ScheduledWriter schedWriter;
    schedWriter.name = "wavefield";
    schedWriter.interval = seissolParams.output.waveFieldParameters.interval;
    auto writer = io::instance::mesh::VtkHdfWriter("wavefield", celllist.size(), 3, order);

    writer.addPointProjector([=](double* target, std::size_t index) {
      const auto& element = meshReader.getElements()[cellIndices[index]];
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
        double point[3] = {init::vtk3d::Values[order][i * 3 - 3 + 0],
                           init::vtk3d::Values[order][i * 3 - 3 + 1],
                           init::vtk3d::Values[order][i * 3 - 3 + 2]};
        seissol::transformations::tetrahedronReferenceToGlobal(
            vertexArray[element.vertices[0]].coords,
            vertexArray[element.vertices[1]].coords,
            vertexArray[element.vertices[2]].coords,
            vertexArray[element.vertices[3]].coords,
            point,
            &target[i * 3]);
      }
    });

    for (std::size_t quantity = 0; quantity < seissol::model::MaterialT::Quantities.size();
         ++quantity) {
      if (seissolParams.output.waveFieldParameters.outputMask[quantity]) {
        writer.addPointData<real>(
            seissol::model::MaterialT::Quantities[quantity],
            {},
            [=](real* target, std::size_t index) {
              const auto* dofsAllQuantities = ltsLut->lookup(lts->dofs, cellIndices[index]);
              const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
              kernel::projectBasisToVtkVolume vtkproj;
              vtkproj.qb = dofsSingleQuantity;
              vtkproj.xv(order) = target;
              vtkproj.collvv(ConvergenceOrder, order) =
                  init::collvv::Values[ConvergenceOrder + (ConvergenceOrder + 1) * order];
              vtkproj.execute(order);
            });
      }
    }
    if (seissolParams.model.plasticity) {
      for (std::size_t quantity = 0; quantity < seissol::model::PlasticityData::Quantities.size();
           ++quantity) {
        if (seissolParams.output.waveFieldParameters.plasticityMask[quantity]) {
          writer.addPointData<real>(
              seissol::model::PlasticityData::Quantities[quantity],
              {},
              [=](real* target, std::size_t index) {
                const auto* dofsAllQuantities = ltsLut->lookup(lts->pstrain, cellIndices[index]);
                const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
                kernel::projectBasisToVtkVolume vtkproj;
                vtkproj.qb = dofsSingleQuantity;
                vtkproj.xv(order) = target;
                vtkproj.collvv(ConvergenceOrder, order) =
                    init::collvv::Values[ConvergenceOrder + (ConvergenceOrder + 1) * order];
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
    // Effectively temporary code for now. To be refactored.
    if (seissolParams.output.freeSurfaceParameters.vtkorder == 0) {
      logError() << "VTK order 0 is currently not supported for the free surface output.";
    }

    auto order = seissolParams.output.freeSurfaceParameters.vtkorder;
    auto& freeSurfaceIntegrator = seissolInstance.freeSurfaceIntegrator();
    auto& meshReader = seissolInstance.meshReader();
    io::writer::ScheduledWriter schedWriter;
    schedWriter.name = "free-surface";
    schedWriter.interval = seissolParams.output.freeSurfaceParameters.interval;
    auto* surfaceMeshIds =
        freeSurfaceIntegrator.surfaceLtsTree.var(freeSurfaceIntegrator.surfaceLts.meshId);
    auto* surfaceMeshSides =
        freeSurfaceIntegrator.surfaceLtsTree.var(freeSurfaceIntegrator.surfaceLts.side);
    auto writer = io::instance::mesh::VtkHdfWriter(
        "free-surface", freeSurfaceIntegrator.surfaceLtsTree.getNumberOfCells(), 2, order);
    writer.addPointProjector([=](double* target, std::size_t index) {
      auto meshId = surfaceMeshIds[index];
      auto side = surfaceMeshSides[index];
      const auto& element = meshReader.getElements()[meshId];
      const auto& vertexArray = meshReader.getVertices();

      // for the very time being, circumvent the bounding box mechanism of Yateto as follows.
      const double zero[2] = {0, 0};
      double xez[3];
      seissol::transformations::chiTau2XiEtaZeta(side, zero, xez);
      seissol::transformations::tetrahedronReferenceToGlobal(
          vertexArray[element.vertices[0]].coords,
          vertexArray[element.vertices[1]].coords,
          vertexArray[element.vertices[2]].coords,
          vertexArray[element.vertices[3]].coords,
          xez,
          &target[0]);
      for (std::size_t i = 1; i < tensor::vtk2d::Shape[order][1]; ++i) {
        double point[2] = {init::vtk2d::Values[order][i * 2 - 2 + 0],
                           init::vtk2d::Values[order][i * 2 - 2 + 1]};
        seissol::transformations::chiTau2XiEtaZeta(side, point, xez);
        seissol::transformations::tetrahedronReferenceToGlobal(
            vertexArray[element.vertices[0]].coords,
            vertexArray[element.vertices[1]].coords,
            vertexArray[element.vertices[2]].coords,
            vertexArray[element.vertices[3]].coords,
            xez,
            &target[i * 3]);
      }
    });
    std::vector<std::string> quantityLabels = {"v1", "v2", "v3", "u1", "u2", "u3"};
    for (std::size_t quantity = 0; quantity < FREESURFACE_NUMBER_OF_COMPONENTS; ++quantity) {
      writer.addPointData<real>(quantityLabels[quantity], {}, [=](real* target, std::size_t index) {
        auto meshId = surfaceMeshIds[index];
        auto side = surfaceMeshSides[index];
        const auto* dofsAllQuantities = ltsLut->lookup(lts->dofs, meshId);
        const auto* dofsSingleQuantity =
            dofsAllQuantities + QDofSizePadded * (6 + quantity); // velocities
        kernel::projectBasisToVtkFaceFromVolume vtkproj;
        vtkproj.qb = dofsSingleQuantity;
        vtkproj.xf(order) = target;
        vtkproj.collvf(ConvergenceOrder, order, side) =
            init::collvf::Values[ConvergenceOrder + (ConvergenceOrder + 1) * (order + 9 * side)];
        vtkproj.execute(order, side);
      });
    }
    for (std::size_t quantity = 0; quantity < FREESURFACE_NUMBER_OF_COMPONENTS; ++quantity) {
      writer.addPointData<real>(
          quantityLabels[quantity + FREESURFACE_NUMBER_OF_COMPONENTS],
          {},
          [=](real* target, std::size_t index) {
            auto meshId = surfaceMeshIds[index];
            auto side = surfaceMeshSides[index];
            const auto* faceDisplacements = ltsLut->lookup(lts->faceDisplacements, meshId);
            const auto* faceDisplacementVariable =
                faceDisplacements[side] + FaceDisplacementPadded * quantity;
            kernel::projectNodalToVtkFace vtkproj;
            vtkproj.pn = faceDisplacementVariable;
            vtkproj.MV2nTo2m = nodal::init::MV2nTo2m::Values;
            vtkproj.xf(order) = target;
            vtkproj.collff(ConvergenceOrder, order) =
                init::collff::Values[ConvergenceOrder + (ConvergenceOrder + 1) * order];
            vtkproj.execute(order);
          });
    }
    schedWriter.planWrite = writer.makeWriter();
    seissolInstance.getOutputManager().addOutput(schedWriter);
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

  seissolInstance.flopCounter().init(seissolParams.output.prefix);

  seissolInstance.analysisWriter().init(&seissolInstance.meshReader(), seissolParams.output.prefix);
}

void initFaultOutputManager(seissol::SeisSol& seissolInstance) {
  const auto& backupTimeStamp = seissolInstance.getBackupTimeStamp();
  seissolInstance.getMemoryManager().initFaultOutputManager(backupTimeStamp);
  auto* faultOutputManager = seissolInstance.getMemoryManager().getFaultOutputManager();
  seissolInstance.timeManager().setFaultOutputManager(faultOutputManager);
}

void enableWaveFieldOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  if (seissolParams.output.waveFieldParameters.enabled &&
      seissolParams.output.waveFieldParameters.vtkorder < 0) {
    seissolInstance.waveFieldWriter().enable();
    seissolInstance.waveFieldWriter().setFilename(seissolParams.output.prefix.c_str());
    seissolInstance.waveFieldWriter().setWaveFieldInterval(
        seissolParams.output.waveFieldParameters.interval);
  }
}

void enableFreeSurfaceOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
  if (seissolParams.output.freeSurfaceParameters.enabled) {
    int refinement = seissolParams.output.freeSurfaceParameters.refinement;
    if (seissolParams.output.freeSurfaceParameters.vtkorder < 0) {
      seissolInstance.freeSurfaceWriter().enable();
    } else {
      refinement = 0;
    }

    seissolInstance.freeSurfaceIntegrator().initialize(refinement,
                                                       memoryManager.getGlobalDataOnHost(),
                                                       memoryManager.getLts(),
                                                       memoryManager.getLtsTree(),
                                                       memoryManager.getLtsLut());
  }
}

void setIntegralMask(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  seissolInstance.postProcessor().setIntegrationMask(
      seissolParams.output.waveFieldParameters.integrationMask);
}

} // namespace

void seissol::initializer::initprocedure::initIO(seissol::SeisSol& seissolInstance) {
  const auto rank = MPI::mpi.rank();
  logInfo() << "Begin init output.";

  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const filesystem::path outputPath(seissolParams.output.prefix);
  const auto outputDir = filesystem::directory_entry(outputPath.parent_path());
  if (!filesystem::exists(outputDir)) {
    logWarning() << "Output directory does not exist yet. We therefore create it now.";
    if (rank == 0) {
      filesystem::create_directory(outputDir);
    }
  }
  seissol::MPI::barrier(MPI::mpi.comm());

  enableWaveFieldOutput(seissolInstance);
  setIntegralMask(seissolInstance);
  enableFreeSurfaceOutput(seissolInstance);
  initFaultOutputManager(seissolInstance);
  setupCheckpointing(seissolInstance);
  setupOutput(seissolInstance);
  logInfo() << "End init output.";
}
