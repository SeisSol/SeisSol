// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitIO.h"
#include "Common/Filesystem.h"
#include "Equations/Datastructures.h"
#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Writer/Writer.h"
#include "Init.h"
#include "Numerical/Transformation.h"
#include "SeisSol.h"
#include <Common/Constants.h>
#include <Geometry/MeshDefinition.h>
#include <IO/Instance/Geometry/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/MemoryAllocator.h>
#include <Memory/Tree/Layer.h>
#include <Model/Plasticity.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include <Solver/MultipleSimulations.h>
#include <algorithm>
#include <array>
#include <cstdint>
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
    std::vector<std::size_t> globalIds(tree->size(seissol::initializer::LayerMask(Ghost)));
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
    std::vector<std::size_t> faceIdentifiers(tree->size(seissol::initializer::LayerMask(Ghost)));
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

  {
    auto* tree = seissolInstance.getMemoryManager().getSurfaceTree();
    auto* surf = seissolInstance.getMemoryManager().getSurface();
    std::vector<std::size_t> faceIdentifiers(tree->size(seissol::initializer::LayerMask(Ghost)));
    const auto* meshIds = tree->var(surf->meshId);
    const auto* sides = tree->var(surf->side);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
      // same as for DR
      faceIdentifiers[i] = meshIds[i] * 4 + static_cast<std::size_t>(sides[i]);
    }
    checkpoint.registerTree("surface", tree, faceIdentifiers);
    surf->registerCheckpointVariables(checkpoint, tree);
  }

  const auto& checkpointFile = seissolInstance.getCheckpointLoadFile();
  if (checkpointFile.has_value()) {
    const double time = seissolInstance.getOutputManager().loadCheckpoint(checkpointFile.value());
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

  // TODO(David): change Yateto/TensorForge interface to make padded sizes more accessible
  constexpr auto QDofSizePadded =
      tensor::Q::Size / tensor::Q::Shape[multisim::BasisFunctionDimension + 1];
  constexpr auto FaceDisplacementPadded =
      tensor::faceDisplacement::Size /
      tensor::faceDisplacement::Shape[multisim::BasisFunctionDimension + 1];

  const auto namewrap = [](const std::string& name, std::size_t sim) {
    if constexpr (multisim::MultisimEnabled) {
      return name + "-" + std::to_string(sim + 1);
    } else {
      return name;
    }
  };

  if (seissolParams.output.waveFieldParameters.enabled) {
    const auto orderIO = seissolParams.output.waveFieldParameters.vtkorder;
    const auto order = std::max(0, orderIO);
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
    auto writer = io::instance::geometry::GeometryWriter(
        "wavefield", celllist.size(), io::instance::geometry::Shape::Tetrahedron, order);

    writer.addPointProjector([=](double* target, std::size_t index) {
      const auto& element = meshReader.getElements()[cellIndices[index]];
      const auto& vertexArray = meshReader.getVertices();

      const auto trueOrder = order > 0 ? order : 1;

      // for the very time being, circumvent the bounding box mechanism of Yateto as follows.
      const double zero[3] = {0, 0, 0};
      seissol::transformations::tetrahedronReferenceToGlobal(
          vertexArray[element.vertices[0]].coords,
          vertexArray[element.vertices[1]].coords,
          vertexArray[element.vertices[2]].coords,
          vertexArray[element.vertices[3]].coords,
          zero,
          &target[0]);
      for (std::size_t i = 1; i < tensor::vtk3d::Shape[trueOrder][1]; ++i) {
        double point[3] = {init::vtk3d::Values[trueOrder][i * 3 - 3 + 0],
                           init::vtk3d::Values[trueOrder][i * 3 - 3 + 1],
                           init::vtk3d::Values[trueOrder][i * 3 - 3 + 2]};
        seissol::transformations::tetrahedronReferenceToGlobal(
            vertexArray[element.vertices[0]].coords,
            vertexArray[element.vertices[1]].coords,
            vertexArray[element.vertices[2]].coords,
            vertexArray[element.vertices[3]].coords,
            point,
            &target[i * 3]);
      }
    });

    writer.addCellData<uint64_t>("clustering", {}, true, [=](uint64_t* target, std::size_t index) {
      target[0] = meshReader.getElements()[index].clusterId;
    });

    writer.addCellData<std::size_t>(
        "global-id", {}, true, [=](std::size_t* target, std::size_t index) {
          target[0] = meshReader.getElements()[index].globalId;
        });

    for (std::size_t sim = 0; sim < seissol::multisim::NumSimulations; ++sim) {
      for (std::size_t quantity = 0; quantity < seissol::model::MaterialT::Quantities.size();
           ++quantity) {
        if (seissolParams.output.waveFieldParameters.outputMask[quantity]) {
          writer.addGeometryOutput<real>(
              namewrap(seissol::model::MaterialT::Quantities[quantity], sim),
              {},
              false,
              [=](real* target, std::size_t index) {
                const auto* dofsAllQuantities = ltsLut->lookup(lts->dofs, cellIndices[index]);
                const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
                kernel::projectBasisToVtkVolume vtkproj{};
                memory::AlignedArray<real, multisim::NumSimulations> simselect;
                simselect[sim] = 1;
                vtkproj.simselect = simselect.data();
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
            writer.addGeometryOutput<real>(
                namewrap(seissol::model::PlasticityData::Quantities[quantity], sim),
                {},
                false,
                [=](real* target, std::size_t index) {
                  const auto* dofsAllQuantities = ltsLut->lookup(lts->pstrain, cellIndices[index]);
                  const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
                  kernel::projectBasisToVtkVolume vtkproj{};
                  memory::AlignedArray<real, multisim::NumSimulations> simselect;
                  simselect[sim] = 1;
                  vtkproj.simselect = simselect.data();
                  vtkproj.qb = dofsSingleQuantity;
                  vtkproj.xv(order) = target;
                  vtkproj.collvv(ConvergenceOrder, order) =
                      init::collvv::Values[ConvergenceOrder + (ConvergenceOrder + 1) * order];
                  vtkproj.execute(order);
                });
          }
        }
      }
    }
    schedWriter.planWrite = writer.makeWriter();
    seissolInstance.getOutputManager().addOutput(schedWriter);
  }

  if (seissolParams.output.freeSurfaceParameters.enabled) {
    const auto orderIO = seissolParams.output.freeSurfaceParameters.vtkorder;
    const auto order = std::max(0, orderIO);

    auto& freeSurfaceIntegrator = seissolInstance.freeSurfaceIntegrator();
    auto& meshReader = seissolInstance.meshReader();
    io::writer::ScheduledWriter schedWriter;
    schedWriter.name = "free-surface";
    schedWriter.interval = seissolParams.output.freeSurfaceParameters.interval;
    auto* surfaceMeshIds =
        freeSurfaceIntegrator.surfaceLtsTree->var(freeSurfaceIntegrator.surfaceLts->meshId);
    auto* surfaceMeshSides =
        freeSurfaceIntegrator.surfaceLtsTree->var(freeSurfaceIntegrator.surfaceLts->side);
    auto* surfaceLocationFlag =
        freeSurfaceIntegrator.surfaceLtsTree->var(freeSurfaceIntegrator.surfaceLts->locationFlag);
    auto writer =
        io::instance::geometry::GeometryWriter("free-surface",
                                               freeSurfaceIntegrator.surfaceLtsTree->size(),
                                               io::instance::geometry::Shape::Triangle,
                                               order);
    writer.addPointProjector([=, &freeSurfaceIntegrator](double* target, std::size_t index) {
      auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
      auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
      const auto& element = meshReader.getElements()[meshId];
      const auto& vertexArray = meshReader.getVertices();

      const auto trueOrder = order > 0 ? order : 1;

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
      for (std::size_t i = 1; i < tensor::vtk2d::Shape[trueOrder][1]; ++i) {
        double point[2] = {init::vtk2d::Values[trueOrder][i * 2 - 2 + 0],
                           init::vtk2d::Values[trueOrder][i * 2 - 2 + 1]};
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

    writer.addCellData<std::uint8_t>(
        "locationFlag",
        {},
        true,
        [=, &freeSurfaceIntegrator](std::uint8_t* target, std::size_t index) {
          target[0] = surfaceLocationFlag[freeSurfaceIntegrator.backmap[index]];
        });

    writer.addCellData<std::size_t>(
        "global-id", {}, true, [=, &freeSurfaceIntegrator](std::size_t* target, std::size_t index) {
          const auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
          const auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
          target[0] = meshReader.getElements()[meshId].globalId * 4 + side;
        });

    std::vector<std::string> quantityLabels = {"v1", "v2", "v3", "u1", "u2", "u3"};
    for (std::size_t sim = 0; sim < seissol::multisim::NumSimulations; ++sim) {
      for (std::size_t quantity = 0;
           quantity < seissol::solver::FreeSurfaceIntegrator::NumComponents;
           ++quantity) {
        writer.addGeometryOutput<real>(
            namewrap(quantityLabels[quantity], sim),
            {},
            false,
            [=, &freeSurfaceIntegrator](real* target, std::size_t index) {
              auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
              auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
              const auto* dofsAllQuantities = ltsLut->lookup(lts->dofs, meshId);
              const auto* dofsSingleQuantity =
                  dofsAllQuantities + QDofSizePadded * (6 + quantity); // velocities
              kernel::projectBasisToVtkFaceFromVolume vtkproj{};
              memory::AlignedArray<real, multisim::NumSimulations> simselect;
              simselect[sim] = 1;
              vtkproj.simselect = simselect.data();
              vtkproj.qb = dofsSingleQuantity;
              vtkproj.xf(order) = target;
              vtkproj.collvf(ConvergenceOrder, order, side) =
                  init::collvf::Values[ConvergenceOrder +
                                       (ConvergenceOrder + 1) * (order + 9 * side)];
              vtkproj.execute(order, side);
            });
      }
      for (std::size_t quantity = 0;
           quantity < seissol::solver::FreeSurfaceIntegrator::NumComponents;
           ++quantity) {
        writer.addGeometryOutput<real>(
            namewrap(
                quantityLabels[quantity + seissol::solver::FreeSurfaceIntegrator::NumComponents],
                sim),
            {},
            false,
            [=, &freeSurfaceIntegrator](real* target, std::size_t index) {
              auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
              auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
              const auto* faceDisplacements = ltsLut->lookup(lts->faceDisplacements, meshId);
              const auto* faceDisplacementVariable =
                  faceDisplacements[side] + FaceDisplacementPadded * quantity;
              kernel::projectNodalToVtkFace vtkproj{};
              memory::AlignedArray<real, multisim::NumSimulations> simselect;
              simselect[sim] = 1;
              vtkproj.simselect = simselect.data();
              vtkproj.pn = faceDisplacementVariable;
              vtkproj.MV2nTo2m = nodal::init::MV2nTo2m::Values;
              vtkproj.xf(order) = target;
              vtkproj.collff(ConvergenceOrder, order) =
                  init::collff::Values[ConvergenceOrder + (ConvergenceOrder + 1) * order];
              vtkproj.execute(order);
            });
      }
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
    receiverWriter.addPoints(
        seissolInstance.meshReader(), *ltsLut, *lts, memoryManager.getGlobalData());
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

void enableFreeSurfaceOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
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

  setIntegralMask(seissolInstance);
  enableFreeSurfaceOutput(seissolInstance);
  initFaultOutputManager(seissolInstance);
  setupCheckpointing(seissolInstance);
  setupOutput(seissolInstance);
  logInfo() << "End init output.";
}
