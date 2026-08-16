// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitIO.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Filesystem.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Points.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Writer/Writer.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "Numerical/Projection.h"
#include "Numerical/Transformation.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utils/logger.h>
#include <vector>

namespace {

using namespace seissol;

namespace projection = seissol::numerical::projection;

// The projection matrices are generated for every convergence order up to the compiled one, so
// that a per-cell order (cf. #1421) only requires selecting a different entry at run time.
constexpr std::size_t MinProjectionOrder = 1;
constexpr std::size_t MaxProjectionOrder = ConvergenceOrder;

/**
 * The padded leading dimension of a generated projection tensor, i.e. the stride between two
 * consecutive basis functions. Yateto stores these matrices as [point][basisFunction] with an
 * aligned stride on the point dimension; we read the padding back off the generated metadata
 * instead of re-deriving it from the alignment.
 */
template <typename TensorT>
std::size_t projectionStride(std::size_t degree) {
  const auto index = TensorT::index(ConvergenceOrder, degree);
  return TensorT::Size[index] / TensorT::Shape[index][1];
}

//! The affine embedding of the reference triangle into the given side of the reference tetrahedron.
seissol::numerical::AffineMap<2, 3> faceEmbedding(std::size_t side) {
  const std::array<std::array<double, 2>, 3> corners = {
      std::array<double, 2>{0, 0}, std::array<double, 2>{1, 0}, std::array<double, 2>{0, 1}};
  std::vector<std::array<double, 3>> vertices;
  vertices.reserve(corners.size());
  for (const auto& chiTau : corners) {
    std::array<double, 3> xez{};
    seissol::transformations::chiTau2XiEtaZeta(
        static_cast<std::uint32_t>(side), chiTau.data(), xez.data());
    vertices.push_back(xez);
  }
  return seissol::numerical::AffineMap<2, 3>::fromVertices(vertices);
}

void setupCheckpointing(seissol::SeisSol& seissolInstance) {
  auto& checkpoint = seissolInstance.getOutputManager().getCheckpointManager();

  {
    auto& storage = seissolInstance.getMemoryManager().getLtsStorage();
    std::vector<std::size_t> globalIds(storage.size(seissol::initializer::LayerMask(Ghost)));
    std::size_t offset = 0;
    for (const auto& layer : storage.leaves(Ghost)) {

#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < layer.size(); ++i) {
        const auto meshId = layer.var<LTS::SecondaryInformation>()[i].meshId;
        globalIds[offset + i] = seissolInstance.meshReader().getElements()[meshId].globalId;
      }
      offset += layer.size();
    }
    checkpoint.registerTree("lts", storage, globalIds);
    LTS::registerCheckpointVariables(checkpoint, storage);
  }

  {
    auto& storage = seissolInstance.getMemoryManager().getDRStorage();
    auto& dynrup = seissolInstance.getMemoryManager().getDynamicRupture();
    std::vector<std::size_t> faceIdentifiers(storage.size(seissol::initializer::LayerMask(Ghost)));
    const auto* drFaceInformation = storage.var<DynamicRupture::FaceInformation>();

#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
      auto faultFace = drFaceInformation[i].meshFace;
      const auto& fault = seissolInstance.meshReader().getFault()[faultFace];
      // take the positive cell and side as fault face identifier
      // (should result in roughly twice as large numbers as when indexing all faces; cf.
      // handshake theorem)
      faceIdentifiers[i] = fault.globalId * 4 + fault.side;
    }
    checkpoint.registerTree("dynrup", storage, faceIdentifiers);
    dynrup.registerCheckpointVariables(checkpoint, storage);
  }

  {
    auto& storage = seissolInstance.getMemoryManager().getSurfaceStorage();
    std::vector<std::size_t> faceIdentifiers(storage.size(seissol::initializer::LayerMask(Ghost)));
    const auto* meshIds = storage.var<SurfaceLTS::MeshId>();
    const auto* sides = storage.var<SurfaceLTS::Side>();

#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
      // same as for DR
      faceIdentifiers[i] = meshIds[i] * 4 + static_cast<std::size_t>(sides[i]);
    }
    checkpoint.registerTree("surface", storage, faceIdentifiers);
    SurfaceLTS::registerCheckpointVariables(checkpoint, storage);
  }

  const auto& checkpointFile = seissolInstance.getCheckpointLoadFile();
  if (checkpointFile.has_value()) {
    const double time = seissolInstance.getOutputManager().loadCheckpoint(checkpointFile.value());
    seissolInstance.simulator().setCurrentTime(time);
  }

  if (seissolInstance.getSeisSolParameters().output.checkpointParameters.enabled) {
    // FIXME: for now, we allow only _one_ checkpoint interval which checkpoints everything
    // existent
    seissolInstance.getOutputManager().setupCheckpoint(
        seissolInstance.getSeisSolParameters().output.checkpointParameters.interval);
  }
}

void setupOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  auto& memoryManager = seissolInstance.getMemoryManager();
  auto& ltsStorage = memoryManager.getLtsStorage();
  auto& backmap = memoryManager.getBackmap();
  auto& drStorage = memoryManager.getDRStorage();
  auto* globalData = memoryManager.getGlobalData().onHost;

  // TODO(David): change Yateto/TensorForge interface to make padded sizes more accessible
  constexpr auto QDofSizePadded =
      tensor::Q::Size / tensor::Q::Shape[multisim::BasisFunctionDimension + 1];
  constexpr auto QDofPointsPadded =
      tensor::QStress::Size / tensor::QStress::Shape[multisim::BasisFunctionDimension + 1];
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
    const auto order = static_cast<uint32_t>(std::max(0, orderIO));
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

    const auto dataOrder = order > 0 ? order : 0;
    const auto trueOrder = order > 0 ? order : 1;
    const auto trueBase = io::instance::geometry::pointsTetrahedron(trueOrder);
    const auto dataBase = io::instance::geometry::pointsTetrahedron(dataOrder);

    auto subcells = io::instance::geometry::unrefined<3>();

    if (seissolParams.output.waveFieldParameters.refinement ==
        seissol::initializer::parameters::VolumeRefinement::Refine4) {
      subcells = io::instance::geometry::subdivideMaps(subcells,
                                                       io::instance::geometry::TetrahedronRefine4);
    }
    if (seissolParams.output.waveFieldParameters.refinement ==
        seissol::initializer::parameters::VolumeRefinement::Refine8) {
      subcells = io::instance::geometry::subdivideMaps(subcells,
                                                       io::instance::geometry::TetrahedronRefine8);
    }
    if (seissolParams.output.waveFieldParameters.refinement ==
        seissol::initializer::parameters::VolumeRefinement::Refine32) {
      // the edge division has to come first; the legacy refinement::DivideTetrahedronBy32
      // subdivided by 8 and then split each of those subcells by its center point, which (as
      // subdivideMaps enumerates input-major) is the ordering 4*i + j the output cells had
      subcells = io::instance::geometry::subdivideMaps(subcells,
                                                       io::instance::geometry::TetrahedronRefine8);
      subcells = io::instance::geometry::subdivideMaps(subcells,
                                                       io::instance::geometry::TetrahedronRefine4);
    }

    const auto truePoints = io::instance::geometry::applyMaps(subcells, trueBase);

    const auto projectionTarget =
        seissolParams.output.projection == seissol::initializer::parameters::ProjectionMethod::L2
            ? projection::Target::Project
            : projection::Target::Interpolate;

    const auto makeVolumeTable = [&](projection::Source source,
                                     std::optional<std::size_t> derivative) {
      projection::Spec spec;
      spec.source = source;
      spec.target = projectionTarget;
      spec.derivative = derivative;
      const auto stride = source == projection::Source::Nodal
                              ? projectionStride<tensor::collnv>(order)
                              : projectionStride<tensor::collvv>(order);
      return std::make_shared<projection::Table<3, 3>>(
          subcells, dataBase, dataOrder, stride, spec, MinProjectionOrder, MaxProjectionOrder);
    };

    const auto proj = makeVolumeTable(projection::Source::Modal, {});

    std::array<std::shared_ptr<projection::Table<3, 3>>, Cell::Dim> projD{};
    if (seissolParams.output.waveFieldParameters.computeStrain ||
        seissolParams.output.waveFieldParameters.computeRotation) {
      for (std::size_t direction = 0; direction < Cell::Dim; ++direction) {
        projD[direction] = makeVolumeTable(projection::Source::Modal, direction);
      }
    }

    std::shared_ptr<projection::Table<3, 3>> projNodal;
    if (seissolParams.model.plasticity) {
      projNodal = makeVolumeTable(projection::Source::Nodal, {});
    }

    const auto config = io::instance::geometry::WriterConfig{
        order,
        orderIO < 0 ? io::instance::geometry::WriterFormat::Xdmf
                    : io::instance::geometry::WriterFormat::Vtk,
        seissolParams.output.xdmfWriterBackend ==
                seissol::initializer::parameters::XdmfBackend::Posix
            ? io::instance::geometry::WriterBackend::Binary
            : io::instance::geometry::WriterBackend::Hdf5,
        io::instance::geometry::WriterGroup::FullSnapshot,
        seissolParams.output.hdfcompress};

    io::writer::ScheduledWriter schedWriter;
    schedWriter.name = "wavefield";
    schedWriter.interval = seissolParams.output.waveFieldParameters.interval;
    auto writer = io::instance::geometry::GeometryWriter("wavefield",
                                                         celllist.size(),
                                                         io::instance::geometry::Shape::Tetrahedron,
                                                         config,
                                                         subcells.size());

    writer.addPointProjector([=](double* target, std::size_t index, std::size_t subcell) {
      const auto& element = meshReader.getElements()[cellIndices[index]];
      const auto& vertexArray = meshReader.getVertices();

      for (std::size_t i = 0; i < truePoints[subcell].size(); ++i) {
        seissol::transformations::tetrahedronReferenceToGlobal(
            vertexArray[element.vertices[0]].coords,
            vertexArray[element.vertices[1]].coords,
            vertexArray[element.vertices[2]].coords,
            vertexArray[element.vertices[3]].coords,
            truePoints[subcell][i].data(),
            &target[i * 3]);
      }
    });

    const auto rank = seissol::Mpi::mpi.rank();
    writer.addCellData<int>(
        "partition", {}, true, [=](int* target, std::size_t /*index*/, std::size_t /*subcell*/) {
          target[0] = rank;
        });

    writer.addCellData<uint64_t>(
        "clustering", {}, true, [=](uint64_t* target, std::size_t index, std::size_t /*subcell*/) {
          target[0] = meshReader.getElements()[index].clusterId;
        });

    writer.addCellData<std::size_t>(
        "global-id",
        {},
        true,
        [=](std::size_t* target, std::size_t index, std::size_t /*subcell*/) {
          target[0] = meshReader.getElements()[index].globalId;
        });

    constexpr std::size_t MaxVtk3dPoints =
        tensor::vtk3d::Shape[(sizeof(tensor::vtk3d::Shape) / sizeof(tensor::vtk3d::Shape[0])) - 1]
                            [1];

    for (std::size_t sim = 0; sim < seissol::multisim::NumSimulations; ++sim) {
      const auto projectVolume =
          [=](real* target, const real* dofsSingleQuantity, const real* collvv) {
            kernel::projectBasisToVtkVolume vtkproj{};
            memory::AlignedArray<real, multisim::NumSimulations> simselect{};
            alignas(Alignment) std::array<real, MaxVtk3dPoints> alignedTarget{};
            simselect[sim] = 1;
            vtkproj.simselect = simselect.data();
            vtkproj.qb = dofsSingleQuantity;
            vtkproj.xv(order) = alignedTarget.data();
            vtkproj.collvv(ConvergenceOrder, order) = collvv;
            vtkproj.execute(order);
            std::copy_n(alignedTarget.data(), dataBase.size(), target);
          };

      const auto projectVolumeDeriv = [=](real* target,
                                          const real* dofsSingleQuantity,
                                          std::size_t dir,
                                          std::size_t index,
                                          std::size_t subcell) {
        std::array<real, MaxVtk3dPoints> dataX{};
        std::array<real, MaxVtk3dPoints> dataY{};
        std::array<real, MaxVtk3dPoints> dataZ{};

        projectVolume(dataX.data(), dofsSingleQuantity, (*projD[0])(subcell, ConvergenceOrder));
        projectVolume(dataY.data(), dofsSingleQuantity, (*projD[1])(subcell, ConvergenceOrder));
        projectVolume(dataZ.data(), dofsSingleQuantity, (*projD[2])(subcell, ConvergenceOrder));

        const auto& element = meshReader.getElements()[cellIndices[index]];
        const auto& vertexArray = meshReader.getVertices();

        std::array<double, Cell::NumVertices> coordsX{};
        std::array<double, Cell::NumVertices> coordsY{};
        std::array<double, Cell::NumVertices> coordsZ{};
        std::array<double, Cell::Dim> gradXi{};
        std::array<double, Cell::Dim> gradEta{};
        std::array<double, Cell::Dim> gradZeta{};

        for (std::size_t i = 0; i < Cell::NumVertices; ++i) {
          coordsX[i] = vertexArray[element.vertices[i]].coords[0];
          coordsY[i] = vertexArray[element.vertices[i]].coords[1];
          coordsZ[i] = vertexArray[element.vertices[i]].coords[2];
        }

        seissol::transformations::tetrahedronGlobalToReferenceJacobian(coordsX.data(),
                                                                       coordsY.data(),
                                                                       coordsZ.data(),
                                                                       gradXi.data(),
                                                                       gradEta.data(),
                                                                       gradZeta.data());

        for (std::size_t i = 0; i < dataBase.size(); ++i) {
          target[i] = dataX[i] * gradXi[dir] + dataY[i] * gradEta[dir] + dataZ[i] * gradZeta[dir];
        }
      };

      for (std::size_t quantity = 0; quantity < seissol::model::MaterialT::Quantities.size();
           ++quantity) {

        if (seissolParams.output.waveFieldParameters.outputMask[quantity]) {
          writer.addGeometryOutput<real>(
              namewrap(seissol::model::MaterialT::Quantities[quantity], sim),
              {},
              false,
              [=, &ltsStorage, &backmap](real* target, std::size_t index, std::size_t subcell) {
                const auto position = backmap.get(cellIndices[index]);
                const auto* dofsAllQuantities = ltsStorage.lookup<LTS::Dofs>(position);
                const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
                projectVolume(target, dofsSingleQuantity, (*proj)(subcell, ConvergenceOrder));
              });
        }

        if (seissolParams.output.waveFieldParameters.integrationMask[quantity]) {
          writer.addGeometryOutput<real>(
              namewrap("int-" + seissol::model::MaterialT::Quantities[quantity], sim),
              {},
              false,
              [=, &ltsStorage, &backmap](real* target, std::size_t index, std::size_t subcell) {
                const auto position = backmap.get(cellIndices[index]);
                const auto* dofsAllQuantities = ltsStorage.lookup<LTS::Integrals>(position);
                const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
                projectVolume(target, dofsSingleQuantity, (*proj)(subcell, ConvergenceOrder));
              });
        }
      }

      using Idx = std::tuple<std::string, std::size_t, std::size_t>;

      if (seissolParams.output.waveFieldParameters.computeStrain) {
        for (const auto& idxp : {Idx{"xx", 0, 0},
                                 Idx{"yy", 1, 1},
                                 Idx{"zz", 2, 2},
                                 Idx{"xy", 0, 1},
                                 Idx{"yz", 1, 2},
                                 Idx{"xz", 0, 2}}) {
          const auto& name = std::get<0>(idxp);

          // need non-reference capture (due to lambda usage)
          const auto idx1 = std::get<1>(idxp);
          const auto idx2 = std::get<2>(idxp);

          // compute (d_i1 v_i2 + d_i2 v_i1) / 2

          writer.addGeometryOutput<real>(
              namewrap("eps" + name, sim),
              {},
              false,
              [=, &ltsStorage, &backmap](real* target, std::size_t index, std::size_t subcell) {
                const auto position = backmap.get(cellIndices[index]);
                const auto* dofsAllQuantities = ltsStorage.lookup<LTS::Integrals>(position);
                const auto* dofsSingleQuantity1 =
                    dofsAllQuantities +
                    QDofSizePadded * (idx1 + model::MaterialT::TractionQuantities);
                projectVolumeDeriv(target, dofsSingleQuantity1, idx2, index, subcell);

                if (idx1 != idx2) {
                  const auto* dofsSingleQuantity2 =
                      dofsAllQuantities +
                      QDofSizePadded * (idx2 + model::MaterialT::TractionQuantities);
                  std::array<real, MaxVtk3dPoints> itarget{};
                  projectVolumeDeriv(itarget.data(), dofsSingleQuantity2, idx1, index, subcell);

                  for (std::size_t i = 0; i < dataBase.size(); ++i) {
                    target[i] = (target[i] + itarget[i]) / 2;
                  }
                }
              });
        }
      }

      if (seissolParams.output.waveFieldParameters.computeRotation) {
        for (const auto& idxp : {Idx{"1", 2, 1}, Idx{"2", 0, 2}, Idx{"3", 1, 0}}) {
          const auto& name = std::get<0>(idxp);

          // need non-reference capture (due to lambda usage)
          const auto idx1 = std::get<1>(idxp);
          const auto idx2 = std::get<2>(idxp);

          // compute d_i2 v_i1 - d_i1 v_i2

          writer.addGeometryOutput<real>(
              namewrap("rot" + name, sim),
              {},
              false,
              [=, &ltsStorage, &backmap](real* target, std::size_t index, std::size_t subcell) {
                const auto position = backmap.get(cellIndices[index]);
                const auto* dofsAllQuantities = ltsStorage.lookup<LTS::Dofs>(position);
                const auto* dofsSingleQuantity1 =
                    dofsAllQuantities +
                    QDofSizePadded * (idx1 + model::MaterialT::TractionQuantities);
                projectVolumeDeriv(target, dofsSingleQuantity1, idx2, index, subcell);

                const auto* dofsSingleQuantity2 =
                    dofsAllQuantities +
                    QDofSizePadded * (idx2 + model::MaterialT::TractionQuantities);
                std::array<real, MaxVtk3dPoints> itarget{};
                projectVolumeDeriv(itarget.data(), dofsSingleQuantity2, idx1, index, subcell);

                for (std::size_t i = 0; i < tensor::vtk3d::Shape[order][1]; ++i) {
                  target[i] -= itarget[i];
                }
              });
        }
      }

      if (seissolParams.model.plasticity) {
        for (std::size_t quantity = 0; quantity < seissol::model::PlasticityData::Quantities.size();
             ++quantity) {
          if (seissolParams.output.waveFieldParameters.plasticityMask[quantity]) {
            constexpr std::size_t MaxVtk3dPoints = tensor::vtk3d::Shape
                [(sizeof(tensor::vtk3d::Shape) / sizeof(tensor::vtk3d::Shape[0])) - 1][1];
            writer.addGeometryOutput<real>(
                namewrap(seissol::model::PlasticityData::Quantities[quantity], sim),
                {},
                false,
                [=, &ltsStorage, &backmap](real* target, std::size_t index, std::size_t subcell) {
                  const auto position = backmap.get(cellIndices[index]);
                  const auto* dofsAllQuantities = ltsStorage.lookup<LTS::PStrain>(position);
                  const auto* pointsSingleQuantity =
                      dofsAllQuantities + QDofPointsPadded * quantity;
                  kernel::projectNodalToVtkVolume vtkproj{};
                  memory::AlignedArray<real, multisim::NumSimulations> simselect{};
                  alignas(Alignment) std::array<real, MaxVtk3dPoints> alignedTarget{};
                  simselect[sim] = 1;
                  vtkproj.simselect = simselect.data();
                  vtkproj.qn = pointsSingleQuantity;
                  vtkproj.xv(order) = alignedTarget.data();
                  vtkproj.collnv(ConvergenceOrder, order) = (*projNodal)(subcell, ConvergenceOrder);
                  vtkproj.execute(order);
                  std::copy_n(alignedTarget.data(), dataBase.size(), target);
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
    const auto order = static_cast<std::uint32_t>(std::max(0, orderIO));

    auto& freeSurfaceIntegrator = seissolInstance.freeSurfaceIntegrator();
    auto& meshReader = seissolInstance.meshReader();
    io::writer::ScheduledWriter schedWriter;
    schedWriter.name = "free-surface";
    schedWriter.interval = seissolParams.output.freeSurfaceParameters.interval;
    auto* surfaceMeshIds = freeSurfaceIntegrator.surfaceStorage->var<SurfaceLTS::MeshId>();
    auto* surfaceMeshSides = freeSurfaceIntegrator.surfaceStorage->var<SurfaceLTS::Side>();
    auto* surfaceLocationFlag =
        freeSurfaceIntegrator.surfaceStorage->var<SurfaceLTS::LocationFlag>();

    const auto trueOrder = order > 0 ? order : 1;
    const auto dataOrder = order > 0 ? order : 0;
    const auto trueBase = io::instance::geometry::pointsTriangle(trueOrder);
    const auto dataBase = io::instance::geometry::pointsTriangle(dataOrder);

    auto subcells = io::instance::geometry::unrefined<2>();

    for (std::size_t i = 0; i < seissolParams.output.freeSurfaceParameters.refinement; ++i) {
      subcells =
          io::instance::geometry::subdivideMaps(subcells, io::instance::geometry::TriangleRefine4);
    }

    const auto truePoints = io::instance::geometry::applyMaps(subcells, trueBase);

    const auto config = io::instance::geometry::WriterConfig{
        order,
        orderIO < 0 ? io::instance::geometry::WriterFormat::Xdmf
                    : io::instance::geometry::WriterFormat::Vtk,
        seissolParams.output.xdmfWriterBackend ==
                seissol::initializer::parameters::XdmfBackend::Posix
            ? io::instance::geometry::WriterBackend::Binary
            : io::instance::geometry::WriterBackend::Hdf5,
        io::instance::geometry::WriterGroup::FullSnapshot,
        seissolParams.output.hdfcompress};

    auto writer = io::instance::geometry::GeometryWriter("free-surface",
                                                         freeSurfaceIntegrator.backmap.size(),
                                                         io::instance::geometry::Shape::Triangle,
                                                         config,
                                                         subcells.size());

    writer.addPointProjector(
        [=, &freeSurfaceIntegrator](double* target, std::size_t index, std::size_t subcell) {
          auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
          auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
          const auto& element = meshReader.getElements()[meshId];
          const auto& vertexArray = meshReader.getVertices();

          double xez[3]{};
          for (std::size_t i = 0; i < truePoints[subcell].size(); ++i) {
            seissol::transformations::chiTau2XiEtaZeta(side, truePoints[subcell][i].data(), xez);
            seissol::transformations::tetrahedronReferenceToGlobal(
                vertexArray[element.vertices[0]].coords,
                vertexArray[element.vertices[1]].coords,
                vertexArray[element.vertices[2]].coords,
                vertexArray[element.vertices[3]].coords,
                xez,
                &target[i * 3]);
          }
        });

    const auto projectionTarget =
        seissolParams.output.projection == seissol::initializer::parameters::ProjectionMethod::L2
            ? projection::Target::Project
            : projection::Target::Interpolate;

    // volume basis -> face points, one table per side of the reference tetrahedron
    std::array<std::shared_ptr<projection::Table<2, 3>>, Cell::NumFaces> proj{};
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto embedding = faceEmbedding(f);
      std::vector<seissol::numerical::AffineMap<2, 3>> embedded;
      embedded.reserve(subcells.size());
      for (const auto& subcell : subcells) {
        embedded.emplace_back(embedding.compose(subcell));
      }

      projection::Spec spec;
      spec.target = projectionTarget;
      proj[f] = std::make_shared<projection::Table<2, 3>>(embedded,
                                                          dataBase,
                                                          dataOrder,
                                                          projectionStride<tensor::collvf>(order),
                                                          spec,
                                                          MinProjectionOrder,
                                                          MaxProjectionOrder);
    }

    // face nodes -> face points (the nodal-to-modal transform is folded in)
    projection::Spec faceSpec;
    faceSpec.source = projection::Source::Nodal;
    faceSpec.target = projectionTarget;
    const auto projf =
        std::make_shared<projection::Table<2, 2>>(subcells,
                                                  dataBase,
                                                  dataOrder,
                                                  projectionStride<tensor::collnf>(order),
                                                  faceSpec,
                                                  MinProjectionOrder,
                                                  MaxProjectionOrder);

    const auto rank = seissol::Mpi::mpi.rank();
    writer.addCellData<int>(
        "partition", {}, true, [=](int* target, std::size_t /*index*/, std::size_t /*subcell*/) {
          target[0] = rank;
        });

    writer.addCellData<std::uint8_t>(
        "locationFlag",
        {},
        true,
        [=,
         &freeSurfaceIntegrator](std::uint8_t* target, std::size_t index, std::size_t /*subcell*/) {
          target[0] = surfaceLocationFlag[freeSurfaceIntegrator.backmap[index]];
        });

    writer.addCellData<std::size_t>(
        "global-id",
        {},
        true,
        [=,
         &freeSurfaceIntegrator](std::size_t* target, std::size_t index, std::size_t /*subcell*/) {
          const auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
          const auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
          target[0] = meshReader.getElements()[meshId].globalId * 4 + side;
        });

    const std::vector<std::string> quantityLabelsVelocities = {"v1", "v2", "v3"};
    const std::vector<std::string> quantityLabelsDisplacement = {"u1", "u2", "u3"};

    for (std::size_t sim = 0; sim < seissol::multisim::NumSimulations; ++sim) {
      for (std::size_t quantity = 0; quantity < quantityLabelsVelocities.size(); ++quantity) {
        constexpr std::size_t MaxVtk2dPoints =
            tensor::vtk2d::Shape[(sizeof(tensor::vtk2d::Shape) / sizeof(tensor::vtk2d::Shape[0])) -
                                 1][1];
        writer.addGeometryOutput<real>(
            namewrap(quantityLabelsVelocities[quantity], sim),
            {},
            false,
            [=, &freeSurfaceIntegrator, &ltsStorage, &backmap](
                real* target, std::size_t index, std::size_t subcell) {
              auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
              auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
              const auto position = backmap.get(meshId);
              const auto* dofsAllQuantities = ltsStorage.lookup<LTS::Dofs>(position);

              // velocities start at model::MaterialT::TractionQuantities
              const auto* dofsSingleQuantity =
                  dofsAllQuantities +
                  QDofSizePadded * (model::MaterialT::TractionQuantities + quantity);
              kernel::projectBasisToVtkFaceFromVolume vtkproj{};
              memory::AlignedArray<real, multisim::NumSimulations> simselect{};
              alignas(Alignment) std::array<real, MaxVtk2dPoints> alignedTarget{};
              simselect[sim] = 1;
              vtkproj.simselect = simselect.data();
              vtkproj.qb = dofsSingleQuantity;
              vtkproj.xf(order) = alignedTarget.data();
              vtkproj.collvf(ConvergenceOrder, order) = (*proj[side])(subcell, ConvergenceOrder);
              vtkproj.execute(order);
              std::copy_n(alignedTarget.data(), dataBase.size(), target);
            });
      }
      for (std::size_t quantity = 0; quantity < quantityLabelsDisplacement.size(); ++quantity) {
        constexpr std::size_t MaxVtk2dPoints =
            tensor::vtk2d::Shape[(sizeof(tensor::vtk2d::Shape) / sizeof(tensor::vtk2d::Shape[0])) -
                                 1][1];
        writer.addGeometryOutput<real>(
            namewrap(quantityLabelsDisplacement[quantity], sim),
            {},
            false,
            [=, &freeSurfaceIntegrator, &ltsStorage, &backmap](
                real* target, std::size_t index, std::size_t subcell) {
              auto meshId = surfaceMeshIds[freeSurfaceIntegrator.backmap[index]];
              auto side = surfaceMeshSides[freeSurfaceIntegrator.backmap[index]];
              const auto position = backmap.get(meshId);
              const auto& faceDisplacements = ltsStorage.lookup<LTS::FaceDisplacements>(position);
              const auto* faceDisplacementVariable =
                  faceDisplacements[side] + FaceDisplacementPadded * quantity;
              kernel::projectNodalToVtkFace vtkproj{};
              memory::AlignedArray<real, multisim::NumSimulations> simselect{};
              alignas(Alignment) std::array<real, MaxVtk2dPoints> alignedTarget{};
              simselect[sim] = 1;
              vtkproj.simselect = simselect.data();
              vtkproj.pn = faceDisplacementVariable;
              vtkproj.xf(order) = alignedTarget.data();
              vtkproj.collnf(ConvergenceOrder, order) = (*projf)(subcell, ConvergenceOrder);
              vtkproj.execute(order);
              std::copy_n(alignedTarget.data(), dataBase.size(), target);
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
    receiverWriter.addPoints(seissolInstance.meshReader(), backmap, memoryManager.getGlobalData());
    seissolInstance.timeManager().setReceiverClusters(receiverWriter);
  }

  if (seissolParams.output.energyParameters.enabled) {
    auto& energyOutput = seissolInstance.energyOutput();

    energyOutput.init(globalData,
                      drStorage,
                      seissolInstance.meshReader(),
                      ltsStorage,
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

void enableFreeSurfaceOutput(seissol::SeisSol& seissolInstance) {}

} // namespace

void seissol::initializer::initprocedure::initIO(seissol::SeisSol& seissolInstance) {
  const auto rank = Mpi::mpi.rank();
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
  seissol::Mpi::barrier(Mpi::mpi.comm());

  enableFreeSurfaceOutput(seissolInstance);
  initFaultOutputManager(seissolInstance);
  setupCheckpointing(seissolInstance);
  setupOutput(seissolInstance);
  logInfo() << "End init output.";
}
