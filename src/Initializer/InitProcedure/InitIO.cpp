// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitIO.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Constants.h"
#include "Common/Filesystem.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Expr/Program.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshDefinition.h"
#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Points.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Writer/Writer.h"
#include "Initializer/InitProcedure/DerivedOutput.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "Numerical/Projection.h"
#include "Parallel/MPI.h"
#include "ResultWriter/MiniSeisSolWriter.h"
#include "SeisSol.h"
#include "Solver/FreeSurfaceIntegrator.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace {
using namespace seissol;

//! @brief The file grouping a time series mode asks the writer for.
io::instance::geometry::WriterGroup
    writerGroupOf(seissol::initializer::parameters::TimeSeriesMode mode) {
  switch (mode) {
  case seissol::initializer::parameters::TimeSeriesMode::Incremental:
    return io::instance::geometry::WriterGroup::IncrementalSnapshot;
  case seissol::initializer::parameters::TimeSeriesMode::Monolith:
    return io::instance::geometry::WriterGroup::Monolith;
  case seissol::initializer::parameters::TimeSeriesMode::Snapshot:
    return io::instance::geometry::WriterGroup::FullSnapshot;
  }
  return io::instance::geometry::WriterGroup::FullSnapshot;
}

namespace projection = seissol::numerical::projection;

// Which nodal points the plastic strain lives on is a build option (PLASTICITY_METHOD): "nb"
// uses a unisolvent warp&blend set, "ip" the conical-product quadrature points. Read that back
// off the generated matrices instead of duplicating the CMake variable.
template <typename Cfg>
constexpr auto PlasticityNodalSet = static_cast<std::size_t>(tensor::vNodes<Cfg>::Shape[0]) ==
                                            projection::modalSize(Cell::Dim, Cfg::ConvergenceOrder)
                                        ? projection::NodalSet::WarpBlend
                                        : projection::NodalSet::Stroud;

#define SEISSOL_CHECK_NODAL_SET(Cfg)                                                               \
  static_assert(                                                                                   \
      projection::nodalSize(Cell::Dim, Cfg::ConvergenceOrder, PlasticityNodalSet<Cfg>) ==          \
          static_cast<std::size_t>(tensor::vNodes<Cfg>::Shape[0]),                                 \
      "The projection module and the generated plasticity matrices disagree about the "            \
      "nodal point set of the volume.");
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CHECK_NODAL_SET)
#undef SEISSOL_CHECK_NODAL_SET

void setupCheckpointing(seissol::SeisSol& seissolInstance) {
  auto& checkpoint = seissolInstance.outputManager().getCheckpointManager();
  // the data of the run is checkpointed in its configuration; a checkpoint holds the cells of one
  // configuration only
  const auto& model = seissolInstance.parameters().model;
  const auto config = model.config;
  if (model.configs().size() > 1 &&
      (seissolInstance.parameters().output.checkpointParameters.enabled ||
       seissolInstance.checkpointLoadFile().has_value())) {
    logError() << "Checkpoints are not supported yet for runs whose cells compute in several "
                  "configurations.";
  }
  checkpoint.setConfig(config);

  {
    auto& storage = seissolInstance.memoryManager().ltsStorage();
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
    LTS::registerCheckpointVariables(checkpoint, storage, config);
  }

  {
    auto& storage = seissolInstance.memoryManager().drStorage();
    auto& dynrup = seissolInstance.memoryManager().drDescriptor();
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
    auto& storage = seissolInstance.memoryManager().surfaceStorage();
    std::vector<std::size_t> faceIdentifiers(storage.size(seissol::initializer::LayerMask(Ghost)));
    const auto* meshIds = storage.var<SurfaceLTS::MeshId>();
    const auto* sides = storage.var<SurfaceLTS::Side>();

    const auto& elements = seissolInstance.meshReader().getElements();
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
      // As for DR, a face is identified by its cell and its side -- by the global id of the cell,
      // since the mesh id is local to a rank: the same one names different cells on different
      // ranks, and a restart may distribute the cells differently.
      faceIdentifiers[i] =
          elements[meshIds[i]].globalId * Cell::NumFaces + static_cast<std::size_t>(sides[i]);
    }
    checkpoint.registerTree("surface", storage, faceIdentifiers);
    SurfaceLTS::registerCheckpointVariables(checkpoint, storage);
  }

  // the state of the derived outputs starts from its initial values, unless the checkpoint has it
  initializer::initializeDerivedState(seissolInstance);
  initializer::registerDerivedStateCheckpoints(checkpoint, seissolInstance);

  const auto& checkpointFile = seissolInstance.checkpointLoadFile();
  if (checkpointFile.has_value()) {
    const double time = seissolInstance.outputManager().loadCheckpoint(checkpointFile.value());
    seissolInstance.simulator().setCurrentTime(time);
  }

  if (seissolInstance.parameters().output.checkpointParameters.enabled) {
    // FIXME: for now, we allow only _one_ checkpoint interval which checkpoints everything
    // existent
    seissolInstance.outputManager().setupCheckpoint(
        seissolInstance.parameters().output.checkpointParameters.interval);
  }
}

/// The output of a quantity at the points of a subcell: called with the values to fill, the index
/// of the cell (or face) among the ones written, and the subcell.
using PointOutput = std::function<void(double*, std::size_t, std::size_t)>;

/// The outputs of the cells of one configuration, with their names, in the order of writing.
using NamedOutputs = std::vector<std::pair<std::string, PointOutput>>;

/**
 * Adds the outputs of the configurations `configs` of the run to `writer`, one per name. The names
 * follow the order of `configs`, and the order of the outputs of each. A cell (or face) writes an
 * output in its configuration, given by `configOf` for its index among the written ones; if its
 * configuration has no output of that name, it writes NaN at its `pointsPerSubcell` points.
 */
void addOutputsByName(io::instance::geometry::GeometryWriter& writer,
                      const std::vector<NamedOutputs>& outputsOfConfigs,
                      const std::vector<ConfigId>& configs,
                      const std::shared_ptr<const std::vector<ConfigId>>& configOf,
                      std::size_t pointsPerSubcell) {
  if (configs.size() == 1) {
    for (const auto& [name, output] : outputsOfConfigs.at(configs[0])) {
      writer.addGeometryOutput<double>(name, {}, false, output);
    }
    return;
  }

  std::vector<std::string> names;
  // by name, then by configuration
  std::vector<std::vector<PointOutput>> outputs;
  for (const auto config : configs) {
    for (const auto& [name, output] : outputsOfConfigs.at(config)) {
      auto found = std::find(names.begin(), names.end(), name);
      if (found == names.end()) {
        names.push_back(name);
        outputs.emplace_back(builtConfigCount());
        found = std::prev(names.end());
      }
      outputs[std::distance(names.begin(), found)][config] = output;
    }
  }

  for (std::size_t i = 0; i < names.size(); ++i) {
    writer.addGeometryOutput<double>(
        names[i],
        {},
        false,
        [configOf, outputsOfName = std::move(outputs[i]), pointsPerSubcell](
            double* target, std::size_t index, std::size_t subcell) {
          const auto& output = outputsOfName[configOf->at(index)];
          if (output) {
            output(target, index, subcell);
          } else {
            std::fill_n(target, pointsPerSubcell, std::numeric_limits<double>::quiet_NaN());
          }
        });
  }
}

/// Adds the position of the configuration of every cell (or face) in the list of configurations
/// of the run, `configs`, given the configuration of each by `configOf`, as the cell data
/// `config`; only if the run has several configurations.
void addConfigData(io::instance::geometry::GeometryWriter& writer,
                   const std::vector<ConfigId>& configs,
                   const std::shared_ptr<const std::vector<ConfigId>>& configOf) {
  if (configs.size() > 1) {
    std::vector<int> positions(builtConfigCount());
    for (std::size_t i = 0; i < configs.size(); ++i) {
      positions[configs[i]] = static_cast<int>(i);
    }
    writer.addCellData<int>(
        "config",
        {},
        true,
        [configOf, positions](int* target, std::size_t index, std::size_t /*subcell*/) {
          target[0] = positions[configOf->at(index)];
        });
  }
}

/// What the wave field output of all configurations shares.
struct WaveFieldOutputSetup {
  // the cells written, by their mesh ids
  std::shared_ptr<const std::vector<std::size_t>> cells;
  std::vector<io::instance::geometry::Subcell<Cell::Dim>> subcells;
  // the output points of a subcell, in its reference coordinates
  std::vector<std::array<double, Cell::Dim>> dataBase;
  std::uint32_t order{};
  std::uint32_t dataOrder{};
  projection::Target projectionTarget{};
  // the configured derived-output program, if any
  std::shared_ptr<const expr::Program> script;
};

/// Adds the outputs of `output`, copied from its buffer when written, to `outputs`, and `output`
/// itself to `derived`.
void collectDerived(NamedOutputs& outputs,
                    std::vector<std::shared_ptr<initializer::DerivedOutput>>& derived,
                    std::shared_ptr<initializer::DerivedOutput> output) {
  if (output == nullptr) {
    return;
  }
  for (std::size_t i = 0; i < output->names().size(); ++i) {
    outputs.emplace_back(output->names()[i],
                         [output, i](double* target, std::size_t index, std::size_t subcell) {
                           output->copy(i, target, index, subcell);
                         });
  }
  derived.push_back(std::move(output));
}

/// All outputs of a program are computed at once, before `writer` pulls them name by name; a
/// program with state follows every time step of the cells instead.
void evaluateDerived(seissol::SeisSol& seissolInstance,
                     io::instance::geometry::GeometryWriter& writer,
                     const std::vector<std::shared_ptr<initializer::DerivedOutput>>& derived) {
  writer.addHook([derived](std::size_t /*counter*/, double time) {
    for (const auto& output : derived) {
      output->write(time);
    }
  });
  for (const auto& output : derived) {
    if (output->accumulates()) {
      seissolInstance.timeManager().addCorrectionHook(
          [output](std::size_t layer, double time) { output->step(layer, time); });
    }
  }
}

/// The outputs of the wave field of the cells of the configuration `Cfg`: the built-in ones and
/// those of the configured program, each set computed by one derived program into a buffer that
/// `derived` collects, and copied from there when written.
template <typename Cfg>
NamedOutputs waveFieldOutputsOf(seissol::SeisSol& seissolInstance,
                                const WaveFieldOutputSetup& setup,
                                std::vector<std::shared_ptr<initializer::DerivedOutput>>& derived) {
  using MaterialT = model::MaterialOf<Cfg>;
  const auto& seissolParams = seissolInstance.parameters();
  const auto& parameters = seissolParams.output.waveFieldParameters;

  initializer::DerivedGeometry geometry;
  geometry.subcells = setup.subcells;
  geometry.dataBase = setup.dataBase;
  geometry.dataOrder = setup.dataOrder;
  geometry.order = Cfg::ConvergenceOrder;
  geometry.target = setup.projectionTarget;
  geometry.nodalSet = PlasticityNodalSet<Cfg>;

  initializer::WaveFieldSelection selection;
  selection.quantities.assign(MaterialT::Quantities.begin(), MaterialT::Quantities.end());
  selection.velocityOffset = MaterialT::VelocityOffset;
  selection.outputMask = parameters.outputMask;
  selection.integrationMask = parameters.integrationMask;
  selection.strain = parameters.computeStrain;
  selection.rotation = parameters.computeRotation;
  if (seissolParams.model.plasticity) {
    selection.plasticQuantities.assign(seissol::model::PlasticityData<Cfg>::Quantities.begin(),
                                       seissol::model::PlasticityData<Cfg>::Quantities.end());
    selection.plasticityMask.assign(parameters.plasticityMask.begin(),
                                    parameters.plasticityMask.end());
  }

  std::vector<const expr::Program*> programs;
  const auto builtIn = initializer::waveFieldProgram(selection);
  programs.push_back(&builtIn);
  if (setup.script != nullptr) {
    programs.push_back(setup.script.get());
  }

  NamedOutputs outputs;
  for (const auto* program : programs) {
    if (!program->outputs().empty()) {
      collectDerived(outputs,
                     derived,
                     initializer::makeDerivedVolumeOutput<Cfg>(
                         seissolInstance, setup.cells, geometry, *program));
    }
  }
  return outputs;
}

/// Sets up the output of the wave field. Every cell is written in its configuration; with several
/// configurations in the run, a quantity is written once for the cells of all of them.
void setupWaveFieldOutput(seissol::SeisSol& seissolInstance,
                          const initializer::OutputRegions& regions) {
  const auto& seissolParams = seissolInstance.parameters();
  const auto& parameters = seissolParams.output.waveFieldParameters;
  auto& memoryManager = seissolInstance.memoryManager();
  const auto* meshReader = &seissolInstance.meshReader();

  const auto orderIO = parameters.vtkorder;
  const auto order = static_cast<uint32_t>(std::max(0, orderIO));

  std::vector<std::size_t> celllist;
  celllist.reserve(meshReader->getElements().size());
  if (parameters.bounds.enabled || !parameters.groups.empty() ||
      regions.restricts(initializer::OutputRegions::WaveField)) {
    const auto& elements = meshReader->getElements();
    const auto& vertexArray = meshReader->getVertices();
    const auto inScriptedRegion = regions.select(
        initializer::OutputRegions::WaveField,
        elements.size(),
        Cell::NumVertices,
        [&](std::size_t cell, std::size_t corner) {
          return vertexArray[elements[cell].vertices[corner]].coords;
        },
        [&](std::size_t cell) { return elements[cell].group; });
    for (std::size_t i = 0; i < elements.size(); ++i) {
      const auto& element = elements[i];
      const auto& vertex0 = vertexArray[element.vertices[0]].coords;
      const auto& vertex1 = vertexArray[element.vertices[1]].coords;
      const auto& vertex2 = vertexArray[element.vertices[2]].coords;
      const auto& vertex3 = vertexArray[element.vertices[3]].coords;
      const bool inGroup = parameters.groups.empty() ||
                           parameters.groups.find(element.group) != parameters.groups.end();
      const bool inRegion = !parameters.bounds.enabled ||
                            (parameters.bounds.contains(vertex0[0], vertex0[1], vertex0[2]) ||
                             parameters.bounds.contains(vertex1[0], vertex1[1], vertex1[2]) ||
                             parameters.bounds.contains(vertex2[0], vertex2[1], vertex2[2]) ||
                             parameters.bounds.contains(vertex3[0], vertex3[1], vertex3[2]));
      if (inGroup && inRegion && inScriptedRegion[i]) {
        celllist.push_back(i);
      }
    }
  } else {
    for (std::size_t i = 0; i < meshReader->getElements().size(); ++i) {
      celllist.push_back(i);
    }
  }

  WaveFieldOutputSetup setup;
  // the output points, as the state of the derived outputs is laid out for them
  const auto points = initializer::waveFieldGeometry(parameters);
  setup.order = order;
  setup.dataOrder = static_cast<std::uint32_t>(points.dataOrder);
  setup.dataBase = points.dataBase;
  setup.subcells = points.subcells;
  setup.projectionTarget = points.target;
  const auto trueOrder = order > 0 ? order : 1;
  const auto trueBase = io::instance::geometry::pointsTetrahedron(trueOrder);
  const auto truePoints = io::instance::geometry::applyMaps(setup.subcells, trueBase);

  const auto format = orderIO < 0 ? io::instance::geometry::WriterFormat::Xdmf
                                  : io::instance::geometry::WriterFormat::Vtk;

  const auto config = io::instance::geometry::WriterConfig{
      order,
      format,
      seissolParams.output.xdmfWriterBackend == seissol::initializer::parameters::XdmfBackend::Posix
          ? io::instance::geometry::WriterBackend::Binary
          : io::instance::geometry::WriterBackend::Hdf5,
      io::instance::geometry::supportedWriterGroup(
          writerGroupOf(parameters.timeSeries), format, "wavefield"),
      seissolParams.output.hdfcompress};

  const auto cells = std::make_shared<const std::vector<std::size_t>>(std::move(celllist));
  setup.cells = cells;

  io::writer::ScheduledWriter schedWriter;
  schedWriter.name = "wavefield";
  schedWriter.interval = parameters.interval;
  auto writer = io::instance::geometry::GeometryWriter(
      "wavefield",
      cells->size(),
      io::instance::geometry::Shape::Tetrahedron,
      config,
      setup.subcells.size(),

      [=](double* target, std::size_t index, std::size_t subcell) {
        const auto transform =
            seissol::geometry::AffineTransform::fromMeshCell(cells->at(index), *meshReader);

        for (std::size_t i = 0; i < truePoints[subcell].size(); ++i) {
          const auto xyz = transform.refToSpace(truePoints[subcell][i]);
          std::copy_n(xyz.begin(), Cell::Dim, &target[i * 3]);
        }
      });

  const auto rank = seissol::Mpi::mpi.rank();
  writer.addCellData<int>(
      "partition", {}, true, [=](int* target, std::size_t /*index*/, std::size_t /*subcell*/) {
        target[0] = rank;
      });

  writer.addCellData<uint64_t>(
      "clustering", {}, true, [=](uint64_t* target, std::size_t index, std::size_t /*subcell*/) {
        target[0] = meshReader->getElements()[cells->at(index)].clusterId;
      });

  writer.addCellData<std::size_t>(
      "global-id", {}, true, [=](std::size_t* target, std::size_t index, std::size_t /*subcell*/) {
        target[0] = meshReader->getElements()[cells->at(index)].globalId;
      });

  // the cells of every configuration of the run
  const auto configs = seissolParams.model.configs();
  auto& ltsStorage = memoryManager.ltsStorage();
  auto& backmap = memoryManager.backmap();
  std::vector<ConfigId> configOfCell(cells->size());
  for (std::size_t index = 0; index < cells->size(); ++index) {
    configOfCell[index] =
        ltsStorage.lookup<LTS::SecondaryInformation>(backmap.get(cells->at(index))).configId;
  }
  const auto configOf = std::make_shared<const std::vector<ConfigId>>(std::move(configOfCell));
  addConfigData(writer, configs, configOf);

  if (!parameters.script.empty()) {
    setup.script =
        std::make_shared<const expr::Program>(initializer::loadDerivedProgram(parameters.script));
  }

  std::vector<std::shared_ptr<initializer::DerivedOutput>> derived;
  std::vector<NamedOutputs> outputs(builtConfigCount());
  for (const auto runConfig : configs) {
    dispatchConfig(runConfig, [&](auto cfg) {
      outputs[runConfig] = waveFieldOutputsOf<decltype(cfg)>(seissolInstance, setup, derived);
    });
  }
  addOutputsByName(writer, outputs, configs, configOf, setup.dataBase.size());
  evaluateDerived(seissolInstance, writer, derived);

  schedWriter.planWrite = writer.makeWriter();
  seissolInstance.outputManager().addOutput(schedWriter);
}

/// What the free surface output of all configurations shares.
struct SurfaceOutputSetup {
  // the faces written, by their index in the free surface integrator
  std::shared_ptr<const std::vector<std::size_t>> faces;
  std::vector<io::instance::geometry::Subcell<2>> subcells;
  // the output points of a subcell, in its reference coordinates
  std::vector<std::array<double, 2>> dataBase;
  std::uint32_t order{};
  std::uint32_t dataOrder{};
  projection::Target projectionTarget{};
  // the configured derived-output program, if any
  std::shared_ptr<const expr::Program> script;
};

/// The outputs of the free surface of the faces of cells of the configuration `Cfg`: the
/// built-in ones and those of the configured program, as for the wave field.
template <typename Cfg>
NamedOutputs surfaceOutputsOf(seissol::SeisSol& seissolInstance,
                              const SurfaceOutputSetup& setup,
                              std::vector<std::shared_ptr<initializer::DerivedOutput>>& derived) {
  using MaterialT = model::MaterialOf<Cfg>;
  const auto& parameters = seissolInstance.parameters().output.freeSurfaceParameters;

  initializer::DerivedSurfaceGeometry geometry;
  geometry.subcells = setup.subcells;
  geometry.dataBase = setup.dataBase;
  geometry.dataOrder = setup.dataOrder;
  geometry.order = Cfg::ConvergenceOrder;
  geometry.target = setup.projectionTarget;
  geometry.nodalSet = PlasticityNodalSet<Cfg>;

  std::vector<const expr::Program*> programs;
  const auto builtIn = initializer::surfaceProgram(
      std::vector<std::string>(MaterialT::Quantities.begin(), MaterialT::Quantities.end()),
      parameters.outputMask);
  programs.push_back(&builtIn);
  if (setup.script != nullptr) {
    programs.push_back(setup.script.get());
  }

  NamedOutputs outputs;
  for (const auto* program : programs) {
    if (!program->outputs().empty()) {
      collectDerived(outputs,
                     derived,
                     initializer::makeDerivedSurfaceOutput<Cfg>(
                         seissolInstance, setup.faces, geometry, *program));
    }
  }
  return outputs;
}

/// Sets up the output of the free surface. Every face is written in the configuration of its
/// cell; with several configurations in the run, a quantity is written once for the faces of all
/// of them.
void setupSurfaceOutput(seissol::SeisSol& seissolInstance,
                        const initializer::OutputRegions& regions) {
  const auto& seissolParams = seissolInstance.parameters();
  const auto& parameters = seissolParams.output.freeSurfaceParameters;
  auto& memoryManager = seissolInstance.memoryManager();

  const auto orderIO = parameters.vtkorder;
  const auto order = static_cast<std::uint32_t>(std::max(0, orderIO));

  auto* freeSurfaceIntegrator = &seissolInstance.freeSurfaceIntegrator();
  const auto* meshReader = &seissolInstance.meshReader();
  io::writer::ScheduledWriter schedWriter;
  schedWriter.name = "surface";
  schedWriter.interval = parameters.interval;
  auto* surfaceMeshIds = freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::MeshId>();
  auto* surfaceMeshSides = freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::Side>();
  auto* surfaceLocationFlag =
      freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::LocationFlag>();

  SurfaceOutputSetup setup;
  // the faces in the region of the output, if a model restricts it
  const auto faceCount = freeSurfaceIntegrator->backmap.size();
  const auto inScriptedRegion = regions.select(
      initializer::OutputRegions::Surface,
      faceCount,
      Face::NumVertices,
      [&](std::size_t index, std::size_t corner) {
        const auto face = freeSurfaceIntegrator->backmap[index];
        const auto transform = seissol::geometry::AffineFaceTransform::fromMeshCell(
            surfaceMeshIds[face], surfaceMeshSides[face], *meshReader);
        // the corners of the reference triangle
        const auto xyz = transform.refToSpace(seissol::geometry::FaceTransform::FaceVectorT(
            corner == 1 ? 1.0 : 0.0, corner == 2 ? 1.0 : 0.0));
        return std::array<double, 3>{xyz(0), xyz(1), xyz(2)};
      },
      [&](std::size_t index) {
        return meshReader->getElements()[surfaceMeshIds[freeSurfaceIntegrator->backmap[index]]]
            .group;
      });
  std::vector<std::size_t> faceList;
  faceList.reserve(faceCount);
  for (std::size_t index = 0; index < faceCount; ++index) {
    if (inScriptedRegion[index]) {
      faceList.push_back(index);
    }
  }
  const auto faces = std::make_shared<const std::vector<std::size_t>>(std::move(faceList));
  setup.faces = faces;

  // the output points, as the state of the derived outputs is laid out for them
  const auto points = initializer::surfaceGeometry(parameters);
  setup.order = order;
  setup.dataOrder = static_cast<std::uint32_t>(points.dataOrder);
  setup.dataBase = points.dataBase;
  setup.subcells = points.subcells;
  setup.projectionTarget = points.target;
  const auto trueOrder = order > 0 ? order : 1;
  const auto trueBase = io::instance::geometry::pointsTriangle(trueOrder);
  const auto truePoints = io::instance::geometry::applyMaps(setup.subcells, trueBase);

  const auto format = orderIO < 0 ? io::instance::geometry::WriterFormat::Xdmf
                                  : io::instance::geometry::WriterFormat::Vtk;

  const auto config = io::instance::geometry::WriterConfig{
      order,
      format,
      seissolParams.output.xdmfWriterBackend == seissol::initializer::parameters::XdmfBackend::Posix
          ? io::instance::geometry::WriterBackend::Binary
          : io::instance::geometry::WriterBackend::Hdf5,
      io::instance::geometry::supportedWriterGroup(
          writerGroupOf(parameters.timeSeries), format, "free surface"),
      seissolParams.output.hdfcompress};

  auto writer = io::instance::geometry::GeometryWriter(
      "surface",
      faces->size(),
      io::instance::geometry::Shape::Triangle,
      config,
      setup.subcells.size(),

      [=](double* target, std::size_t index, std::size_t subcell) {
        auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[faces->at(index)]];
        auto side = surfaceMeshSides[freeSurfaceIntegrator->backmap[faces->at(index)]];
        const auto face =
            seissol::geometry::AffineFaceTransform::fromMeshCell(meshId, side, *meshReader);

        for (std::size_t i = 0; i < truePoints[subcell].size(); ++i) {
          const auto xyz = face.refToSpace(
              seissol::geometry::FaceTransform::FaceVectorT(truePoints[subcell][i].data()));
          for (std::size_t d = 0; d < Cell::Dim; ++d) {
            target[i * 3 + d] = xyz(d);
          }
        }
      });

  const auto rank = seissol::Mpi::mpi.rank();
  writer.addCellData<int>(
      "partition", {}, true, [=](int* target, std::size_t /*index*/, std::size_t /*subcell*/) {
        target[0] = rank;
      });

  // four bytes wide, as it has always been written: readers such as seissolxdmf know no narrower
  // integers
  writer.addCellData<std::uint32_t>(
      "locationFlag",
      {},
      true,
      [=](std::uint32_t* target, std::size_t index, std::size_t /*subcell*/) {
        target[0] = surfaceLocationFlag[freeSurfaceIntegrator->backmap[faces->at(index)]];
      });

  writer.addCellData<std::size_t>(
      "global-id", {}, true, [=](std::size_t* target, std::size_t index, std::size_t /*subcell*/) {
        const auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[faces->at(index)]];
        const auto side = surfaceMeshSides[freeSurfaceIntegrator->backmap[faces->at(index)]];
        target[0] = meshReader->getElements()[meshId].globalId * 4 + side;
      });

  // the faces of every configuration of the run
  const auto configs = seissolParams.model.configs();
  auto& ltsStorage = memoryManager.ltsStorage();
  auto& backmap = memoryManager.backmap();
  std::vector<ConfigId> configOfFace(faces->size());
  for (std::size_t index = 0; index < configOfFace.size(); ++index) {
    const auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[faces->at(index)]];
    configOfFace[index] =
        ltsStorage.lookup<LTS::SecondaryInformation>(backmap.get(meshId)).configId;
  }
  const auto configOf = std::make_shared<const std::vector<ConfigId>>(std::move(configOfFace));
  addConfigData(writer, configs, configOf);

  if (!parameters.script.empty()) {
    setup.script =
        std::make_shared<const expr::Program>(initializer::loadDerivedProgram(parameters.script));
  }

  std::vector<std::shared_ptr<initializer::DerivedOutput>> derived;
  std::vector<NamedOutputs> outputs(builtConfigCount());
  for (const auto runConfig : configs) {
    dispatchConfig(runConfig, [&](auto cfg) {
      outputs[runConfig] = surfaceOutputsOf<decltype(cfg)>(seissolInstance, setup, derived);
    });
  }
  addOutputsByName(writer, outputs, configs, configOf, setup.dataBase.size());
  evaluateDerived(seissolInstance, writer, derived);

  schedWriter.planWrite = writer.makeWriter();
  seissolInstance.outputManager().addOutput(schedWriter);
}

void setupOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.parameters();
  auto& memoryManager = seissolInstance.memoryManager();

  // the regions the mesh outputs are restricted to by a model, if any
  const initializer::OutputRegions regions(seissolParams.output.regionFileName);

  if (seissolParams.output.waveFieldParameters.enabled) {
    setupWaveFieldOutput(seissolInstance, regions);
  }

  if (seissolParams.output.freeSurfaceParameters.enabled) {
    setupSurfaceOutput(seissolInstance, regions);
  }

  if (seissolParams.output.receiverParameters.enabled) {
    auto& receiverWriter = seissolInstance.receiverWriter();
    // Initialize receiver output
    receiverWriter.init(seissolParams.output.prefix,
                        seissolParams.timeStepping.endTime,
                        seissolParams.output.receiverParameters);
    receiverWriter.addPoints(seissolInstance.meshReader(), memoryManager.backmap());
    seissolInstance.timeManager().setReceiverClusters(receiverWriter);
  }

  if (seissolParams.output.energyParameters.enabled) {
    auto& energyOutput = seissolInstance.energyOutput();

    energyOutput.init(memoryManager.drStorage(),
                      seissolInstance.meshReader(),
                      memoryManager.ltsStorage(),
                      seissolParams.model.plasticity,
                      seissolParams.output.prefix,
                      seissolParams.output.energyParameters);
  }

  seissolInstance.flopCounter().init(seissolParams.output.prefix);

  seissolInstance.analysisWriter().init(&seissolInstance.meshReader(), seissolParams.output.prefix);
}

void initFaultOutputManager(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.parameters();

  const auto& backupTimeStamp = seissolInstance.backupTimeStamp();

  auto* faultOutputManager = seissolInstance.memoryManager().faultOutputManager();

  if (seissolParams.drParameters.isDynamicRuptureEnabled) {
    auto& ltsStorage = seissolInstance.memoryManager().ltsStorage();
    auto& backmap = seissolInstance.memoryManager().backmap();
    auto& drStorage = seissolInstance.memoryManager().drStorage();

    faultOutputManager->setInputParam(seissolInstance.meshReader());
    faultOutputManager->setLtsData(ltsStorage, backmap, drStorage);
    faultOutputManager->setBackupTimeStamp(backupTimeStamp);
    faultOutputManager->init();
  }

  seissolInstance.timeManager().setFaultOutputManager(faultOutputManager);
}

void enableFreeSurfaceOutput(seissol::SeisSol& seissolInstance) {}

} // namespace

void seissol::initializer::initprocedure::initIO(seissol::SeisSol& seissolInstance) {
  const auto rank = Mpi::mpi.rank();
  logInfo() << "Begin init output.";

  const auto& seissolParams = seissolInstance.parameters();
  const filesystem::path outputPath(seissolParams.output.prefix);
  const auto outputDir = filesystem::directory_entry(outputPath.parent_path());
  if (!filesystem::exists(outputDir)) {
    logWarning() << "Output directory does not exist yet. We therefore create it now.";
    if (rank == 0) {
      filesystem::create_directory(outputDir);
    }
  }
  seissol::Mpi::barrier(Mpi::mpi.comm());

  // recorded while reading the mesh and laying out the clusters, before there was a directory to
  // write them to
  seissolInstance.miniSeisSolWriter().write(seissolParams.output.prefix);
  seissolInstance.timeManager().writeClustering();

  enableFreeSurfaceOutput(seissolInstance);
  initFaultOutputManager(seissolInstance);
  setupCheckpointing(seissolInstance);
  setupOutput(seissolInstance);
  logInfo() << "End init output.";
}
