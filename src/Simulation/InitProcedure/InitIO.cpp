// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitIO.h"

#include "Alignment.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Constants.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Common/Real.h"
#include "Core/Geometry/CellTransform.h"
#include "Core/Geometry/FaceTransform.h"
#include "Core/Geometry/MeshDefinition.h"
#include "Core/IO/Instance/Geometry/Geometry.h"
#include "Core/IO/Instance/Geometry/Points.h"
#include "Core/IO/Instance/Geometry/Refinement.h"
#include "Core/IO/Instance/Geometry/Typedefs.h"
#include "Core/IO/Writer/Writer.h"
#include "Core/IO/WriterParameters.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Descriptor/Surface.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Monitoring/MiniSeisSolWriter.h"
#include "Core/Numerical/Projection.h"
#include "Core/Parallel/MPI.h"
#include "Core/Parallel/ParallelFS.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "Features/Plasticity/Solver/Storage.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Simulation/Parameters/OutputParameters.h"
#include "Simulation/SeisSol.h"
#include "WavePropagation/Output/FreeSurfaceIntegrator.h"
#include "WavePropagation/Output/Receiver.h"
#include "WavePropagation/Output/ReceiverParameters.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <utils/logger.h>
#include <variant>
#include <vector>

namespace {
using namespace seissol;

//! @brief The file grouping a time series mode asks the writer for; none if there is no mode.
std::optional<io::instance::geometry::WriterGroup>
    writerGroupOf(std::optional<seissol::initializer::parameters::TimeSeriesMode> mode) {
  if (!mode.has_value()) {
    return std::nullopt;
  }
  switch (mode.value()) {
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

// The projection matrices are generated for every convergence order up to the one of the
// configuration, so that a per-cell order (cf. #1421) only requires selecting a different entry at
// run time.
constexpr std::size_t MinProjectionOrder = 1;

// Which nodal points the plastic strain lives on is a build option (PLASTICITY_METHOD): "nb"
// uses a unisolvent warp&blend set, "ip" the conical-product quadrature points, "sq" the points of
// a fully symmetric quadrature rule. Read that back off the generated matrices instead of
// duplicating the CMake variable: it is the set whose points are the generated vNodes. (Their
// number alone does not tell, e.g. "nb" and "sq" both have four points for order 2.)
template <typename Cfg>
projection::NodalSet plasticityNodalSet() {
  const auto nodes = init::vNodes<Cfg>::view::create(init::vNodes<Cfg>::Values);
  std::vector<std::array<double, Cell::Dim>> points(
      static_cast<std::size_t>(tensor::vNodes<Cfg>::Shape[0]));
  for (std::size_t p = 0; p < points.size(); ++p) {
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      points[p][d] = nodes.isInRange(p, d) ? nodes(p, d) : 0.0;
    }
  }
  // the generated values are stored in the precision of the configuration
  constexpr double Tolerance = std::max(1e-10, 10.0 * std::numeric_limits<Real<Cfg>>::epsilon());
  const auto set = projection::identifyNodalSet3D(Cfg::ConvergenceOrder, points, Tolerance);
  if (!set.has_value()) {
    logError() << "The generated plasticity points (vNodes) are none of the nodal point sets of "
                  "the projection module.";
  }
  return set.value_or(projection::NodalSet::WarpBlend);
}

/**
 * The padded leading dimension of a generated projection tensor, i.e. the stride between two
 * consecutive basis functions. Yateto stores these matrices as [point][basisFunction] with an
 * aligned stride on the point dimension; we read the padding back off the generated metadata
 * instead of re-deriving it from the alignment.
 */
template <typename Cfg, typename TensorT>
std::size_t projectionStride(std::size_t degree) {
  const auto index = TensorT::index(Cfg::ConvergenceOrder, degree);
  return TensorT::Size[index] / TensorT::Shape[index][1];
}

//! The affine embedding of the reference triangle into the given side of the reference tetrahedron.
seissol::numerical::AffineMap<2, 3> faceEmbedding(std::size_t side) {
  const std::array<std::array<double, 2>, 3> corners = {
      std::array<double, 2>{0, 0}, std::array<double, 2>{1, 0}, std::array<double, 2>{0, 1}};
  const auto faceMap = seissol::geometry::ReferenceFaceMap(side);
  std::vector<std::array<double, 3>> vertices;
  vertices.reserve(corners.size());
  for (const auto& chiTau : corners) {
    const auto xez =
        faceMap.faceToCell(seissol::geometry::ReferenceFaceMap::FaceVectorT(chiTau.data()));
    vertices.push_back({xez(0), xez(1), xez(2)});
  }
  return seissol::numerical::AffineMap<2, 3>::fromVertices(vertices);
}

void setupCheckpointing(seissol::SeisSol& seissolInstance) {
  auto& checkpoint = seissolInstance.outputManager().getCheckpointManager();
  // the data of the run is checkpointed in its configuration; a checkpoint holds the cells of one
  // configuration only
  const auto& runConfigs = seissolInstance.parameters().runConfigs;
  const auto config = runConfigs.config;
  // (writing them is rejected with the parameters already)
  if (runConfigs.configs().size() > 1 && seissolInstance.checkpointLoadFile().has_value()) {
    logError() << "Checkpoints are not supported yet for runs whose cells compute in several "
                  "configurations.";
  }
  checkpoint.setConfig(config);

  // every copy of a cell or face is a data point of its own
  const auto& simulationSettings = seissolInstance.memoryManager().simulationSettings();
  checkpoint.setSimulations(simulationSettings.simulationCount, simulationSettings.fusedCount);
  const auto copyCount = simulationSettings.copyCount();

  {
    auto& storage = seissolInstance.memoryManager().ltsStorage();
    std::vector<std::size_t> globalIds(storage.size(seissol::initializer::LayerMask(Ghost)));
    std::size_t offset = 0;
    for (const auto& layer : storage.leaves(Ghost)) {

#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < layer.size(); ++i) {
        const auto& secondaryInformation = layer.var<LTS::SecondaryInformation>()[i];
        const auto globalId =
            seissolInstance.meshReader().getElements()[secondaryInformation.meshId].globalId;
        globalIds[offset + i] = globalId * copyCount + secondaryInformation.copyId;
      }
      offset += layer.size();
    }
    checkpoint.registerTree("lts", storage, globalIds);
    LTS::registerCheckpointVariables(checkpoint, storage, config);
    LTSPlasticity::registerCheckpointVariables(checkpoint, storage);
  }

  {
    auto& storage = seissolInstance.memoryManager().drStorage();
    const auto& dynrup = seissolInstance.memoryManager().drDescriptor();
    std::vector<std::size_t> faceIdentifiers(storage.size(seissol::initializer::LayerMask(Ghost)));
    const auto* drFaceInformation = storage.var<DynamicRupture::FaceInformation>();

#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
      const auto faultFace = drFaceInformation[i].meshFace;
      const auto& fault = seissolInstance.meshReader().getFault()[faultFace];
      // take the positive cell and side as fault face identifier
      // (should result in roughly twice as large numbers as when indexing all faces; cf.
      // handshake theorem)
      faceIdentifiers[i] =
          (fault.globalId * 4 + fault.side) * copyCount + drFaceInformation[i].copy;
    }
    checkpoint.registerTree("dynrup", storage, faceIdentifiers);
    dynrup.registerCheckpointVariables(checkpoint, storage);
  }

  {
    auto& storage = seissolInstance.memoryManager().surfaceStorage();
    std::vector<std::size_t> faceIdentifiers(storage.size(seissol::initializer::LayerMask(Ghost)));
    const auto* meshIds = storage.var<SurfaceLTS::MeshId>();
    const auto* copyIds = storage.var<SurfaceLTS::CopyId>();
    const auto* sides = storage.var<SurfaceLTS::Side>();

    const auto& elements = seissolInstance.meshReader().getElements();
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
      // As for DR, a face is identified by its cell and its side -- by the global id of the cell,
      // since the mesh id is local to a rank: the same one names different cells on different
      // ranks, and a restart may distribute the cells differently.
      faceIdentifiers[i] =
          (elements[meshIds[i]].globalId * Cell::NumFaces + static_cast<std::size_t>(sides[i])) *
              copyCount +
          copyIds[i];
    }
    checkpoint.registerTree("surface", storage, faceIdentifiers);
    SurfaceLTS::registerCheckpointVariables(checkpoint, storage);
  }

  const auto& checkpointFile = seissolInstance.checkpointLoadFile();
  if (checkpointFile.has_value()) {
    const double time = seissolInstance.outputManager().loadCheckpoint(checkpointFile.value());
    seissolInstance.simulator().setCurrentTime(time);
  }

  if (const auto& checkpoint = seissolInstance.parameters().output.checkpoint;
      checkpoint.has_value()) {
    // FIXME: for now, we allow only _one_ checkpoint interval which checkpoints everything
    // existent
    seissolInstance.outputManager().setupCheckpoint(checkpoint->interval);
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
      auto found = std::ranges::find(names, name);
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
};

/// The outputs of the wave field of the cells of the configuration `Cfg`.
template <typename Cfg>
NamedOutputs waveFieldOutputsOf(
    seissol::SeisSol& seissolInstance,
    const WaveFieldOutputSetup& setup,
    const seissol::initializer::parameters::WaveFieldOutputParameters& waveField) {
  using real = Real<Cfg>;
  using MaterialT = model::MaterialOf<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();
  // the projection matrices of the configuration
  constexpr std::size_t MaxProjectionOrder = Cfg::ConvergenceOrder;
  const auto& seissolParams = seissolInstance.parameters();
  auto& memoryManager = seissolInstance.memoryManager();
  const auto* ltsStorage = &memoryManager.ltsStorage();
  const auto* backmap = &memoryManager.backmap();
  const auto* meshReader = &seissolInstance.meshReader();

  // TODO(David): change Yateto/TensorForge interface to make padded sizes more accessible
  constexpr auto QDofSizePadded =
      tensor::Q<Cfg>::Size / tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  constexpr auto QDofPointsPadded = tensor::QStressNodal<Cfg>::Size /
                                    tensor::QStressNodal<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];

  const auto simulationCount = memoryManager.simulationSettings().simulationCount;
  const auto namewrap = [simulationCount](const std::string& name, std::size_t sim) {
    if (simulationCount > 1) {
      return name + "-" + std::to_string(sim + 1);
    } else {
      return name;
    }
  };

  const auto order = setup.order;
  const auto& cells = setup.cells;
  const auto& dataBase = setup.dataBase;
  const auto pointsPerSubcell = dataBase.size();

  const auto makeVolumeTable = [&](projection::Source source,
                                   std::optional<std::size_t> derivative) {
    projection::Spec spec;
    spec.source = source;
    spec.target = setup.projectionTarget;
    spec.nodalSet = plasticityNodalSet<Cfg>();
    spec.derivative = derivative;
    const auto stride = source == projection::Source::Nodal
                            ? projectionStride<Cfg, tensor::collnv<Cfg>>(order)
                            : projectionStride<Cfg, tensor::collvv<Cfg>>(order);
    return std::make_shared<projection::Table<3, 3, real>>(setup.subcells,
                                                           dataBase,
                                                           setup.dataOrder,
                                                           stride,
                                                           spec,
                                                           MinProjectionOrder,
                                                           MaxProjectionOrder);
  };

  const auto proj = makeVolumeTable(projection::Source::Modal, {});

  std::array<std::shared_ptr<projection::Table<3, 3, real>>, Cell::Dim> projD{};
  if (waveField.computeStrain || waveField.computeRotation) {
    for (std::size_t direction = 0; direction < Cell::Dim; ++direction) {
      projD[direction] = makeVolumeTable(projection::Source::Modal, direction);
    }
  }

  std::shared_ptr<projection::Table<3, 3, real>> projNodal;
  if (seissolParams.plasticity.enabled) {
    projNodal = makeVolumeTable(projection::Source::Nodal, {});
  }

  constexpr std::size_t MaxVtk3dPoints = tensor::vtk3d<Cfg>::Shape
      [(sizeof(tensor::vtk3d<Cfg>::Shape) / sizeof(tensor::vtk3d<Cfg>::Shape[0])) - 1][1];

  NamedOutputs outputs;

  for (std::size_t sim = 0; sim < simulationCount; ++sim) {
    // the copy of the cells which holds the simulation, and the fused simulation there
    const auto copy = sim / Cfg::NumSimulations;
    const auto fused = sim % Cfg::NumSimulations;

    const auto projectVolume =
        [=](double* target, const real* dofsSingleQuantity, const real* collvv) {
          runtime::kernel::projectBasisToVtkVolume vtkproj{};
          memory::AlignedArray<real, Cfg::NumSimulations> simselect{};
          alignas(Alignment) std::array<real, MaxVtk3dPoints> alignedTarget{};
          simselect[fused] = 1;
          vtkproj.simselect = runtime::init::simselect::view(Variant, simselect.data());
          vtkproj.qb = runtime::init::qb::view(Variant, dofsSingleQuantity);
          vtkproj.xv(order) = runtime::init::xv::view(Variant, order, alignedTarget.data());
          vtkproj.collvv(Cfg::ConvergenceOrder, order) =
              runtime::init::collvv::view(Variant, Cfg::ConvergenceOrder, order, collvv);
          vtkproj.execute(Variant, order);
          std::copy_n(alignedTarget.data(), pointsPerSubcell, target);
        };

    const auto projectVolumeDeriv = [=](double* target,
                                        const real* dofsSingleQuantity,
                                        std::size_t dir,
                                        std::size_t index,
                                        std::size_t subcell) {
      std::array<double, MaxVtk3dPoints> dataX{};
      std::array<double, MaxVtk3dPoints> dataY{};
      std::array<double, MaxVtk3dPoints> dataZ{};

      projectVolume(dataX.data(), dofsSingleQuantity, (*projD[0])(subcell, Cfg::ConvergenceOrder));
      projectVolume(dataY.data(), dofsSingleQuantity, (*projD[1])(subcell, Cfg::ConvergenceOrder));
      projectVolume(dataZ.data(), dofsSingleQuantity, (*projD[2])(subcell, Cfg::ConvergenceOrder));

      const auto transform =
          seissol::geometry::AffineTransform::fromMeshCell(cells->at(index), *meshReader);

      // IMPORTANT NOTE: we rely on the linearity of the cell transform in this place.
      // (the rows of the inverse Jacobian are grad xi, grad eta, grad zeta)
      const auto grad = transform.refToSpaceJacobianInverse(
          seissol::geometry::CellTransform::VectorEigenT(Cell::ReferenceBarycenter.data()));

      for (std::size_t i = 0; i < pointsPerSubcell; ++i) {
        target[i] = dataX[i] * grad(0, dir) + dataY[i] * grad(1, dir) + dataZ[i] * grad(2, dir);
      }
    };

    for (std::size_t quantity = 0; quantity < MaterialT::Quantities.size(); ++quantity) {

      if (waveField.outputMask[Variant][quantity]) {
        outputs.emplace_back(
            namewrap(MaterialT::Quantities[quantity], sim),
            [=](double* target, std::size_t index, std::size_t subcell) {
              const auto position = backmap->get(cells->at(index), copy);
              const auto* dofsAllQuantities = ltsStorage->lookup<LTS::Dofs>(Cfg(), position);
              const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
              projectVolume(target, dofsSingleQuantity, (*proj)(subcell, Cfg::ConvergenceOrder));
            });
      }

      if (waveField.integrationMask[Variant][quantity]) {
        outputs.emplace_back(
            namewrap("int-" + MaterialT::Quantities[quantity], sim),
            [=](double* target, std::size_t index, std::size_t subcell) {
              const auto position = backmap->get(cells->at(index), copy);
              const auto* dofsAllQuantities = ltsStorage->lookup<LTS::Integrals>(Cfg(), position);
              const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
              projectVolume(target, dofsSingleQuantity, (*proj)(subcell, Cfg::ConvergenceOrder));
            });
      }
    }

    using Idx = std::tuple<std::string, std::size_t, std::size_t>;

    if (waveField.computeStrain) {
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

        outputs.emplace_back(
            namewrap("eps" + name, sim),
            [=](double* target, std::size_t index, std::size_t subcell) {
              const auto position = backmap->get(cells->at(index), copy);
              const auto* dofsAllQuantities = ltsStorage->lookup<LTS::Integrals>(Cfg(), position);
              const auto* dofsSingleQuantity1 =
                  dofsAllQuantities + QDofSizePadded * (idx1 + MaterialT::VelocityOffset);
              projectVolumeDeriv(target, dofsSingleQuantity1, idx2, index, subcell);

              if (idx1 != idx2) {
                const auto* dofsSingleQuantity2 =
                    dofsAllQuantities + QDofSizePadded * (idx2 + MaterialT::VelocityOffset);
                std::array<double, MaxVtk3dPoints> itarget{};
                projectVolumeDeriv(itarget.data(), dofsSingleQuantity2, idx1, index, subcell);

                for (std::size_t i = 0; i < pointsPerSubcell; ++i) {
                  target[i] = (target[i] + itarget[i]) / 2;
                }
              }
            });
      }
    }

    if (waveField.computeRotation) {
      for (const auto& idxp : {Idx{"1", 2, 1}, Idx{"2", 0, 2}, Idx{"3", 1, 0}}) {
        const auto& name = std::get<0>(idxp);

        // need non-reference capture (due to lambda usage)
        const auto idx1 = std::get<1>(idxp);
        const auto idx2 = std::get<2>(idxp);

        // compute d_i2 v_i1 - d_i1 v_i2

        outputs.emplace_back(
            namewrap("rot" + name, sim),
            [=](double* target, std::size_t index, std::size_t subcell) {
              const auto position = backmap->get(cells->at(index), copy);
              const auto* dofsAllQuantities = ltsStorage->lookup<LTS::Dofs>(Cfg(), position);
              const auto* dofsSingleQuantity1 =
                  dofsAllQuantities + QDofSizePadded * (idx1 + MaterialT::VelocityOffset);
              projectVolumeDeriv(target, dofsSingleQuantity1, idx2, index, subcell);

              const auto* dofsSingleQuantity2 =
                  dofsAllQuantities + QDofSizePadded * (idx2 + MaterialT::VelocityOffset);
              std::array<double, MaxVtk3dPoints> itarget{};
              projectVolumeDeriv(itarget.data(), dofsSingleQuantity2, idx1, index, subcell);

              for (std::size_t i = 0; i < tensor::vtk3d<Cfg>::Shape[order][1]; ++i) {
                target[i] -= itarget[i];
              }
            });
      }
    }

    if (seissolParams.plasticity.enabled) {
      for (std::size_t quantity = 0;
           quantity < seissol::model::PlasticityData<Cfg>::Quantities.size();
           ++quantity) {
        if (waveField.plasticityMask[quantity]) {
          outputs.emplace_back(
              namewrap(seissol::model::PlasticityData<Cfg>::Quantities[quantity], sim),
              [=](double* target, std::size_t index, std::size_t subcell) {
                const auto position = backmap->get(cells->at(index), copy);
                const auto* dofsAllQuantities =
                    ltsStorage->lookup<LTSPlasticity::PStrain>(Cfg(), position);
                const auto* pointsSingleQuantity = dofsAllQuantities + QDofPointsPadded * quantity;
                runtime::kernel::projectNodalToVtkVolume vtkproj{};
                memory::AlignedArray<real, Cfg::NumSimulations> simselect{};
                alignas(Alignment) std::array<real, MaxVtk3dPoints> alignedTarget{};
                simselect[fused] = 1;
                vtkproj.simselect = runtime::init::simselect::view(Variant, simselect.data());
                vtkproj.qn = runtime::init::qn::view(Variant, pointsSingleQuantity);
                vtkproj.xv(order) = runtime::init::xv::view(Variant, order, alignedTarget.data());
                vtkproj.collnv(Cfg::ConvergenceOrder, order) =
                    runtime::init::collnv::view(Variant,
                                                Cfg::ConvergenceOrder,
                                                order,
                                                (*projNodal)(subcell, Cfg::ConvergenceOrder));
                vtkproj.execute(Variant, order);
                std::copy_n(alignedTarget.data(), pointsPerSubcell, target);
              });
        }
      }
    }
  }
  return outputs;
}

/// Sets up the output of the wave field. Every cell is written in its configuration; with several
/// configurations in the run, a quantity is written once for the cells of all of them.
void setupWaveFieldOutput(
    seissol::SeisSol& seissolInstance,
    const seissol::initializer::parameters::WaveFieldOutputParameters& parameters,
    const std::string& name) {
  const auto& seissolParams = seissolInstance.parameters();
  auto& memoryManager = seissolInstance.memoryManager();
  const auto* meshReader = &seissolInstance.meshReader();

  const auto orderIO = parameters.vtkorder;
  const auto order = static_cast<uint32_t>(std::max(0, orderIO));

  std::vector<std::size_t> celllist;
  celllist.reserve(meshReader->getElements().size());
  if (parameters.bounds.enabled || !parameters.groups.empty()) {
    const auto& vertexArray = meshReader->getVertices();
    for (std::size_t i = 0; i < meshReader->getElements().size(); ++i) {
      const auto& element = meshReader->getElements()[i];
      const auto& vertex0 = vertexArray[element.vertices[0]].coords;
      const auto& vertex1 = vertexArray[element.vertices[1]].coords;
      const auto& vertex2 = vertexArray[element.vertices[2]].coords;
      const auto& vertex3 = vertexArray[element.vertices[3]].coords;
      const bool inGroup = parameters.groups.empty() || parameters.groups.contains(element.group);
      const bool inRegion = !parameters.bounds.enabled ||
                            (parameters.bounds.contains(vertex0[0], vertex0[1], vertex0[2]) ||
                             parameters.bounds.contains(vertex1[0], vertex1[1], vertex1[2]) ||
                             parameters.bounds.contains(vertex2[0], vertex2[1], vertex2[2]) ||
                             parameters.bounds.contains(vertex3[0], vertex3[1], vertex3[2]));
      if (inGroup && inRegion) {
        celllist.push_back(i);
      }
    }
  } else {
    for (std::size_t i = 0; i < meshReader->getElements().size(); ++i) {
      celllist.push_back(i);
    }
  }

  WaveFieldOutputSetup setup;
  setup.order = order;
  setup.dataOrder = order > 0 ? order : 0;
  const auto trueOrder = order > 0 ? order : 1;
  const auto trueBase = io::instance::geometry::pointsTetrahedron(trueOrder);
  setup.dataBase = io::instance::geometry::pointsTetrahedron(setup.dataOrder);

  setup.subcells = io::instance::geometry::unrefined<3>();

  if (parameters.refinement == seissol::initializer::parameters::VolumeRefinement::Refine4) {
    setup.subcells = io::instance::geometry::subdivideMaps(
        setup.subcells, io::instance::geometry::TetrahedronRefine4);
  }
  if (parameters.refinement == seissol::initializer::parameters::VolumeRefinement::Refine8) {
    setup.subcells = io::instance::geometry::subdivideMaps(
        setup.subcells, io::instance::geometry::TetrahedronRefine8);
  }
  if (parameters.refinement == seissol::initializer::parameters::VolumeRefinement::Refine32) {
    // the edge division has to come first; the legacy refinement::DivideTetrahedronBy32
    // subdivided by 8 and then split each of those subcells by its center point, which (as
    // subdivideMaps enumerates input-major) is the ordering 4*i + j the output cells had
    setup.subcells = io::instance::geometry::subdivideMaps(
        setup.subcells, io::instance::geometry::TetrahedronRefine8);
    setup.subcells = io::instance::geometry::subdivideMaps(
        setup.subcells, io::instance::geometry::TetrahedronRefine4);
  }

  const auto truePoints = io::instance::geometry::applyMaps(setup.subcells, trueBase);

  setup.projectionTarget =
      parameters.projection == seissol::initializer::parameters::ProjectionMethod::L2
          ? projection::Target::Project
          : projection::Target::Interpolate;

  const auto format = orderIO < 0 ? io::instance::geometry::WriterFormat::Xdmf
                                  : io::instance::geometry::WriterFormat::Vtk;

  const auto config = io::instance::geometry::WriterConfig{
      .order = order,
      .format = format,
      .backend = seissolParams.output.xdmfWriterBackend ==
                         seissol::initializer::parameters::XdmfBackend::Posix
                     ? io::instance::geometry::WriterBackend::Binary
                     : io::instance::geometry::WriterBackend::Hdf5,
      .time = io::instance::geometry::supportedWriterGroup(
          writerGroupOf(parameters.timeSeries), format, "wavefield"),
      .compress = seissolParams.output.hdfcompress};

  const auto cells = std::make_shared<const std::vector<std::size_t>>(std::move(celllist));
  setup.cells = cells;

  io::writer::ScheduledWriter schedWriter;
  schedWriter.name = name;
  schedWriter.interval = parameters.interval;
  // the unknowns of the cells, their time integrals and their plastic strain
  schedWriter.hostData = {};
  schedWriter.hostData.add<LTS::Dofs>(memoryManager.ltsStorage())
      .add<LTS::Integrals>(memoryManager.ltsStorage())
      .add<LTSPlasticity::PStrain>(memoryManager.ltsStorage());
  auto writer = io::instance::geometry::GeometryWriter(
      name,
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
  const auto configs = seissolParams.runConfigs.configs();
  auto& ltsStorage = memoryManager.ltsStorage();
  const auto& backmap = memoryManager.backmap();
  std::vector<ConfigId> configOfCell(cells->size());
  for (std::size_t index = 0; index < cells->size(); ++index) {
    configOfCell[index] =
        ltsStorage.lookup<LTS::SecondaryInformation>(backmap.get(cells->at(index), 0)).configId;
  }
  const auto configOf = std::make_shared<const std::vector<ConfigId>>(std::move(configOfCell));
  addConfigData(writer, configs, configOf);

  std::vector<NamedOutputs> outputs(builtConfigCount());
  for (const auto runConfig : configs) {
    dispatchConfig(runConfig, [&](auto cfg) {
      outputs[runConfig] = waveFieldOutputsOf<decltype(cfg)>(seissolInstance, setup, parameters);
    });
  }
  addOutputsByName(writer, outputs, configs, configOf, setup.dataBase.size());

  schedWriter.planWrite = writer.makeWriter();
  seissolInstance.outputManager().addOutput(schedWriter);
}

/// What the free surface output of all configurations shares.
struct SurfaceOutputSetup {
  std::vector<io::instance::geometry::Subcell<2>> subcells;
  // the output points of a subcell, in its reference coordinates
  std::vector<std::array<double, 2>> dataBase;
  std::uint32_t order{};
  std::uint32_t dataOrder{};
  projection::Target projectionTarget{};
};

/// The outputs of the free surface of the faces of cells of the configuration `Cfg`.
template <typename Cfg>
NamedOutputs
    surfaceOutputsOf(seissol::SeisSol& seissolInstance,
                     const SurfaceOutputSetup& setup,
                     const seissol::initializer::parameters::FreeSurfaceOutputParameters& surface) {
  using real = Real<Cfg>;
  using MaterialT = model::MaterialOf<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();
  // the projection matrices of the configuration
  constexpr std::size_t MaxProjectionOrder = Cfg::ConvergenceOrder;
  auto& memoryManager = seissolInstance.memoryManager();
  const auto* ltsStorage = &memoryManager.ltsStorage();
  const auto* backmap = &memoryManager.backmap();
  auto* freeSurfaceIntegrator = &seissolInstance.freeSurfaceIntegrator();
  const auto* surfaceMeshIds = freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::MeshId>();
  const auto* surfaceMeshSides = freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::Side>();

  constexpr auto QDofSizePadded =
      tensor::Q<Cfg>::Size / tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  constexpr auto FaceDisplacementPadded =
      tensor::faceDisplacement<Cfg>::Size /
      tensor::faceDisplacement<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];

  const auto simulationCount = memoryManager.simulationSettings().simulationCount;
  const auto namewrap = [simulationCount](const std::string& name, std::size_t sim) {
    if (simulationCount > 1) {
      return name + "-" + std::to_string(sim + 1);
    } else {
      return name;
    }
  };

  const auto order = setup.order;
  const auto pointsPerSubcell = setup.dataBase.size();

  // volume basis -> face points, one table per side of the reference tetrahedron
  std::array<std::shared_ptr<projection::Table<2, 3, real>>, Cell::NumFaces> proj{};
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    const auto embedding = faceEmbedding(f);
    std::vector<seissol::numerical::AffineMap<2, 3>> embedded;
    embedded.reserve(setup.subcells.size());
    for (const auto& subcell : setup.subcells) {
      embedded.emplace_back(embedding.compose(subcell));
    }

    projection::Spec spec;
    spec.target = setup.projectionTarget;
    proj[f] = std::make_shared<projection::Table<2, 3, real>>(
        embedded,
        setup.dataBase,
        setup.dataOrder,
        projectionStride<Cfg, tensor::collvf<Cfg>>(order),
        spec,
        MinProjectionOrder,
        MaxProjectionOrder);
  }

  // face nodes -> face points (the nodal-to-modal transform is folded in)
  projection::Spec faceSpec;
  faceSpec.source = projection::Source::Nodal;
  faceSpec.target = setup.projectionTarget;
  // the face displacement is stored at the nodes2D points, which are always warp&blend
  faceSpec.nodalSet = projection::NodalSet::WarpBlend;
  const auto projf = std::make_shared<projection::Table<2, 2, real>>(
      setup.subcells,
      setup.dataBase,
      setup.dataOrder,
      projectionStride<Cfg, tensor::collnf<Cfg>>(order),
      faceSpec,
      MinProjectionOrder,
      MaxProjectionOrder);

  constexpr std::size_t MaxVtk2dPoints = tensor::vtk2d<Cfg>::Shape
      [(sizeof(tensor::vtk2d<Cfg>::Shape) / sizeof(tensor::vtk2d<Cfg>::Shape[0])) - 1][1];

  const std::vector<std::string> quantityLabelsDisplacement = {"u1", "u2", "u3"};

  NamedOutputs outputs;
  for (std::size_t sim = 0; sim < simulationCount; ++sim) {
    // the copy of the cells which holds the simulation, and the fused simulation there
    const auto copy = sim / Cfg::NumSimulations;
    const auto fused = sim % Cfg::NumSimulations;

    for (std::size_t quantity = 0; quantity < MaterialT::Quantities.size(); ++quantity) {
      if (surface.outputMask[Variant][quantity]) {
        outputs.emplace_back(
            namewrap(MaterialT::Quantities[quantity], sim),
            [=](double* target, std::size_t index, std::size_t subcell) {
              const auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[index]];
              const auto side = surfaceMeshSides[freeSurfaceIntegrator->backmap[index]];
              const auto position = backmap->get(meshId, copy);
              const auto* dofsAllQuantities = ltsStorage->lookup<LTS::Dofs>(Cfg(), position);
              const auto* dofsSingleQuantity = dofsAllQuantities + QDofSizePadded * quantity;
              runtime::kernel::projectBasisToVtkFaceFromVolume vtkproj{};
              memory::AlignedArray<real, Cfg::NumSimulations> simselect{};
              alignas(Alignment) std::array<real, MaxVtk2dPoints> alignedTarget{};
              simselect[fused] = 1;
              vtkproj.simselect = runtime::init::simselect::view(Variant, simselect.data());
              vtkproj.qb = runtime::init::qb::view(Variant, dofsSingleQuantity);
              vtkproj.xf(order) = runtime::init::xf::view(Variant, order, alignedTarget.data());
              vtkproj.collvf(Cfg::ConvergenceOrder, order) =
                  runtime::init::collvf::view(Variant,
                                              Cfg::ConvergenceOrder,
                                              order,
                                              (*proj[side])(subcell, Cfg::ConvergenceOrder));
              vtkproj.execute(Variant, order);
              std::copy_n(alignedTarget.data(), pointsPerSubcell, target);
            });
      }
    }
    for (std::size_t quantity = 0; quantity < quantityLabelsDisplacement.size(); ++quantity) {
      outputs.emplace_back(
          namewrap(quantityLabelsDisplacement[quantity], sim),
          [=](double* target, std::size_t index, std::size_t subcell) {
            const auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[index]];
            const auto side = surfaceMeshSides[freeSurfaceIntegrator->backmap[index]];
            const auto position = backmap->get(meshId, copy);
            const auto& faceDisplacements =
                ltsStorage->lookup<LTS::FaceDisplacements>(Cfg(), position);
            const auto* faceDisplacementVariable =
                faceDisplacements[side] + FaceDisplacementPadded * quantity;
            runtime::kernel::projectNodalToVtkFace vtkproj{};
            memory::AlignedArray<real, Cfg::NumSimulations> simselect{};
            alignas(Alignment) std::array<real, MaxVtk2dPoints> alignedTarget{};
            simselect[fused] = 1;
            vtkproj.simselect = runtime::init::simselect::view(Variant, simselect.data());
            vtkproj.pn = runtime::init::pn::view(Variant, faceDisplacementVariable);
            vtkproj.xf(order) = runtime::init::xf::view(Variant, order, alignedTarget.data());
            vtkproj.collnf(Cfg::ConvergenceOrder, order) = runtime::init::collnf::view(
                Variant, Cfg::ConvergenceOrder, order, (*projf)(subcell, Cfg::ConvergenceOrder));
            vtkproj.execute(Variant, order);
            std::copy_n(alignedTarget.data(), pointsPerSubcell, target);
          });
    }
  }
  return outputs;
}

/// Sets up the output of the free surface. Every face is written in the configuration of its
/// cell; with several configurations in the run, a quantity is written once for the faces of all
/// of them.
void setupSurfaceOutput(
    seissol::SeisSol& seissolInstance,
    const seissol::initializer::parameters::FreeSurfaceOutputParameters& parameters,
    const std::string& name) {
  const auto& seissolParams = seissolInstance.parameters();
  auto& memoryManager = seissolInstance.memoryManager();

  const auto orderIO = parameters.vtkorder;
  const auto order = static_cast<std::uint32_t>(std::max(0, orderIO));

  auto* freeSurfaceIntegrator = &seissolInstance.freeSurfaceIntegrator();
  const auto* meshReader = &seissolInstance.meshReader();
  io::writer::ScheduledWriter schedWriter;
  schedWriter.name = name;
  schedWriter.interval = parameters.interval;
  // the unknowns of the cells at the surface, and its displacements
  schedWriter.hostData = {};
  schedWriter.hostData.add<LTS::Dofs>(memoryManager.ltsStorage())
      .add<SurfaceLTS::DisplacementDofs>(memoryManager.surfaceStorage());
  const auto* surfaceMeshIds = freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::MeshId>();
  const auto* surfaceMeshSides = freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::Side>();
  const auto* surfaceLocationFlag =
      freeSurfaceIntegrator->surfaceStorage->var<SurfaceLTS::LocationFlag>();

  SurfaceOutputSetup setup;
  setup.order = order;
  const auto trueOrder = order > 0 ? order : 1;
  setup.dataOrder = order > 0 ? order : 0;
  const auto trueBase = io::instance::geometry::pointsTriangle(trueOrder);
  setup.dataBase = io::instance::geometry::pointsTriangle(setup.dataOrder);

  setup.subcells = io::instance::geometry::unrefined<2>();

  for (std::size_t i = 0; i < parameters.refinement; ++i) {
    setup.subcells = io::instance::geometry::subdivideMaps(setup.subcells,
                                                           io::instance::geometry::TriangleRefine4);
  }

  const auto truePoints = io::instance::geometry::applyMaps(setup.subcells, trueBase);

  const auto format = orderIO < 0 ? io::instance::geometry::WriterFormat::Xdmf
                                  : io::instance::geometry::WriterFormat::Vtk;

  const auto config = io::instance::geometry::WriterConfig{
      .order = order,
      .format = format,
      .backend = seissolParams.output.xdmfWriterBackend ==
                         seissol::initializer::parameters::XdmfBackend::Posix
                     ? io::instance::geometry::WriterBackend::Binary
                     : io::instance::geometry::WriterBackend::Hdf5,
      .time = io::instance::geometry::supportedWriterGroup(
          writerGroupOf(parameters.timeSeries), format, "free surface"),
      .compress = seissolParams.output.hdfcompress};

  auto writer = io::instance::geometry::GeometryWriter(
      name,
      freeSurfaceIntegrator->backmap.size(),
      io::instance::geometry::Shape::Triangle,
      config,
      setup.subcells.size(),

      [=](double* target, std::size_t index, std::size_t subcell) {
        const auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[index]];
        const auto side = surfaceMeshSides[freeSurfaceIntegrator->backmap[index]];
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

  setup.projectionTarget =
      parameters.projection == seissol::initializer::parameters::ProjectionMethod::L2
          ? projection::Target::Project
          : projection::Target::Interpolate;

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
        target[0] = surfaceLocationFlag[freeSurfaceIntegrator->backmap[index]];
      });

  writer.addCellData<std::size_t>(
      "global-id", {}, true, [=](std::size_t* target, std::size_t index, std::size_t /*subcell*/) {
        const auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[index]];
        const auto side = surfaceMeshSides[freeSurfaceIntegrator->backmap[index]];
        target[0] = meshReader->getElements()[meshId].globalId * 4 + side;
      });

  // the faces of every configuration of the run
  const auto configs = seissolParams.runConfigs.configs();
  auto& ltsStorage = memoryManager.ltsStorage();
  const auto& backmap = memoryManager.backmap();
  std::vector<ConfigId> configOfFace(freeSurfaceIntegrator->backmap.size());
  for (std::size_t index = 0; index < configOfFace.size(); ++index) {
    const auto meshId = surfaceMeshIds[freeSurfaceIntegrator->backmap[index]];
    configOfFace[index] =
        ltsStorage.lookup<LTS::SecondaryInformation>(backmap.get(meshId, 0)).configId;
  }
  const auto configOf = std::make_shared<const std::vector<ConfigId>>(std::move(configOfFace));
  addConfigData(writer, configs, configOf);

  std::vector<NamedOutputs> outputs(builtConfigCount());
  for (const auto runConfig : configs) {
    dispatchConfig(runConfig, [&](auto cfg) {
      outputs[runConfig] = surfaceOutputsOf<decltype(cfg)>(seissolInstance, setup, parameters);
    });
  }
  addOutputsByName(writer, outputs, configs, configOf, setup.dataBase.size());

  schedWriter.planWrite = writer.makeWriter();
  seissolInstance.outputManager().addOutput(schedWriter);
}

void setupOutput(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.parameters();
  auto& memoryManager = seissolInstance.memoryManager();

  for (const auto& entry : seissolParams.output.outputs) {
    if (const auto* waveField =
            std::get_if<seissol::initializer::parameters::WaveFieldOutputParameters>(
                &entry.settings)) {
      setupWaveFieldOutput(seissolInstance, *waveField, entry.name);
    }
    if (const auto* surface =
            std::get_if<seissol::initializer::parameters::FreeSurfaceOutputParameters>(
                &entry.settings)) {
      setupSurfaceOutput(seissolInstance, *surface, entry.name);
    }
    if (const auto* receivers =
            std::get_if<seissol::initializer::parameters::ReceiverOutputParameters>(
                &entry.settings)) {
      // the receivers of the output, which the cell clusters sample
      auto& receiverWriter = seissolInstance.addReceiverWriter();
      receiverWriter.init(seissolInstance.outputPrefix(),
                          entry.name,
                          seissolParams.timeStepping.endTime,
                          *receivers);
      receiverWriter.addPoints(
          seissolInstance.meshReader(),
          memoryManager.backmap(),
          kernels::ReceiverCells{.storage = memoryManager.ltsStorage(),
                                 .globalData = memoryManager.globalDataStore(),
                                 .settings = memoryManager.simulationSettings(),
                                 .gravitationalAcceleration =
                                     seissolInstance.gravitationSetup().acceleration});
      seissolInstance.timeManager().addReceivers(receiverWriter);
    }
  }

  if (const auto& energy = seissolParams.output.energy; energy.has_value()) {
    auto& energyOutput = seissolInstance.energyOutput();

    energyOutput.init(memoryManager.drStorage(),
                      seissolInstance.meshReader(),
                      memoryManager.ltsStorage(),
                      seissolParams.plasticity.enabled,
                      seissolInstance.outputPrefix(),
                      energy.value());
  }

  seissolInstance.flopCounter().init(seissolInstance.outputPrefix(), seissolInstance.runId());

  seissolInstance.epochSynchronization().init(seissolInstance.modules(),
                                              seissolParams.output.epochInterval);

  seissolInstance.analysisWriter().init(&seissolInstance.meshReader(),
                                        seissolInstance.outputPrefix(),
                                        seissolParams.initialization,
                                        memoryManager.initialConditions(),
                                        memoryManager.ltsStorage(),
                                        memoryManager.simulationSettings().simulationCount);
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

void enableFreeSurfaceOutput(const seissol::SeisSol& seissolInstance) {}

} // namespace

void seissol::initializer::initprocedure::createOutputDirectory(seissol::SeisSol& seissolInstance) {
  const auto& outputParams = seissolInstance.parameters().output;

  // an empty directory denotes the working directory
  auto directory = outputParams.outputPath.empty() ? std::filesystem::path(".")
                                                   : std::filesystem::path(outputParams.outputPath);

  if (outputParams.outputDirPerSim) {
    // the ID is the same on all ranks, and so is the directory name; unlike a time stamp, it
    // differs between runs that start at the same time, e.g. those of a job array
    directory /= outputParams.outputName + "-" + seissolInstance.runId();
  }

  createDirectoryGlobally(directory);

  logInfo() << "Output directory:" << directory.string();
  logInfo() << "Output file prefix:" << outputParams.outputName;

  // TODO: (after #1180) properly split prefix into dir and simulation name
  seissolInstance.setOutputPrefix((directory / outputParams.outputName).string());
}

void seissol::initializer::initprocedure::initIO(seissol::SeisSol& seissolInstance) {
  logInfo() << "Begin init output.";

  // recorded while reading the mesh and laying out the clusters, before the writers were set up
  seissolInstance.miniSeisSolWriter().write(
      seissolInstance.outputPrefix(), seissolInstance.runId(), seissolInstance.systemInfo());
  seissolInstance.timeManager().writeClustering();

  enableFreeSurfaceOutput(seissolInstance);
  initFaultOutputManager(seissolInstance);
  setupCheckpointing(seissolInstance);
  setupOutput(seissolInstance);
  logInfo() << "End init output.";
}
