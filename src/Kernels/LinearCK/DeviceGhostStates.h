// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_DEVICEGHOSTSTATES_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_DEVICEGHOSTSTATES_H_

#ifdef ACL_DEVICE

#include "Common/Constants.h"
#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Expr/Backend.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/Typedefs.h"
#include "Kernels/DeviceAux/TimePointsAux.h"
#include "Kernels/GhostStateBatch.h"
#include "Kernels/NonlinearDirichletScratch.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/MemoryAllocator.h"
#include "Parallel/Runtime/Stream.h"
#include "Physics/GhostState.h"
#include "Physics/NonlinearDirichlet.h"

#include <Device/device.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <utils/logger.h>
#include <vector>
#include <yateto.h>

namespace seissol::kernels::solver::linearck {

/// The ghost states of the boundaries a script gives -- the analytical one and the nonlinear
/// Dirichlet one -- for the faces of the cells of one layer, on the device: for the analytical
/// boundary, one kernel of the script for the faces of all sides of the cells
/// (kernels::GhostStateBatch); for the nonlinear Dirichlet boundary, one per side, and the kernels
/// around it, which take the state inside from the time derivatives to the nodes of the faces and
/// the ghost state back to global coordinates. Both write the ghost state of a face into its slot
/// of LTS::AnalyticScratch, which the nodal flux of its side takes: the recorder gives every face
/// a slot of its own, the analytical faces of all sides first, side by side.
///
/// The kernels read the time of the step from the clock of the cluster on the device, through the
/// times of the points of the quadrature that beginStep() computes there: the work of a step
/// replayed from a graph (a super-timestep, say) then evaluates the boundaries at the time of the
/// step it runs in. Everything else a kernel takes by value depends on the width of the step only.
/// The kernels are made by prepare(), before the time stepping, so that no recording of a graph
/// makes them.
template <typename Cfg>
class DeviceGhostStates {
  public:
  using real = Real<Cfg>;

  static constexpr std::size_t NodeCount = tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg>];

  explicit DeviceGhostStates(const GlobalData<Cfg>* global) : global_(global) {}

  /// Makes the kernels of the boundaries of `layer`: of the analytical one if `state` is given, of
  /// the nonlinear Dirichlet one if `condition` is.
  void prepare(const physics::AnalyticalGhostState* state,
               const physics::NonlinearDirichlet* condition,
               recording::ConditionalPointersToRealsTable& dataTable,
               recording::ConditionalIndicesTable& indicesTable,
               LTS::Layer& layer,
               void* stream) {
    using namespace seissol::recording;
    if (state != nullptr) {
      makeAnalytical(*state, dataTable, indicesTable, layer, stream);
    }
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const ConditionalKey nonlinearKey(
          *KernelNames::BoundaryConditions, *ComputationKind::NonlinearDirichlet, face);
      if (condition != nullptr && indicesTable.contains(nonlinearKey)) {
        makeNonlinearDirichlet(*condition, dataTable, indicesTable, layer, face, stream);
      }
    }
  }

  /// The times of the points of the quadrature of the step that starts at `*clock`, on the device
  /// (solver::ClusterClock::device()), and of the width `timeStepWidth`, for the kernels that
  /// follow on `stream`: once per step, before the boundaries of its faces.
  void beginStep(const double* clock, double timeStepWidth, void* stream) {
    if (timePointCount_ == 0) {
      return;
    }
    kernels::timepoints::TimeOffsets offsets;
    offsets.count = timePointCount_;
    offsetsOf_->timeOffsets(timeStepWidth, offsets.values);
    kernels::timepoints::timePoints(clock, offsets, times_.data(), stream);
  }

  /// The ghost state of the analytical boundary at the faces of all sides of the cells of
  /// `layer`, integrated over the step of the width `timeStepWidth` that beginStep() began: per
  /// simulation, the integral of the script of its field (`state`), with one kernel per program.
  /// What the scripts do not give is zero, which the kernels write as well.
  void analytical(const physics::AnalyticalGhostState& state,
                  recording::ConditionalPointersToRealsTable& dataTable,
                  recording::ConditionalIndicesTable& indicesTable,
                  LTS::Layer& layer,
                  double timeStepWidth,
                  void* stream) {
    auto& all = makeAnalytical(state, dataTable, indicesTable, layer, stream);
    if (all.faces == 0) {
      return;
    }
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      all.batches[all.batchOf[s]]->evaluate(
          timeStepWidth, s, static_cast<double>(state.scriptSimulation(s)), stream);
    }
  }

  /// The ghost state of the nonlinear Dirichlet boundary at the faces `face` of the cells of
  /// `layer`, integrated over the step of the width `timeStepWidth` that beginStep() began: the
  /// integral of the state inside (projected from the integrated state), and where the script gives
  /// a quantity, the integral of the script (condition.deviceIntegral()), from the state inside at
  /// the time points of the step if it reads it.
  void nonlinearDirichlet(const physics::NonlinearDirichlet& condition,
                          recording::ConditionalPointersToRealsTable& dataTable,
                          recording::ConditionalIndicesTable& indicesTable,
                          LTS::Layer& layer,
                          std::size_t face,
                          double timeStepWidth,
                          parallel::runtime::StreamRuntime& runtime) {
#ifdef SEISSOL_KERNELS_LINEARCK
    if constexpr (Cfg::Solver == SolverType::LinearCK) {
      using namespace seissol::recording;
      using Scratch = kernels::NonlinearDirichletScratch<Cfg>;
      const ConditionalKey key(
          *KernelNames::BoundaryConditions, *ComputationKind::NonlinearDirichlet, face);
      auto& entry = dataTable[key];
      auto* cells = indicesTable[key].get<unsigned>(inner_keys::Indices::Id::Cells);
      const std::size_t faces = cells->getSize();
      const auto* integral = condition.deviceIntegral();
      if (integral == nullptr) {
        logError() << "The nonlinear Dirichlet boundary has no program for the device.";
        return;
      }
      const bool aligned = condition.faceAligned();
      void* stream = runtime.stream();

      auto& side = makeNonlinearDirichlet(condition, dataTable, indicesTable, layer, face, stream);
      // the slots of the side
      auto* ghost = side.ghost;

      auto** integrated = entry.get<real*>(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr();
      auto** blocks =
          entry.get<real*>(inner_keys::Wp::Id::NonlinearDirichletScratch)->getDeviceDataPtr();
      auto** ghosts = entry.get<real*>(inner_keys::Wp::Id::Analytical)->getDeviceDataPtr();
      auto** derivatives = entry.get<real*>(inner_keys::Wp::Id::Derivatives)->getDeviceDataPtr();
      auto** rotations = entry.get<real*>(inner_keys::Wp::Id::T)->getDeviceDataPtr();
      auto** inverseRotations = entry.get<real*>(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();

      // the integral of the state inside: the ghost state of what the script does not give
      if (aligned) {
        // the columns beyond the quantities of the Riemann problem stay zero
        device::DeviceInstance::instance().algorithms().fillArray(
            ghost, static_cast<real>(0), faces * tensor::INodal<Cfg>::size(), stream);
        kernel::gpu_projectToFaceNodesAligned<Cfg> project;
        project.bindGlobals(*global_);
        project.I = const_cast<const real**>(integrated);
        project.Tinv = const_cast<const real**>(inverseRotations);
        project.INodalAligned = blocks;
        project.extraOffset_INodalAligned = Scratch::Aligned;
        launch(project, faces, runtime, face);
      } else {
        kernel::gpu_projectToFaceNodes<Cfg> project;
        project.bindGlobals(*global_);
        project.I = const_cast<const real**>(integrated);
        project.INodal = ghosts;
        launch(project, faces, runtime, face);
      }

      // the state inside at the time points of the step, where the script reads it: the Taylor
      // series at the offsets of the points into the step, which depend on its width alone
      if (integral->innerCount() > 0) {
        std::array<double, Cfg::ConvergenceOrder> offsets{};
        integral->timeOffsets(timeStepWidth, offsets.data());
        if (aligned) {
          kernel::gpu_projectToFaceNodesAlignedAtTimes<Cfg> project;
          project.bindGlobals(*global_);
          project.Tinv = const_cast<const real**>(inverseRotations);
          std::size_t offset = Scratch::Inner;
          for (std::size_t c = 0; c < Cfg::ConvergenceOrder; ++c) {
            project.INodalAlignedAt(c) = blocks;
            project.extraOffset_INodalAlignedAt(c) = static_cast<int>(offset);
            checkTimePoint(offset - Scratch::Inner, c * tensor::INodalAligned<Cfg>::size());
            offset += tensor::INodalAlignedAt<Cfg>::size(c);
          }
          stateAt(project, offsets, derivatives, blocks);
          launch(project, faces, runtime, face);
        } else {
          kernel::gpu_projectToFaceNodesAtTimes<Cfg> project;
          project.bindGlobals(*global_);
          std::size_t offset = Scratch::Inner;
          for (std::size_t c = 0; c < Cfg::ConvergenceOrder; ++c) {
            project.INodalAt(c) = blocks;
            project.extraOffset_INodalAt(c) = static_cast<int>(offset);
            checkTimePoint(offset - Scratch::Inner, c * tensor::INodal<Cfg>::size());
            offset += tensor::INodalAt<Cfg>::size(c);
          }
          stateAt(project, offsets, derivatives, blocks);
          launch(project, faces, runtime, face);
        }
      }

      // the script, where it gives a quantity
      for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
        if (aligned) {
          side.aligned->evaluate(timeStepWidth, s, static_cast<double>(s), stream);
        } else {
          side.batches.front()->evaluate(timeStepWidth, s, static_cast<double>(s), stream);
        }
      }

      if (aligned) {
        kernel::gpu_rotateNodalToGlobal<Cfg> rotate;
        rotate.bindGlobals(*global_);
        rotate.INodalAligned = const_cast<const real**>(blocks);
        rotate.extraOffset_INodalAligned = Scratch::Aligned;
        rotate.T = const_cast<const real**>(rotations);
        rotate.INodal = ghosts;
        launchOnce(rotate, faces, runtime);
      }
    }
#else
    static_cast<void>(condition);
    static_cast<void>(dataTable);
    static_cast<void>(indicesTable);
    static_cast<void>(layer);
    static_cast<void>(face);
    static_cast<void>(timeStepWidth);
    static_cast<void>(runtime);
    logError() << "The nonlinear Dirichlet boundary needs the Taylor series of the state.";
#endif
  }

  private:
  using GlobalBatch = kernels::GhostStateBatch<Cfg, false>;
  using AlignedBatch = kernels::GhostStateBatch<Cfg, true>;

  struct Side;

  struct Analytical;

  /// The kernels of the analytical boundary at the faces of all sides, made if they are not yet.
  Analytical& makeAnalytical(const physics::AnalyticalGhostState& state,
                             recording::ConditionalPointersToRealsTable& dataTable,
                             recording::ConditionalIndicesTable& indicesTable,
                             LTS::Layer& layer,
                             void* stream) {
    using namespace seissol::recording;
    auto& all = analytical_;
    if (all.made) {
      return all;
    }
    all.made = true;
    auto* ghost = reinterpret_cast<real*>(
        layer.var<LTS::AnalyticScratch>(Cfg(), initializer::AllocationPlace::Device));
    // the faces of all sides, side by side, in the slots the recorder gave them
    std::vector<unsigned> cells;
    std::vector<real> nodes;
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const ConditionalKey key(
          *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
      if (!indicesTable.contains(key)) {
        continue;
      }
      const auto& sideCells =
          indicesTable[key].get<unsigned>(inner_keys::Indices::Id::Cells)->getHostData();
      const auto& slots = dataTable[key].get<real*>(inner_keys::Wp::Id::Analytical)->getHostData();
      for (std::size_t e = 0; e < sideCells.size(); ++e) {
        if (slots[e] != ghost + (cells.size() + e) * tensor::INodal<Cfg>::size()) {
          logError() << "The ghost state of an analytical face does not lie in the slot that "
                        "follows the one of the face before it.";
        }
      }
      const auto sideNodes = hostNodesOf(layer, sideCells, face);
      cells.insert(cells.end(), sideCells.begin(), sideCells.end());
      nodes.insert(nodes.end(), sideNodes.begin(), sideNodes.end());
    }
    all.faces = cells.size();
    if (all.faces == 0) {
      return all;
    }
    all.nodes = memory::MemkindArray<real>(nodes, memory::Memkind::DeviceGlobalMemory);
    all.cells = memory::MemkindArray<unsigned>(cells, memory::Memkind::DeviceGlobalMemory);

    typename GlobalBatch::Data data;
    data.faces = all.faces;
    data.nodes = all.nodes.data();
    data.parameters = reinterpret_cast<const real*>(
        layer.var<LTS::MaterialParameters>(Cfg(), initializer::AllocationPlace::Device));
    data.cells = all.cells.data();
    data.ghost = ghost;
    data.ghostStride = tensor::INodal<Cfg>::size();
    // what the scripts do not give is zero
    data.zeroOthers = true;
    // one kernel per program: fused simulations of one script share theirs
    std::vector<const physics::GhostStateIntegral*> integrals;
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      const auto* integral = &state.integral(s);
      const auto found = std::find(integrals.begin(), integrals.end(), integral);
      all.batchOf.push_back(static_cast<std::size_t>(found - integrals.begin()));
      if (found == integrals.end()) {
        integrals.push_back(integral);
        data.times = timesFor(*integral);
        all.batches.push_back(
            std::make_unique<GlobalBatch>(*integral, data, optionsFor(stream, quiet_)));
        checkKind(all.batches.back()->kind(), "the analytical boundary");
        quiet_ = true;
      }
    }
    return all;
  }

  /// The kernels of the nonlinear Dirichlet boundary at the faces `face`, made if they are not
  /// yet.
  Side& makeNonlinearDirichlet(const physics::NonlinearDirichlet& condition,
                               recording::ConditionalPointersToRealsTable& dataTable,
                               recording::ConditionalIndicesTable& indicesTable,
                               LTS::Layer& layer,
                               std::size_t face,
                               void* stream) {
    using namespace seissol::recording;
    using Scratch = kernels::NonlinearDirichletScratch<Cfg>;
    auto& side = nonlinear_[face];
    const auto* integral = condition.deviceIntegral();
    if (side.made || integral == nullptr) {
      return side;
    }
    side.made = true;
    const ConditionalKey key(
        *KernelNames::BoundaryConditions, *ComputationKind::NonlinearDirichlet, face);
    auto* cells = indicesTable[key].get<unsigned>(inner_keys::Indices::Id::Cells);
    // the slots of the ghost states of the side, which the recorder gave them one after the other
    side.ghost = dataTable[key].get<real*>(inner_keys::Wp::Id::Analytical)->getHostData().front();
    auto* scratch = reinterpret_cast<real*>(
        layer.var<LTS::NonlinearDirichletScratch>(Cfg(), initializer::AllocationPlace::Device));
    side.nodes = nodesOf(layer, cells->getHostData(), face);
    auto data = dataOf(layer, *cells, side.nodes);
    data.inner = scratch + Scratch::Inner;
    data.innerStride = Scratch::PerFace;
    data.times = timesFor(*integral);
    if (condition.faceAligned()) {
      typename AlignedBatch::Data alignedData;
      alignedData.faces = data.faces;
      alignedData.nodes = data.nodes;
      alignedData.parameters = data.parameters;
      alignedData.cells = data.cells;
      alignedData.inner = data.inner;
      alignedData.innerStride = data.innerStride;
      alignedData.times = data.times;
      alignedData.ghost = scratch + Scratch::Aligned;
      alignedData.ghostStride = Scratch::PerFace;
      side.aligned =
          std::make_unique<AlignedBatch>(*integral, alignedData, optionsFor(stream, quiet_));
      checkKind(side.aligned->kind(), "the nonlinear Dirichlet boundary");
    } else {
      data.ghost = side.ghost;
      data.ghostStride = tensor::INodal<Cfg>::size();
      side.batches.push_back(
          std::make_unique<GlobalBatch>(*integral, data, optionsFor(stream, quiet_)));
      checkKind(side.batches.back()->kind(), "the nonlinear Dirichlet boundary");
    }
    quiet_ = true;
    return side;
  }

  /// The times on the device that the kernels of `integral` read (see beginStep()). All integrals
  /// of the boundaries of a layer integrate over the same points of the step, which they share.
  const double* timesFor(const physics::GhostStateIntegral& integral) {
    const std::size_t count = integral.timePointCount();
    if (count > kernels::timepoints::MaxTimeOffsets) {
      logError() << "The ghost state of a boundary integrates over" << count
                 << "points of a step; the device takes at most"
                 << kernels::timepoints::MaxTimeOffsets << ".";
    }
    if (timePointCount_ == 0) {
      timePointCount_ = count;
      offsetsOf_ = &integral;
      times_ = memory::MemkindArray<double>(count, memory::Memkind::DeviceGlobalMemory);
    } else {
      std::array<double, kernels::timepoints::MaxTimeOffsets> mine{};
      std::array<double, kernels::timepoints::MaxTimeOffsets> theirs{};
      integral.timeOffsets(1.0, mine.data());
      offsetsOf_->timeOffsets(1.0, theirs.data());
      if (count != timePointCount_ || mine != theirs) {
        logError() << "The ghost states of the boundaries of a layer integrate over other points "
                      "of a step.";
      }
    }
    return times_.data();
  }

  static_assert(std::is_same_v<unsigned, std::uint32_t>,
                "the cells of a batch are the cell index of a column of a program");

  /// The kernels of the analytical boundary, for the faces of all sides of the cells of the layer.
  struct Analytical {
    bool made{false};
    std::size_t faces{0};
    /// the nodes of the faces and the cells they belong to, on the device
    /// (GhostStateBatch::Data::nodes and cells)
    memory::MemkindArray<real> nodes{memory::Memkind::DeviceGlobalMemory};
    memory::MemkindArray<unsigned> cells{memory::Memkind::DeviceGlobalMemory};
    /// one per program, of the fields
    std::vector<std::unique_ptr<GlobalBatch>> batches;
    /// per simulation, its program
    std::vector<std::size_t> batchOf;
  };

  /// The kernels of the nonlinear Dirichlet boundary at the faces of one side of the cells of the
  /// layer.
  struct Side {
    bool made{false};
    /// the first slot of the ghost states of the faces (LTS::AnalyticScratch)
    real* ghost{nullptr};
    /// the nodes of the faces, on the device (GhostStateBatch::Data::nodes)
    memory::MemkindArray<real> nodes{memory::Memkind::DeviceGlobalMemory};
    /// in global coordinates
    std::vector<std::unique_ptr<GlobalBatch>> batches;
    /// in the face-aligned basis
    std::unique_ptr<AlignedBatch> aligned;
  };

  /// The nodes of the faces `face` of the cells `cells` of `layer`, on the host.
  static std::vector<real>
      hostNodesOf(LTS::Layer& layer, const std::vector<unsigned>& cells, std::size_t face) {
    const auto* mappings = layer.var<LTS::BoundaryMapping>(Cfg());
    std::vector<real> nodes(cells.size() * NodeCount * 3);
    for (std::size_t e = 0; e < cells.size(); ++e) {
      const auto* source = mappings[cells[e]][face].nodes;
      std::copy_n(source, NodeCount * 3, nodes.data() + e * NodeCount * 3);
    }
    return nodes;
  }

  /// The same, on the device.
  static memory::MemkindArray<real>
      nodesOf(LTS::Layer& layer, const std::vector<unsigned>& cells, std::size_t face) {
    return memory::MemkindArray<real>(hostNodesOf(layer, cells, face),
                                      memory::Memkind::DeviceGlobalMemory);
  }

  /// What a batch of the faces of `cells` reads, besides the state inside.
  static typename GlobalBatch::Data dataOf(LTS::Layer& layer,
                                           recording::GenericTableEntry<unsigned>& cells,
                                           const memory::MemkindArray<real>& nodes) {
    typename GlobalBatch::Data data;
    data.faces = cells.getSize();
    data.nodes = nodes.data();
    data.parameters = reinterpret_cast<const real*>(
        layer.var<LTS::MaterialParameters>(Cfg(), initializer::AllocationPlace::Device));
    data.cells = cells.getDeviceDataPtr();
    return data;
  }

  static expr::BackendOptions optionsFor(void* stream, bool quiet) {
    expr::BackendOptions options;
    options.preferred = expr::deviceBackendOfBuild();
    options.deviceQueue = stream;
    options.quiet = quiet;
    return options;
  }

  /// A device build evaluates the boundaries on the device, without a fallback to the host.
  static void checkKind(expr::BackendKind kind, const char* boundary) {
    if (kind != expr::deviceBackendOfBuild()) {
      logError() << "The script of" << boundary << "cannot run on the device (see above), and a "
                 << "device build does not evaluate it on the host.";
    }
  }

  /// The time points of the state inside lie one tensor after the other, as the program reads
  /// them.
  static void checkTimePoint(std::size_t offset, std::size_t expected) {
    if (offset != expected) {
      logError() << "The state inside at a time point lies at" << offset << "in the scratch of a "
                 << "face, not at" << expected << ".";
    }
  }

  /// The Taylor series of the state of the cells of the batch at the offsets `offsets` into the
  /// step, as the host evaluates it (kernel::derivativeTaylorExpansion).
  template <typename Kernel>
  static void stateAt(Kernel& project,
                      const std::array<double, Cfg::ConvergenceOrder>& offsets,
                      real** derivatives,
                      real** blocks) {
    std::size_t offset = 0;
    for (std::size_t d = 0; d < Cfg::ConvergenceOrder; ++d) {
      project.dQ(d) = const_cast<const real**>(derivatives);
      project.extraOffset_dQ(d) = static_cast<int>(offset);
      offset += tensor::dQ<Cfg>::size(d);
    }
    for (std::size_t c = 0; c < Cfg::ConvergenceOrder; ++c) {
      const double tau = offsets[c];
      double power = 1;
      for (std::size_t d = 0; d < Cfg::ConvergenceOrder; ++d) {
        project.taylorCoefficient(c * Cfg::ConvergenceOrder + d) = static_cast<real>(power);
        power *= tau / static_cast<double>(d + 1);
      }
    }
    // the state at each time point, in a tensor of its own
    for (std::size_t c = 0; c < Cfg::ConvergenceOrder; ++c) {
      project.faceStateI(c) = blocks;
      project.extraOffset_faceStateI(c) = static_cast<int>(
          kernels::NonlinearDirichletScratch<Cfg>::State + c * tensor::faceStateI<Cfg>::size(c));
    }
  }

  template <typename Kernel>
  static void launch(Kernel& kernel,
                     std::size_t faces,
                     parallel::runtime::StreamRuntime& runtime,
                     std::size_t face) {
    auto tmpMem =
        runtime.memoryHandle<real>((yateto::getMaxTmpMemRequired(kernel) * faces) / sizeof(real));
    kernel.linearAllocator.initialize(tmpMem.get());
    kernel.numElements = faces;
    kernel.streamPtr = runtime.stream();
    kernel.execute(face);
  }

  template <typename Kernel>
  static void
      launchOnce(Kernel& kernel, std::size_t faces, parallel::runtime::StreamRuntime& runtime) {
    auto tmpMem =
        runtime.memoryHandle<real>((yateto::getMaxTmpMemRequired(kernel) * faces) / sizeof(real));
    kernel.linearAllocator.initialize(tmpMem.get());
    kernel.numElements = faces;
    kernel.streamPtr = runtime.stream();
    kernel.execute();
  }

  const GlobalData<Cfg>* global_;
  Analytical analytical_;
  std::array<Side, Cell::NumFaces> nonlinear_;
  /// whether a kernel was made already, which said how it went
  bool quiet_{false};
  /// the times of the points of the quadrature of the current step, on the device (beginStep()),
  /// and the integral whose points they are; none before a kernel is made
  memory::MemkindArray<double> times_{memory::Memkind::DeviceGlobalMemory};
  std::size_t timePointCount_{0};
  const physics::GhostStateIntegral* offsetsOf_{nullptr};
};

} // namespace seissol::kernels::solver::linearck

#endif // ACL_DEVICE

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DEVICEGHOSTSTATES_H_
