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
/// Dirichlet one -- for the faces of the cells of one layer, on the device: per side of the cells,
/// one kernel of the script for all faces of the side (kernels::GhostStateBatch), and for the
/// nonlinear Dirichlet boundary the kernels around it, which take the state inside from the time
/// derivatives to the nodes of the faces and the ghost state back to global coordinates. Both
/// write the ghost state of a face into its slot of LTS::AnalyticScratch, which the nodal flux
/// takes. The kernels of a side are made at its first step.
template <typename Cfg>
class DeviceGhostStates {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  static constexpr std::size_t NodeCount = tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg>];

  explicit DeviceGhostStates(const GlobalData<Cfg>* global) : global_(global) {}

  /// The ghost state of the analytical boundary at the faces `face` of the cells of `layer`,
  /// integrated over the step [time, time + timeStepWidth]: per simulation, the integral of the
  /// script of its field (`state`). What the scripts do not give is zero.
  void analytical(const physics::AnalyticalGhostState& state,
                  recording::ConditionalIndicesTable& indicesTable,
                  LTS::Layer& layer,
                  std::size_t face,
                  double time,
                  double timeStepWidth,
                  void* stream) {
    using namespace seissol::recording;
    const ConditionalKey key(*KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    auto* cells = indicesTable[key].get(inner_keys::Indices::Id::Cells);
    const std::size_t faces = cells->getSize();
    auto* ghost = reinterpret_cast<real*>(
        layer.var<LTS::AnalyticScratch>(Cfg(), initializer::AllocationPlace::Device));

    auto& side = analytical_[face];
    if (!side.made) {
      side.made = true;
      side.nodes = nodesOf(layer, cells->getHostData(), face);
      typename GlobalBatch::Data data = dataOf(layer, *cells, side.nodes);
      data.ghost = ghost;
      data.ghostStride = tensor::INodal<Cfg>::size();
      // one kernel per program: fused simulations of one script share theirs
      std::vector<const physics::GhostStateIntegral*> integrals;
      for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
        const auto* integral = &state.integral(s);
        const auto found = std::find(integrals.begin(), integrals.end(), integral);
        side.batchOf.push_back(static_cast<std::size_t>(found - integrals.begin()));
        if (found == integrals.end()) {
          integrals.push_back(integral);
          side.batches.push_back(
              std::make_unique<GlobalBatch>(*integral, data, optionsFor(stream, quiet_)));
          checkKind(side.batches.back()->kind(), "the analytical boundary");
          quiet_ = true;
        }
      }
    }

    // what the scripts do not give is zero
    device::DeviceInstance::instance().algorithms().fillArray(
        ghost, static_cast<real>(0), faces * tensor::INodal<Cfg>::size(), stream);
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      side.batches[side.batchOf[s]]->evaluate(
          time, timeStepWidth, s, static_cast<double>(state.scriptSimulation(s)), stream);
    }
  }

  /// The ghost state of the nonlinear Dirichlet boundary at the faces `face` of the cells of
  /// `layer`, integrated over the step [time, time + timeStepWidth]: the integral of the state
  /// inside (projected from the integrated state), and where the script gives a quantity, the
  /// integral of the script (condition.deviceIntegral()), from the state inside at the time points
  /// of the step if it reads it.
  void nonlinearDirichlet(const physics::NonlinearDirichlet& condition,
                          recording::ConditionalPointersToRealsTable& dataTable,
                          recording::ConditionalIndicesTable& indicesTable,
                          LTS::Layer& layer,
                          std::size_t face,
                          double time,
                          double timeStepWidth,
                          parallel::runtime::StreamRuntime& runtime) {
#ifdef SEISSOL_KERNELS_LINEARCK
    if constexpr (Cfg::Solver == SolverType::LinearCK) {
      using namespace seissol::recording;
      using Scratch = kernels::NonlinearDirichletScratch<Cfg>;
      const ConditionalKey key(
          *KernelNames::BoundaryConditions, *ComputationKind::NonlinearDirichlet, face);
      auto& entry = dataTable[key];
      auto* cells = indicesTable[key].get(inner_keys::Indices::Id::Cells);
      const std::size_t faces = cells->getSize();
      const auto* integral = condition.deviceIntegral();
      if (integral == nullptr) {
        logError() << "The nonlinear Dirichlet boundary has no program for the device.";
        return;
      }
      const bool aligned = condition.faceAligned();
      void* stream = runtime.stream();
      auto* ghost = reinterpret_cast<real*>(
          layer.var<LTS::AnalyticScratch>(Cfg(), initializer::AllocationPlace::Device));
      auto* scratch = reinterpret_cast<real*>(
          layer.var<LTS::NonlinearDirichletScratch>(Cfg(), initializer::AllocationPlace::Device));

      auto& side = nonlinear_[face];
      if (!side.made) {
        side.made = true;
        side.nodes = nodesOf(layer, cells->getHostData(), face);
        auto data = dataOf(layer, *cells, side.nodes);
        data.inner = scratch + Scratch::Inner;
        data.innerStride = Scratch::PerFace;
        if (aligned) {
          typename AlignedBatch::Data alignedData;
          alignedData.faces = data.faces;
          alignedData.nodes = data.nodes;
          alignedData.parameters = data.parameters;
          alignedData.cells = data.cells;
          alignedData.inner = data.inner;
          alignedData.innerStride = data.innerStride;
          alignedData.ghost = scratch + Scratch::Aligned;
          alignedData.ghostStride = Scratch::PerFace;
          side.aligned =
              std::make_unique<AlignedBatch>(*integral, alignedData, optionsFor(stream, quiet_));
          checkKind(side.aligned->kind(), "the nonlinear Dirichlet boundary");
        } else {
          data.ghost = ghost;
          data.ghostStride = tensor::INodal<Cfg>::size();
          side.batches.push_back(
              std::make_unique<GlobalBatch>(*integral, data, optionsFor(stream, quiet_)));
          checkKind(side.batches.back()->kind(), "the nonlinear Dirichlet boundary");
        }
        quiet_ = true;
      }

      auto** integrated = entry.get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr();
      auto** blocks = entry.get(inner_keys::Wp::Id::NonlinearDirichletScratch)->getDeviceDataPtr();
      auto** ghosts = entry.get(inner_keys::Wp::Id::Analytical)->getDeviceDataPtr();
      auto** derivatives = entry.get(inner_keys::Wp::Id::Derivatives)->getDeviceDataPtr();
      auto** rotations = entry.get(inner_keys::Wp::Id::T)->getDeviceDataPtr();
      auto** inverseRotations = entry.get(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();

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

      // the state inside at the time points of the step, where the script reads it
      if (integral->innerCount() > 0) {
        std::array<double, Cfg::ConvergenceOrder> times{};
        integral->timePoints(time, timeStepWidth, times.data());
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
          stateAt(project, times, time, derivatives, blocks);
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
          stateAt(project, times, time, derivatives, blocks);
          launch(project, faces, runtime, face);
        }
      }

      // the script, where it gives a quantity
      for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
        if (aligned) {
          side.aligned->evaluate(time, timeStepWidth, s, static_cast<double>(s), stream);
        } else {
          side.batches.front()->evaluate(time, timeStepWidth, s, static_cast<double>(s), stream);
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
    static_cast<void>(time);
    static_cast<void>(timeStepWidth);
    static_cast<void>(runtime);
    logError() << "The nonlinear Dirichlet boundary needs the Taylor series of the state.";
#endif
  }

  private:
  using GlobalBatch = kernels::GhostStateBatch<Cfg, false>;
  using AlignedBatch = kernels::GhostStateBatch<Cfg, true>;

  static_assert(std::is_same_v<unsigned, std::uint32_t>,
                "the cells of a batch are the cell index of a column of a program");

  /// The kernels of the faces of one side of the cells of the layer.
  struct Side {
    bool made{false};
    /// the nodes of the faces, on the device (GhostStateBatch::Data::nodes)
    memory::MemkindArray<real> nodes{memory::Memkind::DeviceGlobalMemory};
    /// in global coordinates: one per program (the analytical boundary: of the fields)
    std::vector<std::unique_ptr<GlobalBatch>> batches;
    /// per simulation, its program
    std::vector<std::size_t> batchOf;
    /// in the face-aligned basis (the nonlinear Dirichlet boundary)
    std::unique_ptr<AlignedBatch> aligned;
  };

  /// The nodes of the faces `face` of the cells `cells` of `layer`, on the device.
  static memory::MemkindArray<real>
      nodesOf(LTS::Layer& layer, const std::vector<unsigned>& cells, std::size_t face) {
    const auto* mappings = layer.var<LTS::BoundaryMapping>(Cfg());
    std::vector<real> nodes(cells.size() * NodeCount * 3);
    for (std::size_t e = 0; e < cells.size(); ++e) {
      const auto* source = mappings[cells[e]][face].nodes;
      std::copy_n(source, NodeCount * 3, nodes.data() + e * NodeCount * 3);
    }
    return memory::MemkindArray<real>(nodes, memory::Memkind::DeviceGlobalMemory);
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

  /// The Taylor series of the state of the cells of the batch at the times `times` of the step
  /// that starts at `start`, as the host evaluates it (kernel::derivativeTaylorExpansion).
  template <typename Kernel>
  static void stateAt(Kernel& project,
                      const std::array<double, Cfg::ConvergenceOrder>& times,
                      double start,
                      real** derivatives,
                      real** blocks) {
    std::size_t offset = 0;
    for (std::size_t d = 0; d < Cfg::ConvergenceOrder; ++d) {
      project.dQ(d) = const_cast<const real**>(derivatives);
      project.extraOffset_dQ(d) = static_cast<int>(offset);
      offset += tensor::dQ<Cfg>::size(d);
    }
    for (std::size_t c = 0; c < Cfg::ConvergenceOrder; ++c) {
      const double tau = times[c] - start;
      double power = 1;
      for (std::size_t d = 0; d < Cfg::ConvergenceOrder; ++d) {
        project.taylorCoefficient(c * Cfg::ConvergenceOrder + d) = static_cast<real>(power);
        power *= tau / static_cast<double>(d + 1);
      }
    }
    project.faceStateI = blocks;
    project.extraOffset_faceStateI =
        static_cast<int>(kernels::NonlinearDirichletScratch<Cfg>::State);
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
  std::array<Side, Cell::NumFaces> analytical_;
  std::array<Side, Cell::NumFaces> nonlinear_;
  /// whether a kernel was made already, which said how it went
  bool quiet_{false};
};

} // namespace seissol::kernels::solver::linearck

#endif // ACL_DEVICE

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DEVICEGHOSTSTATES_H_
