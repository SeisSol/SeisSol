// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_GHOSTSTATEBATCH_H_
#define SEISSOL_SRC_KERNELS_GHOSTSTATEBATCH_H_

#include "Common/Real.h"
#include "Equations/Datastructures.h"
#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Model/ScriptParameters.h"
#include "Physics/GhostState.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <utils/logger.h>
#include <vector>

namespace seissol::kernels {

/// The ghost state of a boundary whose script is integrated over a time step at once
/// (physics::GhostStateIntegral), for all faces of a batch in one kernel -- the faces of one side
/// of the cells of a layer, as the recorder batches them for a device. The kernel reads the nodes
/// of the faces, the parameters of the materials of their cells and, where the script reads it,
/// the state inside at the time points of the step, all where they are; it writes the integral of
/// each component the program gives into the tensor of the face, shaped like INodal (or, with
/// `Aligned`, like INodalAligned, in the face-aligned basis). What the program does not give, it
/// leaves as it is.
///
/// The points of the kernel are the nodes of the faces, face by face. A simulation of a fused run
/// is one call, which moves the bases of the state inside and of the ghost state to its slice.
template <typename Cfg, bool Aligned>
class GhostStateBatch {
  public:
  using real = Real<Cfg>;
  using TargetTensor = std::conditional_t<Aligned, tensor::INodalAligned<Cfg>, tensor::INodal<Cfg>>;
  using TargetInit = std::conditional_t<Aligned, init::INodalAligned<Cfg>, init::INodal<Cfg>>;

  static constexpr std::size_t NodeCount = TargetTensor::Shape[multisim::BasisDim<Cfg>];
  static constexpr std::size_t Columns = TargetTensor::Shape[multisim::BasisDim<Cfg> + 1];
  static constexpr std::size_t TargetSize = TargetTensor::size();
  static constexpr std::size_t ParameterCount = model::ScriptParameterCount<model::MaterialOf<Cfg>>;

  /// Where the kernel reads and writes, in the memory it runs in.
  struct Data {
    std::size_t faces{0};
    /// coordinate d of node p of face e at nodes[(e * NodeCount + p) * 3 + d]
    const real* nodes{nullptr};
    /// the parameters of the material of the cell of face e at parameters[cells[e] *
    /// ParameterCount], in the order of model::scriptParameterNames
    const real* parameters{nullptr};
    const std::uint32_t* cells{nullptr};
    /// the state inside of face e at time point t, a tensor like the target, at
    /// inner + e * innerStride + t * TargetSize; read only where the script reads the state inside
    const real* inner{nullptr};
    std::size_t innerStride{0};
    /// the ghost state of face e, a tensor like the target, at ghost + e * ghostStride
    real* ghost{nullptr};
    std::size_t ghostStride{TargetSize};
  };

  /// The kernel of `integral` for the faces of `data`, on the backend `options` asks for; it falls
  /// back as expr::makeKernel does, which kind() tells.
  GhostStateBatch(const physics::GhostStateIntegral& integral,
                  const Data& data,
                  const expr::BackendOptions& options)
      : integral_(&integral), data_(data), preferred_(options.preferred),
        table_(std::make_unique<reader::scripting::DataTable>(data.faces * NodeCount)) {
    using reader::scripting::Direction;
    using Input = physics::GhostStateIntegral::Input;
    const auto& program = integral.program();
    const std::size_t innerCount = integral.innerCount();
    if (innerCount > 0 && data.innerStride < integral.timePointCount() * TargetSize) {
      logError() << "The state inside of a face at the time points of a step needs"
                 << integral.timePointCount() * TargetSize << "reals per face, it has"
                 << data.innerStride << ".";
    }
    const std::size_t nodeStride = nodeStrideOf();

    shifts_.assign(Cfg::NumSimulations * (program.inputs().size() + program.outputs().size()), 0);
    const auto shiftsOf = [&](std::size_t column, std::size_t component) {
      for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
        shifts_[s * columnCount() + column] =
            static_cast<std::ptrdiff_t>(offsetOf(0, component, s)) -
            static_cast<std::ptrdiff_t>(offsetOf(0, component, 0));
      }
    };

    auto& table = *table_;
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      const auto& name = program.inputs()[i].name;
      const auto index = integral.inputIndex(i);
      switch (integral.input(i)) {
      case Input::X:
        table.bindCellPointViewConst<real>(name, data.nodes, NodeCount, NodeCount * 3, 3, 0);
        break;
      case Input::Y:
        table.bindCellPointViewConst<real>(name, data.nodes, NodeCount, NodeCount * 3, 3, 1);
        break;
      case Input::Z:
        table.bindCellPointViewConst<real>(name, data.nodes, NodeCount, NodeCount * 3, 3, 2);
        break;
      case Input::Start:
      case Input::Width:
      case Input::Simulation:
        // the same at all points, and passed by value per call
        table.bindConstant<double>(name, 0.0);
        break;
      case Input::Parameter:
        table.bindCellView<real>(
            name, data.parameters, NodeCount, ParameterCount, index, data.cells);
        break;
      case Input::Inner: {
        const std::size_t point = index / innerCount;
        const std::size_t component = index % innerCount;
        checkColumn(component);
        table.bindCellPointViewConst<real>(name,
                                           data.inner,
                                           NodeCount,
                                           data.innerStride,
                                           nodeStride,
                                           point * TargetSize + offsetOf(0, component, 0));
        shiftsOf(i, component);
        break;
      }
      }
    }
    for (std::size_t j = 0; j < program.outputs().size(); ++j) {
      const auto component = integral.outputComponent(j);
      checkColumn(component);
      table.bindCellPointView<real>(program.outputs()[j].name,
                                    Direction::Out,
                                    data.ghost,
                                    NodeCount,
                                    data.ghostStride,
                                    nodeStride,
                                    offsetOf(0, component, 0));
      shiftsOf(program.inputs().size() + j, component);
    }

    inputs_.resize(program.inputs().size());
    outputs_.resize(program.outputs().size());
    if (program.outputs().empty() || data.faces == 0) {
      // nothing to compute: the script gives none of the components
      return;
    }
    binding_ = std::make_unique<expr::Binding>(expr::Binding::bind(program, table));
    kernel_ = expr::makeKernel(program, *binding_, reader::datafield::sharedGridStore(), options);
    kernel_->precompute(table);
  }

  /// The integral over the step [start, start + width] for the simulation `simulation` of the
  /// run, whose script reads `scriptSimulation` as `sim`, on `stream`.
  void evaluate(
      double start, double width, std::size_t simulation, double scriptSimulation, void* stream) {
    using Input = physics::GhostStateIntegral::Input;
    if (kernel_ == nullptr) {
      return;
    }
    start_ = start;
    width_ = width;
    simulation_ = scriptSimulation;
    const auto& program = integral_->program();
    const auto* shifts = shifts_.data() + simulation * columnCount();
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      switch (integral_->input(i)) {
      case Input::Start:
        inputs_[i] = &start_;
        break;
      case Input::Width:
        inputs_[i] = &width_;
        break;
      case Input::Simulation:
        inputs_[i] = &simulation_;
        break;
      case Input::Inner:
        inputs_[i] = data_.inner + shifts[i];
        break;
      default:
        // where it is bound
        inputs_[i] = nullptr;
        break;
      }
    }
    for (std::size_t j = 0; j < program.outputs().size(); ++j) {
      outputs_[j] = data_.ghost + shifts[program.inputs().size() + j];
    }
    expr::KernelArgs args;
    args.inputs = inputs_.data();
    args.inputCount = inputs_.size();
    args.outputs = outputs_.data();
    args.outputCount = outputs_.size();
    args.first = 0;
    args.count = data_.faces * NodeCount;
    args.stream = stream;
    kernel_->run(args);
  }

  /// The backend of the kernel; the one asked for if there is nothing to compute.
  [[nodiscard]] expr::BackendKind kind() const {
    return kernel_ == nullptr ? preferred_ : kernel_->kind();
  }

  private:
  [[nodiscard]] std::size_t columnCount() const {
    return integral_->program().inputs().size() + integral_->program().outputs().size();
  }

  /// The element of node `node`, component `component` and simulation `simulation` in a tensor
  /// like the target.
  static std::size_t offsetOf(std::size_t node, std::size_t component, std::size_t simulation) {
    static std::array<real, TargetSize> probe{};
    auto view = TargetInit::view::create(probe.data());
    auto slice = multisim::simtensor<Cfg>(view, simulation);
    return static_cast<std::size_t>(&slice(node, component) - probe.data());
  }

  /// The elements from one node of a component to the next: the same for every node, component
  /// and simulation, as a kernel per point of a face reads the tensor.
  static std::size_t nodeStrideOf() {
    const std::size_t stride = NodeCount > 1 ? offsetOf(1, 0, 0) - offsetOf(0, 0, 0) : 1;
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      for (std::size_t k = 0; k < Columns; ++k) {
        for (std::size_t p = 0; p < NodeCount; ++p) {
          if (offsetOf(p, k, s) != offsetOf(0, k, s) + p * stride) {
            logError() << "The nodes of a face are no column with a fixed stride in its tensor.";
          }
        }
      }
    }
    return stride;
  }

  static void checkColumn(std::size_t component) {
    if (component >= Columns) {
      logError() << "The ghost state of a face has" << Columns << "components; a script gives"
                 << "component" << component << ".";
    }
  }

  const physics::GhostStateIntegral* integral_;
  Data data_;
  expr::BackendKind preferred_;
  std::unique_ptr<reader::scripting::DataTable> table_;
  std::unique_ptr<expr::Binding> binding_;
  std::unique_ptr<expr::Kernel> kernel_;
  /// per simulation and column (the inputs, then the outputs), how far its base moves, in reals
  std::vector<std::ptrdiff_t> shifts_;
  std::vector<const void*> inputs_;
  std::vector<void*> outputs_;
  double start_{0};
  double width_{0};
  double simulation_{0};
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_GHOSTSTATEBATCH_H_
