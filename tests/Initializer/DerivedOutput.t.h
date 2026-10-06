// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/SderivFrontend.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "IO/Instance/Geometry/Points.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "Initializer/InitProcedure/DerivedOutput.h"
#include "Memory/MemoryAllocator.h"
#include "Numerical/Projection.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/LuaTracer.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace seissol::unit_test::derived_output {

using namespace seissol::initializer;
namespace projection = seissol::numerical::projection;
using reader::scripting::DataTable;
using reader::scripting::Direction;

using Cfg = Config;
using RealT = Real<Cfg>;
using MaterialT = model::MaterialOf<Cfg>;

// the stride between two quantities of the coefficients of a cell
constexpr std::size_t QuantityStride =
    tensor::Q<Cfg>::Size / tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
constexpr std::size_t CellStride = tensor::Q<Cfg>::size();
constexpr std::size_t Simulations = Cfg::NumSimulations;

/// Cells with random coefficients, their time integrals, and a random straight-sided shape each.
struct Cells {
  std::size_t count{0};
  std::vector<RealT> dofs;
  std::vector<RealT> integrals;
  std::vector<double> jacobians; // 9 per cell, row-major d xi_k / d x_d
  std::vector<seissol::geometry::AffineTransform> shapes;
  // 12 per cell: the origin, then the images of the reference unit vectors less the origin
  std::vector<double> transforms;

  Cells(std::size_t count, unsigned seed) : count(count) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> value(-1.0, 1.0);
    dofs.resize(count * CellStride);
    integrals.resize(count * CellStride);
    for (auto& entry : dofs) {
      entry = static_cast<RealT>(value(rng));
    }
    for (auto& entry : integrals) {
      entry = static_cast<RealT>(value(rng));
    }
    const auto barycenter =
        seissol::geometry::CellTransform::VectorEigenT(Cell::ReferenceBarycenter.data());
    for (std::size_t cell = 0; cell < count; ++cell) {
      std::array<std::array<double, 3>, 4> vertices{};
      const std::array<double, 3> origin{value(rng), value(rng), value(rng)};
      vertices[0] = origin;
      for (std::size_t k = 1; k < 4; ++k) {
        for (std::size_t d = 0; d < 3; ++d) {
          vertices[k][d] = origin[d] + (k - 1 == d ? 1.0 : 0.0) + 0.3 * value(rng);
        }
      }
      const seissol::geometry::AffineTransform transform(vertices);
      const auto inverse = transform.refToSpaceJacobianInverse(barycenter);
      for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t d = 0; d < 3; ++d) {
          jacobians.push_back(inverse(k, d));
        }
      }
      shapes.push_back(transform);
      transforms.insert(transforms.end(), origin.begin(), origin.end());
      for (std::size_t k = 1; k < 4; ++k) {
        for (std::size_t d = 0; d < 3; ++d) {
          transforms.push_back(vertices[k][d] - origin[d]);
        }
      }
    }
  }

  void perturb(unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> value(-1.0, 1.0);
    for (auto& entry : dofs) {
      entry = static_cast<RealT>(value(rng));
    }
  }
};

/// The quantities and their time integrals, as a configuration offers them.
std::vector<DerivedSource> sources() {
  std::vector<DerivedSource> result;
  for (const auto& name : MaterialT::Quantities) {
    result.push_back(DerivedSource{name, false});
  }
  for (const auto& name : MaterialT::Quantities) {
    result.push_back(DerivedSource{"int_" + name, false});
  }
  return result;
}

DerivedGeometry refinedGeometry(std::size_t degree) {
  DerivedGeometry geometry;
  geometry.subcells = io::instance::geometry::subdivideMaps(
      io::instance::geometry::unrefined<3>(), io::instance::geometry::TetrahedronRefine4);
  geometry.dataBase = io::instance::geometry::pointsTetrahedron(degree);
  geometry.dataOrder = degree;
  geometry.order = Cfg::ConvergenceOrder;
  return geometry;
}

/// A program bound to cells, as the volume output binds it: one block per quantity, the inverse
/// Jacobian per cell, and one column per output.
struct Evaluation {
  const DerivedProgram* derived;
  std::size_t numPoints;
  std::vector<double> values;
  DataTable table;
  expr::Binding binding;

  Evaluation(const DerivedProgram& program, Cells& cells, std::size_t simulation = 0)
      : derived(&program), numPoints(cells.count * program.pointsPerCell()),
        values(program.program().outputs().size() * numPoints,
               std::numeric_limits<double>::quiet_NaN()),
        table(numPoints), binding(bindAll(program, cells, simulation)) {}

  expr::Binding bindAll(const DerivedProgram& program, Cells& cells, std::size_t simulation) {
    const auto all = sources();
    program.bindMatrices(table);
    for (std::size_t b = 0; b < program.usedSources().size(); ++b) {
      const std::size_t source = program.usedSources()[b];
      const std::size_t quantity = source % MaterialT::Quantities.size();
      const auto& storage = source < MaterialT::Quantities.size() ? cells.dofs : cells.integrals;
      table.bindBlock<RealT>(program.program().blocks()[b].name,
                             storage.data() + quantity * QuantityStride + simulation,
                             program.program().blocks()[b].length,
                             CellStride,
                             Simulations);
    }
    program.bindGeometry(table, cells.transforms.data());
    if (program.readsJacobian()) {
      for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t d = 0; d < 3; ++d) {
          table.bindCellView<double>("jinv" + std::to_string(k) + std::to_string(d),
                                     cells.jacobians.data(),
                                     program.pointsPerCell(),
                                     9,
                                     k * 3 + d);
        }
      }
    }
    for (std::size_t j = 0; j < program.program().outputs().size(); ++j) {
      table.bindView<double>(
          program.program().outputs()[j].name, Direction::Out, values.data() + j * numPoints);
    }
    return expr::Binding::bind(program.program(), table);
  }

  [[nodiscard]] double value(const std::string& output, std::size_t point) const {
    const auto& outputs = derived->program().outputs();
    for (std::size_t j = 0; j < outputs.size(); ++j) {
      if (outputs[j].name == output) {
        return values[j * numPoints + point];
      }
    }
    FAIL("no output " << output);
    return 0;
  }
};

std::unique_ptr<expr::Kernel> kernelFor(Evaluation& evaluation,
                                        expr::BackendKind backend,
                                        reader::datafield::GridStore& store) {
  expr::BackendOptions options;
  options.preferred = backend;
  auto kernel = expr::makeKernel(evaluation.derived->program(), evaluation.binding, store, options);
  kernel->precompute(evaluation.table);
  return kernel;
}

/// The wave field outputs as master computed them, one projection per call through the generated
/// kernel, and the chain rule per output. Kept verbatim apart from the plumbing, so that the
/// comparison is against the code the program replaces.
class HandWritten {
  public:
  HandWritten(const DerivedGeometry& geometry, std::size_t degree)
      : degree_(degree), pointsPerSubcell_(geometry.dataBase.size()) {
    const auto makeTable = [&](std::optional<std::size_t> derivative) {
      projection::Spec spec;
      spec.target = geometry.target;
      spec.derivative = derivative;
      const auto index = tensor::collvv<Cfg>::index(Cfg::ConvergenceOrder, degree);
      const std::size_t stride =
          tensor::collvv<Cfg>::Size[index] / tensor::collvv<Cfg>::Shape[index][1];
      return std::make_shared<projection::Table<3, 3, RealT>>(geometry.subcells,
                                                              geometry.dataBase,
                                                              geometry.dataOrder,
                                                              stride,
                                                              spec,
                                                              1,
                                                              Cfg::ConvergenceOrder);
    };
    proj_ = makeTable(std::nullopt);
    for (std::size_t direction = 0; direction < 3; ++direction) {
      projD_[direction] = makeTable(direction);
    }
  }

  void projectVolume(double* target, const RealT* dofsSingleQuantity, const RealT* collvv) const {
    constexpr auto Variant = configIdOf<Cfg>();
    runtime::kernel::projectBasisToVtkVolume vtkproj{};
    memory::AlignedArray<RealT, Cfg::NumSimulations> simselect{};
    alignas(Alignment) std::array<RealT, MaxVtk3dPoints> alignedTarget{};
    simselect[0] = 1;
    vtkproj.simselect = runtime::init::simselect::view(Variant, simselect.data());
    vtkproj.qb = runtime::init::qb::view(Variant, dofsSingleQuantity);
    vtkproj.xv(degree_) = runtime::init::xv::view(Variant, degree_, alignedTarget.data());
    vtkproj.collvv(Cfg::ConvergenceOrder, degree_) =
        runtime::init::collvv::view(Variant, Cfg::ConvergenceOrder, degree_, collvv);
    vtkproj.execute(Variant, degree_);
    std::copy_n(alignedTarget.data(), pointsPerSubcell_, target);
  }

  void projectVolumeDeriv(double* target,
                          const RealT* dofsSingleQuantity,
                          std::size_t dir,
                          const double* grad,
                          std::size_t subcell) const {
    std::array<double, MaxVtk3dPoints> dataX{};
    std::array<double, MaxVtk3dPoints> dataY{};
    std::array<double, MaxVtk3dPoints> dataZ{};
    projectVolume(dataX.data(), dofsSingleQuantity, (*projD_[0])(subcell, Cfg::ConvergenceOrder));
    projectVolume(dataY.data(), dofsSingleQuantity, (*projD_[1])(subcell, Cfg::ConvergenceOrder));
    projectVolume(dataZ.data(), dofsSingleQuantity, (*projD_[2])(subcell, Cfg::ConvergenceOrder));
    for (std::size_t i = 0; i < pointsPerSubcell_; ++i) {
      target[i] = dataX[i] * grad[0 * 3 + dir] + dataY[i] * grad[1 * 3 + dir] +
                  dataZ[i] * grad[2 * 3 + dir];
    }
  }

  void value(double* target, const RealT* dofsSingleQuantity, std::size_t subcell) const {
    projectVolume(target, dofsSingleQuantity, (*proj_)(subcell, Cfg::ConvergenceOrder));
  }

  void strain(double* target,
              const RealT* integrals,
              std::size_t idx1,
              std::size_t idx2,
              const double* grad,
              std::size_t subcell) const {
    const auto* dofsSingleQuantity1 =
        integrals + QuantityStride * (idx1 + MaterialT::VelocityOffset);
    projectVolumeDeriv(target, dofsSingleQuantity1, idx2, grad, subcell);
    if (idx1 != idx2) {
      const auto* dofsSingleQuantity2 =
          integrals + QuantityStride * (idx2 + MaterialT::VelocityOffset);
      std::array<double, MaxVtk3dPoints> itarget{};
      projectVolumeDeriv(itarget.data(), dofsSingleQuantity2, idx1, grad, subcell);
      for (std::size_t i = 0; i < pointsPerSubcell_; ++i) {
        target[i] = (target[i] + itarget[i]) / 2;
      }
    }
  }

  void rotation(double* target,
                const RealT* dofs,
                std::size_t idx1,
                std::size_t idx2,
                const double* grad,
                std::size_t subcell) const {
    const auto* dofsSingleQuantity1 = dofs + QuantityStride * (idx1 + MaterialT::VelocityOffset);
    projectVolumeDeriv(target, dofsSingleQuantity1, idx2, grad, subcell);
    const auto* dofsSingleQuantity2 = dofs + QuantityStride * (idx2 + MaterialT::VelocityOffset);
    std::array<double, MaxVtk3dPoints> itarget{};
    projectVolumeDeriv(itarget.data(), dofsSingleQuantity2, idx1, grad, subcell);
    for (std::size_t i = 0; i < pointsPerSubcell_; ++i) {
      target[i] -= itarget[i];
    }
  }

  static constexpr std::size_t MaxVtk3dPoints = tensor::vtk3d<Cfg>::Shape
      [(sizeof(tensor::vtk3d<Cfg>::Shape) / sizeof(tensor::vtk3d<Cfg>::Shape[0])) - 1][1];

  private:
  std::size_t degree_;
  std::size_t pointsPerSubcell_;
  std::shared_ptr<projection::Table<3, 3, RealT>> proj_;
  std::array<std::shared_ptr<projection::Table<3, 3, RealT>>, 3> projD_;
};

WaveFieldSelection fullSelection() {
  WaveFieldSelection selection;
  selection.quantities.assign(MaterialT::Quantities.begin(), MaterialT::Quantities.end());
  selection.velocityOffset = MaterialT::VelocityOffset;
  selection.outputMask.assign(selection.quantities.size(), true);
  selection.integrationMask.assign(selection.quantities.size(), true);
  selection.strain = true;
  selection.rotation = true;
  return selection;
}

using Index = std::tuple<std::string, std::size_t, std::size_t>;
const std::array<Index, 6> StrainIndices{Index{"xx", 0, 0},
                                         Index{"yy", 1, 1},
                                         Index{"zz", 2, 2},
                                         Index{"xy", 0, 1},
                                         Index{"yz", 1, 2},
                                         Index{"xz", 0, 2}};
const std::array<Index, 3> RotationIndices{Index{"1", 2, 1}, Index{"2", 0, 2}, Index{"3", 1, 0}};

// What the generated kernel computes in the precision of the build, and possibly with fused
// multiply-adds, so the comparison cannot be bitwise.
const double Tolerance = 1000 * std::numeric_limits<RealT>::epsilon();

TEST_CASE("DerivedOutput: the built-in program reproduces the hand-written wave field outputs") {
  constexpr std::size_t Degree = 2;
  const auto geometry = refinedGeometry(Degree);
  Cells cells(7, 1);

  const DerivedProgram derived(waveFieldProgram(fullSelection()), sources(), geometry);
  REQUIRE(derived.pointsPerCell() == geometry.subcells.size() * geometry.dataBase.size());
  Evaluation evaluation(derived, cells);
  reader::datafield::GridStore store;
  kernelFor(evaluation, expr::BackendKind::Interpreter, store)->run(evaluation.table);

  const HandWritten master(geometry, Degree);
  const std::size_t pointsPerSubcell = derived.pointsPerSubcell();
  std::vector<double> reference(pointsPerSubcell);
  const auto check = [&](const std::string& output, std::size_t cell, std::size_t subcell) {
    for (std::size_t i = 0; i < pointsPerSubcell; ++i) {
      const std::size_t point = cell * derived.pointsPerCell() + subcell * pointsPerSubcell + i;
      REQUIRE(evaluation.value(output, point) == doctest::Approx(reference[i]).epsilon(Tolerance));
    }
  };

  for (std::size_t cell = 0; cell < cells.count; ++cell) {
    const RealT* dofs = cells.dofs.data() + cell * CellStride;
    const RealT* integrals = cells.integrals.data() + cell * CellStride;
    const double* grad = cells.jacobians.data() + cell * 9;
    for (std::size_t subcell = 0; subcell < geometry.subcells.size(); ++subcell) {
      for (std::size_t q = 0; q < MaterialT::Quantities.size(); ++q) {
        master.value(reference.data(), dofs + q * QuantityStride, subcell);
        check(MaterialT::Quantities[q], cell, subcell);
        master.value(reference.data(), integrals + q * QuantityStride, subcell);
        check("int-" + MaterialT::Quantities[q], cell, subcell);
      }
      for (const auto& [name, idx1, idx2] : StrainIndices) {
        master.strain(reference.data(), integrals, idx1, idx2, grad, subcell);
        check("eps" + name, cell, subcell);
      }
      for (const auto& [name, idx1, idx2] : RotationIndices) {
        master.rotation(reference.data(), dofs, idx1, idx2, grad, subcell);
        check("rot" + name, cell, subcell);
      }
    }
  }
}

namespace {
// Products and sums through volatile, so that the reference rounds every operation on its own,
// like the interpreter, whatever this translation unit is compiled with.
double multiply(double a, double b) {
  const volatile double product = a * b;
  return product;
}
double add(double a, double b) {
  const volatile double sum = a + b;
  return sum;
}
} // namespace

TEST_CASE("DerivedOutput: contractions, chain rule and stacking are exact") {
  constexpr std::size_t Degree = 1;
  const auto geometry = refinedGeometry(Degree);
  Cells cells(5, 2);

  const DerivedProgram derived(waveFieldProgram(fullSelection()), sources(), geometry);
  Evaluation interpreted(derived, cells);
  Evaluation compiled(derived, cells);
  reader::datafield::GridStore store;
  kernelFor(interpreted, expr::BackendKind::Interpreter, store)->run(interpreted.table);
  auto compiledKernel = kernelFor(compiled, expr::BackendKind::RtcCpu, store);
  compiledKernel->run(compiled.table);

  // the compiled kernel against the interpreter, where there is a compiler
  if (compiledKernel->kind() == expr::BackendKind::RtcCpu) {
    CHECK(std::memcmp(interpreted.values.data(),
                      compiled.values.data(),
                      interpreted.values.size() * sizeof(double)) == 0);
  }

  // and the interpreter against the formula, in double, summed in ascending order
  const std::size_t modes = projection::modalSize(3, Cfg::ConvergenceOrder);
  const std::size_t pointsPerSubcell = derived.pointsPerSubcell();
  std::array<numerical::DenseMatrix, 3> matrices;
  for (std::size_t subcell = 0; subcell < geometry.subcells.size(); ++subcell) {
    for (std::size_t k = 0; k < 3; ++k) {
      projection::Spec spec;
      spec.order = Cfg::ConvergenceOrder;
      spec.derivative = k;
      matrices[k] = projection::build<3, 3>(
          geometry.dataBase, geometry.dataOrder, geometry.subcells[subcell], spec);
    }
    for (std::size_t cell = 0; cell < cells.count; ++cell) {
      const double* grad = cells.jacobians.data() + cell * 9;
      const auto derivative = [&](const RealT* coefficients, std::size_t point, std::size_t axis) {
        double sum = 0;
        for (std::size_t k = 0; k < 3; ++k) {
          double reference = 0;
          for (std::size_t m = 0; m < modes; ++m) {
            reference = add(reference,
                            multiply(matrices[k](point, m),
                                     static_cast<double>(coefficients[m * Simulations])));
          }
          const double term = multiply(reference, grad[k * 3 + axis]);
          sum = k == 0 ? term : add(sum, term);
        }
        return sum;
      };
      const RealT* dofs = cells.dofs.data() + cell * CellStride;
      const RealT* integrals = cells.integrals.data() + cell * CellStride;
      const auto velocity = [&](const RealT* base, std::size_t i) {
        return base + (MaterialT::VelocityOffset + i) * QuantityStride;
      };
      for (std::size_t p = 0; p < pointsPerSubcell; ++p) {
        const std::size_t point = cell * derived.pointsPerCell() + subcell * pointsPerSubcell + p;
        for (const auto& [name, i, j] : StrainIndices) {
          const double first = derivative(velocity(integrals, i), p, j);
          const double expected =
              i == j ? first : add(first, derivative(velocity(integrals, j), p, i)) / 2.0;
          REQUIRE(interpreted.value("eps" + name, point) == expected);
        }
        for (const auto& [name, i, j] : RotationIndices) {
          const double expected =
              derivative(velocity(dofs, i), p, j) - derivative(velocity(dofs, j), p, i);
          REQUIRE(interpreted.value("rot" + name, point) == expected);
        }
      }
    }
  }
}

TEST_CASE("DerivedOutput: a program written in sderiv gives the built-in outputs bitwise") {
  constexpr std::size_t Degree = 1;
  const auto geometry = refinedGeometry(Degree);
  Cells cells(3, 3);

  auto selection = fullSelection();
  selection.outputMask.assign(selection.quantities.size(), false);
  selection.integrationMask.assign(selection.quantities.size(), false);
  const DerivedProgram builtIn(waveFieldProgram(selection), sources(), geometry);

  // the strain and rotation as a model author would write them, with the chain rule spelled out
  std::string source;
  const auto v = [](std::size_t i) { return MaterialT::Quantities[MaterialT::VelocityOffset + i]; };
  const std::string axes = "xyz";
  for (const auto& prefix : {std::string("int_"), std::string()}) {
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t d = 0; d < 3; ++d) {
        source += "def d" + std::string(1, axes[d]) + "_" + prefix + v(i) + " = ";
        for (std::size_t k = 0; k < 3; ++k) {
          source += (k == 0 ? "" : " + ") + prefix + v(i) + "_r" + std::to_string(k) + " * jinv" +
                    std::to_string(k) + std::to_string(d);
        }
        source += "\n";
      }
    }
  }
  for (const auto& [name, i, j] : StrainIndices) {
    const std::string first = std::string("d") + axes[j] + "_int_" + v(i);
    source +=
        "out def eps" + name + " = " +
        (i == j ? first
                : "(" + first + " + d" + std::string(1, axes[i]) + "_int_" + v(j) + ") / 2.0") +
        "\n";
  }
  for (const auto& [name, i, j] : RotationIndices) {
    source += "out def rot" + name + " = d" + std::string(1, axes[j]) + "_" + v(i) + " - d" +
              std::string(1, axes[i]) + "_" + v(j) + "\n";
  }
  const DerivedProgram written(expr::compileSderivModule(source), sources(), geometry);

  Evaluation first(builtIn, cells);
  Evaluation second(written, cells);
  reader::datafield::GridStore store;
  kernelFor(first, expr::BackendKind::Interpreter, store)->run(first.table);
  kernelFor(second, expr::BackendKind::Interpreter, store)->run(second.table);
  REQUIRE(first.values.size() == second.values.size());
  CHECK(std::memcmp(
            first.values.data(), second.values.data(), first.values.size() * sizeof(double)) == 0);
}

TEST_CASE("DerivedOutput: a maximum over time follows a hand-kept reference") {
  constexpr std::size_t Degree = 1;
  const auto geometry = refinedGeometry(Degree);
  Cells cells(4, 4);

  const std::string v1 = MaterialT::Quantities[MaterialT::VelocityOffset];
  const std::string v2 = MaterialT::Quantities[MaterialT::VelocityOffset + 1];
  const std::string v3 = MaterialT::Quantities[MaterialT::VelocityOffset + 2];
  const std::string speed =
      "sqrt(" + v1 + "*" + v1 + " + " + v2 + "*" + v2 + " + " + v3 + "*" + v3 + ")";
  const DerivedProgram derived(expr::compileSderivModule("state pgv = 0.0\n"
                                                         "out def speed = " +
                                                         speed +
                                                         "\n"
                                                         "out def pgv = max(pgv, speed)\n"),
                               sources(),
                               geometry);
  REQUIRE(derived.program().state().size() == 1);

  Evaluation evaluation(derived, cells);
  reader::datafield::GridStore store;
  auto kernel = kernelFor(evaluation, expr::BackendKind::Interpreter, store);

  std::vector<double> reference(evaluation.numPoints, 0.0);
  for (unsigned call = 0; call < 6; ++call) {
    cells.perturb(100 + call);
    kernel->run(evaluation.table);
    for (std::size_t point = 0; point < evaluation.numPoints; ++point) {
      reference[point] = std::max(reference[point], evaluation.value("speed", point));
      REQUIRE(evaluation.value("pgv", point) == reference[point]);
    }
  }
}

TEST_CASE("DerivedOutput: a Lua program computes what its sderiv counterpart does") {
  constexpr std::size_t Degree = 1;
  const auto geometry = refinedGeometry(Degree);
  Cells cells(3, 6);

  const std::string v1 = MaterialT::Quantities[MaterialT::VelocityOffset];
  const std::string v2 = MaterialT::Quantities[MaterialT::VelocityOffset + 1];
  const std::string v3 = MaterialT::Quantities[MaterialT::VelocityOffset + 2];
  const std::string lua = "local M = {}\n"
                          "M.state = { pgv = 0.0 }\n"
                          "function M.evaluate(fields, " +
                          v1 + ", " + v2 + ", " + v3 + ", dx_" + v1 + ", dy_" + v2 + ", dz_" + v3 +
                          ", pgv)\n"
                          "  return { pgv = math.max(pgv, math.sqrt(" +
                          v1 + "*" + v1 + " + " + v2 + "*" + v2 + " + " + v3 + "*" + v3 +
                          ")),\n"
                          "           divv = dx_" +
                          v1 + " + dy_" + v2 + " + dz_" + v3 +
                          " }\n"
                          "end\n"
                          "return M\n";
  const std::string sderiv = "state pgv = 0.0\n"
                             "out def pgv = max(pgv, sqrt(" +
                             v1 + "*" + v1 + " + " + v2 + "*" + v2 + " + " + v3 + "*" + v3 +
                             "))\n"
                             "out def divv = dx_" +
                             v1 + " + dy_" + v2 + " + dz_" + v3 + "\n";

  reader::scripting::TraceFailure failure;
  auto traced = reader::scripting::traceLuaModule(lua, {}, failure);
  REQUIRE_MESSAGE(traced.has_value(), failure.reason);
  const DerivedProgram fromLua(std::move(*traced), sources(), geometry);
  const DerivedProgram fromSderiv(expr::compileSderivModule(sderiv), sources(), geometry);

  Evaluation first(fromLua, cells);
  Evaluation second(fromSderiv, cells);
  reader::datafield::GridStore store;
  auto firstKernel = kernelFor(first, expr::BackendKind::Interpreter, store);
  auto secondKernel = kernelFor(second, expr::BackendKind::Interpreter, store);
  for (unsigned call = 0; call < 3; ++call) {
    cells.perturb(200 + call);
    firstKernel->run(first.table);
    secondKernel->run(second.table);
    for (const auto* output : {"pgv", "divv"}) {
      for (std::size_t point = 0; point < first.numPoints; ++point) {
        const double a = first.value(output, point);
        const double b = second.value(output, point);
        REQUIRE(std::memcmp(&a, &b, sizeof(double)) == 0);
      }
    }
  }
}

TEST_CASE("DerivedOutput: shared reference derivatives are contracted once") {
  const auto geometry = refinedGeometry(1);
  auto selection = fullSelection();
  selection.outputMask.assign(selection.quantities.size(), false);
  selection.integrationMask.assign(selection.quantities.size(), false);

  const auto contractions = [&](bool strain, bool rotation) {
    selection.strain = strain;
    selection.rotation = rotation;
    const DerivedProgram derived(waveFieldProgram(selection), sources(), geometry);
    const auto lowered = expr::lower(derived.program());
    return std::count_if(lowered.run().code.begin(),
                         lowered.run().code.end(),
                         [](const expr::Instruction& instruction) {
                           return instruction.op == expr::Opcode::Contract;
                         });
  };
  // the hand-written outputs project each velocity along all three reference directions for
  // every derivative they take: 27 projections for the strain, 18 for the rotation
  CHECK(contractions(true, false) == 9);
  CHECK(contractions(false, true) == 9);
  // the strain reads the time integral and the rotation the solution, so they share nothing
  CHECK(contractions(true, true) == 18);
}

TEST_CASE("DerivedOutput: the vocabulary is checked") {
  const auto geometry = refinedGeometry(1);
  const auto unknown = expr::compileSderivModule("out def a = v7 * 2.0\n");
  CHECK_THROWS_AS(DerivedProgram(unknown, sources(), geometry), std::invalid_argument);

  const auto known = expr::compileSderivModule("out def a = dz_" + MaterialT::Quantities[0] +
                                               " + jinv12 + x * y * z + t + dt\n");
  const DerivedProgram derived(known, sources(), geometry);
  CHECK(derived.readsJacobian());
  CHECK(derived.readsCoordinates());
  CHECK(derived.readsTime());
  CHECK(derived.readsTimeStep());
  CHECK(derived.usedSources() == std::vector<std::size_t>{0});
  CHECK(derived.referencePoints().size() == derived.pointsPerCell());

  CHECK(readsTimeIntegral("int_v1"));
  CHECK(readsTimeIntegral("dy_int_v1"));
  CHECK(readsTimeIntegral("int_v1_r2"));
  CHECK_FALSE(readsTimeIntegral("dy_v1"));
  CHECK_FALSE(readsTimeIntegral("v1"));
}

TEST_CASE("DerivedOutput: the coordinates are the affine map of the cell at its points") {
  const auto program =
      expr::compileSderivModule("out def px = x\nout def py = y\nout def pz = z\n");
  const auto geometry = refinedGeometry(2);
  Cells cells(5, 8);
  const DerivedProgram derived(program, sources(), geometry);
  CHECK(derived.readsCoordinates());
  CHECK(derived.program().inputs().empty());
  Evaluation evaluation(derived, cells);
  reader::datafield::GridStore store;
  kernelFor(evaluation, expr::BackendKind::Interpreter, store)->run(evaluation.table);
  const auto& reference = derived.referencePoints();
  for (std::size_t cell = 0; cell < cells.count; ++cell) {
    for (std::size_t p = 0; p < reference.size(); ++p) {
      const auto expected = cells.shapes[cell].refToSpace(reference[p]);
      const std::size_t point = cell * reference.size() + p;
      REQUIRE(evaluation.value("px", point) == doctest::Approx(expected[0]).epsilon(1e-14));
      REQUIRE(evaluation.value("py", point) == doctest::Approx(expected[1]).epsilon(1e-14));
      REQUIRE(evaluation.value("pz", point) == doctest::Approx(expected[2]).epsilon(1e-14));
    }
  }
}

TEST_CASE("DerivedOutput: the strain and rotation cost less than by hand" * doctest::skip(true)) {
  // A measurement, not a check: run with --no-skip. Both evaluate the full strain and rotation
  // set at every output point of every cell, as the writer pulls it.
  constexpr std::size_t Degree = 1;
  constexpr std::size_t Repetitions = 5;
  DerivedGeometry geometry;
  geometry.dataBase = io::instance::geometry::pointsTetrahedron(Degree);
  geometry.dataOrder = Degree;
  geometry.order = Cfg::ConvergenceOrder;
  Cells cells(20000, 5);

  auto selection = fullSelection();
  selection.outputMask.assign(selection.quantities.size(), false);
  selection.integrationMask.assign(selection.quantities.size(), false);
  const DerivedProgram derived(waveFieldProgram(selection), sources(), geometry);

  const HandWritten master(geometry, Degree);
  const std::size_t pointsPerSubcell = derived.pointsPerSubcell();
  std::vector<double> target(9 * cells.count * pointsPerSubcell);
  const auto handWritten = [&]() {
    for (std::size_t cell = 0; cell < cells.count; ++cell) {
      const RealT* dofs = cells.dofs.data() + cell * CellStride;
      const RealT* integrals = cells.integrals.data() + cell * CellStride;
      const double* grad = cells.jacobians.data() + cell * 9;
      std::size_t output = 0;
      for (const auto& [name, idx1, idx2] : StrainIndices) {
        master.strain(target.data() + (output++ * cells.count + cell) * pointsPerSubcell,
                      integrals,
                      idx1,
                      idx2,
                      grad,
                      0);
      }
      for (const auto& [name, idx1, idx2] : RotationIndices) {
        master.rotation(target.data() + (output++ * cells.count + cell) * pointsPerSubcell,
                        dofs,
                        idx1,
                        idx2,
                        grad,
                        0);
      }
    }
  };

  const auto measure = [&](const auto& run) {
    run();
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < Repetitions; ++i) {
      run();
    }
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() /
           Repetitions;
  };

  const double byHand = measure(handWritten);
  reader::datafield::GridStore store;
  for (const auto backend : {expr::BackendKind::Interpreter, expr::BackendKind::RtcCpu}) {
    Evaluation evaluation(derived, cells);
    auto kernel = kernelFor(evaluation, backend, store);
    const double program = measure([&]() { kernel->run(evaluation.table); });
    MESSAGE("strain and rotation of " << cells.count << " cells: by hand " << byHand * 1e3
                                      << " ms, " << expr::name(kernel->kind()) << " "
                                      << program * 1e3 << " ms, ratio " << byHand / program);
  }
}

} // namespace seissol::unit_test::derived_output
