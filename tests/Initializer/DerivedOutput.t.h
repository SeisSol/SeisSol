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
#include "Geometry/FaceTransform.h"
#include "IO/Instance/Geometry/Points.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "Initializer/InitProcedure/DerivedOutput.h"
#include "Memory/MemoryAllocator.h"
#include "Numerical/Projection.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/LuaTracer.h"
#include "Solver/MultipleSimulations.h"
#include "TestHelper.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
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
inline std::vector<DerivedSource> sources() {
  std::vector<DerivedSource> result;
  result.reserve(2 * MaterialT::Quantities.size());
  for (const auto& name : MaterialT::Quantities) {
    result.push_back(DerivedSource{name, Representation::Modal});
  }
  for (const auto& name : MaterialT::Quantities) {
    result.push_back(DerivedSource{"int_" + name, Representation::Modal});
  }
  return result;
}

inline DerivedGeometry refinedGeometry(std::size_t degree) {
  DerivedGeometry geometry;
  geometry.subcells = io::instance::geometry::subdivideMaps(
      io::instance::geometry::unrefined<3>(), io::instance::geometry::TetrahedronRefine4);
  geometry.dataBase = io::instance::geometry::pointsTetrahedron(degree);
  geometry.dataOrder = degree;
  geometry.order = Cfg::ConvergenceOrder;
  return geometry;
}

// the coefficients of the displacement of a face: from one face to the next, and from one
// component to the next
constexpr std::size_t FaceStride = tensor::faceDisplacement<Cfg>::size();
constexpr std::size_t ComponentStride =
    tensor::faceDisplacement<Cfg>::Size /
    tensor::faceDisplacement<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];

/// Faces of the free surface: face i lies on side i % 4 of cell i, and has a random displacement.
struct Faces {
  Cells cells;
  std::vector<RealT> displacements;
  std::vector<seissol::geometry::AffineFaceTransform> shapes;
  // 12 per face: the origin, then the images of the two reference unit vectors of the face less
  // the origin, and zeros
  std::vector<double> transforms;

  Faces(std::size_t count, unsigned seed) : cells(count, seed) {
    std::mt19937 rng(seed + 1000);
    std::uniform_real_distribution<double> value(-1.0, 1.0);
    displacements.resize(count * FaceStride);
    for (auto& entry : displacements) {
      entry = static_cast<RealT>(value(rng));
    }
    using FaceVectorT = seissol::geometry::FaceTransform::FaceVectorT;
    for (std::size_t face = 0; face < count; ++face) {
      shapes.emplace_back(cells.shapes[face], seissol::geometry::ReferenceFaceMap(side(face)));
      const auto origin = shapes.back().refToSpace(FaceVectorT(0.0, 0.0));
      const auto first = shapes.back().refToSpace(FaceVectorT(1.0, 0.0));
      const auto second = shapes.back().refToSpace(FaceVectorT(0.0, 1.0));
      for (std::size_t d = 0; d < 3; ++d) {
        transforms.push_back(origin(d));
      }
      for (std::size_t d = 0; d < 3; ++d) {
        transforms.push_back(first(d) - origin(d));
      }
      for (std::size_t d = 0; d < 3; ++d) {
        transforms.push_back(second(d) - origin(d));
      }
      for (std::size_t d = 0; d < 3; ++d) {
        transforms.push_back(0.0);
      }
    }
  }

  [[nodiscard]] static std::size_t side(std::size_t face) { return face % Cell::NumFaces; }
};

/// The quantities a face offers: those of its cell, and its displacement.
inline std::vector<DerivedSource> surfaceSources() {
  auto result = sources();
  for (const auto* name : {"u1", "u2", "u3"}) {
    result.push_back(DerivedSource{name, Representation::FaceNodal});
  }
  return result;
}

inline DerivedSurfaceGeometry surfaceGeometry(std::size_t degree, std::size_t refinement) {
  DerivedSurfaceGeometry geometry;
  for (std::size_t i = 0; i < refinement; ++i) {
    geometry.subcells = io::instance::geometry::subdivideMaps(
        geometry.subcells, io::instance::geometry::TriangleRefine4);
  }
  geometry.dataBase = io::instance::geometry::pointsTriangle(degree);
  geometry.dataOrder = degree;
  geometry.order = Cfg::ConvergenceOrder;
  return geometry;
}

/// A program bound to cells, as the volume output binds it: one block per quantity, the inverse
/// Jacobian per cell, and one column per output. With `faces`, bound to the faces as the free
/// surface output binds them, and run face by face with the matrices of the side of the face.
struct Evaluation {
  const DerivedProgram* derived;
  const Faces* faces;
  std::size_t numPoints;
  std::vector<double> values;
  DataTable table;
  expr::Binding binding;

  Evaluation(const DerivedProgram& program,
             Cells& cells,
             std::size_t simulation = 0,
             const Faces* faces = nullptr)
      : derived(&program), faces(faces), numPoints(cells.count * program.pointsPerElement()),
        values(program.program().outputs().size() * numPoints,
               std::numeric_limits<double>::quiet_NaN()),
        table(numPoints), binding(bindAll(program, cells, simulation)) {}

  expr::Binding bindAll(const DerivedProgram& program, Cells& cells, std::size_t simulation) {
    constexpr std::size_t Quantities = MaterialT::Quantities.size();
    program.bindMatrices(table);
    for (std::size_t b = 0; b < program.usedSources().size(); ++b) {
      const std::size_t source = program.usedSources()[b];
      const auto& block = program.program().blocks()[b];
      if (source < 2 * Quantities) {
        const std::size_t quantity = source % Quantities;
        const auto& storage = source < Quantities ? cells.dofs : cells.integrals;
        table.bindBlock<RealT>(block.name,
                               storage.data() + quantity * QuantityStride + simulation,
                               block.length,
                               CellStride,
                               Simulations);
      } else {
        REQUIRE(faces != nullptr);
        const std::size_t component = source - 2 * Quantities;
        table.bindBlock<RealT>(block.name,
                               faces->displacements.data() + component * ComponentStride +
                                   simulation,
                               block.length,
                               FaceStride,
                               Simulations);
      }
    }
    program.bindGeometry(table,
                         faces != nullptr ? faces->transforms.data() : cells.transforms.data());
    if (program.readsJacobian()) {
      for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t d = 0; d < 3; ++d) {
          table.bindCellView<double>("jinv" + std::to_string(k) + std::to_string(d),
                                     cells.jacobians.data(),
                                     program.pointsPerElement(),
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

  void run(expr::Kernel& kernel) const {
    if (faces == nullptr) {
      kernel.run(table);
      return;
    }
    const std::size_t pointsPerElement = derived->pointsPerElement();
    for (std::size_t face = 0; face < numPoints / pointsPerElement; ++face) {
      const auto matrices = derived->matrixBases(Faces::side(face));
      expr::KernelArgs args;
      args.matrices = matrices.data();
      args.matrixCount = matrices.size();
      args.first = face * pointsPerElement;
      args.count = pointsPerElement;
      kernel.run(args);
    }
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

inline std::unique_ptr<expr::Kernel> kernelFor(Evaluation& evaluation,
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
      target[i] = dataX[i] * grad[dir] + dataY[i] * grad[3 + dir] + dataZ[i] * grad[6 + dir];
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

/// The affine embedding of the reference triangle into side `side` of the reference tetrahedron,
/// as the free surface output of master built it.
inline numerical::AffineMap<2, 3> faceEmbedding(std::size_t side) {
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
  return numerical::AffineMap<2, 3>::fromVertices(vertices);
}

/// The free surface outputs as master computed them, through the generated kernels, kept
/// verbatim apart from the plumbing.
class HandWrittenSurface {
  public:
  HandWrittenSurface(const DerivedSurfaceGeometry& geometry, std::size_t degree)
      : degree_(degree), pointsPerSubcell_(geometry.dataBase.size()) {
    // volume basis -> face points, one table per side of the reference tetrahedron
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto embedding = faceEmbedding(f);
      std::vector<numerical::AffineMap<2, 3>> embedded;
      embedded.reserve(geometry.subcells.size());
      for (const auto& subcell : geometry.subcells) {
        embedded.emplace_back(embedding.compose(subcell));
      }
      projection::Spec spec;
      spec.target = geometry.target;
      proj_[f] = std::make_shared<projection::Table<2, 3, RealT>>(embedded,
                                                                  geometry.dataBase,
                                                                  geometry.dataOrder,
                                                                  stride<tensor::collvf<Cfg>>(),
                                                                  spec,
                                                                  1,
                                                                  Cfg::ConvergenceOrder);
    }
    // face nodes -> face points
    projection::Spec faceSpec;
    faceSpec.source = projection::Source::Nodal;
    faceSpec.target = geometry.target;
    faceSpec.nodalSet = projection::NodalSet::WarpBlend;
    projf_ = std::make_shared<projection::Table<2, 2, RealT>>(geometry.subcells,
                                                              geometry.dataBase,
                                                              geometry.dataOrder,
                                                              stride<tensor::collnf<Cfg>>(),
                                                              faceSpec,
                                                              1,
                                                              Cfg::ConvergenceOrder);
  }

  void value(double* target,
             const RealT* dofsSingleQuantity,
             std::size_t side,
             std::size_t subcell) const {
    constexpr auto Variant = configIdOf<Cfg>();
    runtime::kernel::projectBasisToVtkFaceFromVolume vtkproj{};
    memory::AlignedArray<RealT, Cfg::NumSimulations> simselect{};
    alignas(Alignment) std::array<RealT, MaxVtk2dPoints> alignedTarget{};
    simselect[0] = 1;
    vtkproj.simselect = runtime::init::simselect::view(Variant, simselect.data());
    vtkproj.qb = runtime::init::qb::view(Variant, dofsSingleQuantity);
    vtkproj.xf(degree_) = runtime::init::xf::view(Variant, degree_, alignedTarget.data());
    vtkproj.collvf(Cfg::ConvergenceOrder, degree_) = runtime::init::collvf::view(
        Variant, Cfg::ConvergenceOrder, degree_, (*proj_[side])(subcell, Cfg::ConvergenceOrder));
    vtkproj.execute(Variant, degree_);
    std::copy_n(alignedTarget.data(), pointsPerSubcell_, target);
  }

  void displacement(double* target,
                    const RealT* faceDisplacementVariable,
                    std::size_t subcell) const {
    constexpr auto Variant = configIdOf<Cfg>();
    runtime::kernel::projectNodalToVtkFace vtkproj{};
    memory::AlignedArray<RealT, Cfg::NumSimulations> simselect{};
    alignas(Alignment) std::array<RealT, MaxVtk2dPoints> alignedTarget{};
    simselect[0] = 1;
    vtkproj.simselect = runtime::init::simselect::view(Variant, simselect.data());
    vtkproj.pn = runtime::init::pn::view(Variant, faceDisplacementVariable);
    vtkproj.xf(degree_) = runtime::init::xf::view(Variant, degree_, alignedTarget.data());
    vtkproj.collnf(Cfg::ConvergenceOrder, degree_) = runtime::init::collnf::view(
        Variant, Cfg::ConvergenceOrder, degree_, (*projf_)(subcell, Cfg::ConvergenceOrder));
    vtkproj.execute(Variant, degree_);
    std::copy_n(alignedTarget.data(), pointsPerSubcell_, target);
  }

  static constexpr std::size_t MaxVtk2dPoints = tensor::vtk2d<Cfg>::Shape
      [(sizeof(tensor::vtk2d<Cfg>::Shape) / sizeof(tensor::vtk2d<Cfg>::Shape[0])) - 1][1];

  private:
  template <typename TensorT>
  [[nodiscard]] std::size_t stride() const {
    const auto index = TensorT::index(Cfg::ConvergenceOrder, degree_);
    return TensorT::Size[index] / TensorT::Shape[index][1];
  }

  std::size_t degree_;
  std::size_t pointsPerSubcell_;
  std::array<std::shared_ptr<projection::Table<2, 3, RealT>>, Cell::NumFaces> proj_;
  std::shared_ptr<projection::Table<2, 2, RealT>> projf_;
};

inline WaveFieldSelection fullSelection() {
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
  REQUIRE(derived.pointsPerElement() == geometry.subcells.size() * geometry.dataBase.size());
  Evaluation evaluation(derived, cells);
  reader::datafield::GridStore store;
  kernelFor(evaluation, expr::BackendKind::Interpreter, store)->run(evaluation.table);

  const HandWritten master(geometry, Degree);
  const std::size_t pointsPerSubcell = derived.pointsPerSubcell();
  std::vector<double> reference(pointsPerSubcell);
  const auto check = [&](const std::string& output, std::size_t cell, std::size_t subcell) {
    for (std::size_t i = 0; i < pointsPerSubcell; ++i) {
      const std::size_t point = cell * derived.pointsPerElement() + subcell * pointsPerSubcell + i;
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

// Products and sums through volatile, so that the reference rounds every operation on its own,
// like the interpreter, whatever this translation unit is compiled with.
inline double multiply(double a, double b) {
  const volatile double product = a * b;
  return product;
}
inline double add(double a, double b) {
  const volatile double sum = a + b;
  return sum;
}

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
    CHECK(
        bitwiseEqual(interpreted.values.data(), compiled.values.data(), interpreted.values.size()));
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
        const std::size_t point =
            cell * derived.pointsPerElement() + subcell * pointsPerSubcell + p;
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
  CHECK(bitwiseEqual(first.values.data(), second.values.data(), first.values.size()));
}

TEST_CASE("DerivedOutput: a maximum over time follows a hand-kept reference") {
  constexpr std::size_t Degree = 1;
  const auto geometry = refinedGeometry(Degree);
  Cells cells(4, 4);

  const std::string& v1 = MaterialT::Quantities[MaterialT::VelocityOffset];
  const std::string& v2 = MaterialT::Quantities[MaterialT::VelocityOffset + 1];
  const std::string& v3 = MaterialT::Quantities[MaterialT::VelocityOffset + 2];
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

  const std::string& v1 = MaterialT::Quantities[MaterialT::VelocityOffset];
  const std::string& v2 = MaterialT::Quantities[MaterialT::VelocityOffset + 1];
  const std::string& v3 = MaterialT::Quantities[MaterialT::VelocityOffset + 2];
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
        REQUIRE(bitwiseEqual(first.value(output, point), second.value(output, point)));
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
  CHECK(VolumePoints(geometry).referencePoints().size() == derived.pointsPerElement());

  CHECK(readsTimeIntegral("int_v1"));
  CHECK(readsTimeIntegral("dy_int_v1"));
  CHECK(readsTimeIntegral("int_v1_r2"));
  CHECK_FALSE(readsTimeIntegral("dy_v1"));
  CHECK_FALSE(readsTimeIntegral("v1"));
}

TEST_CASE("DerivedOutput: the built-in surface program reproduces the hand-written outputs") {
  // the unrefined default of the free surface output, and a refined higher-order one
  for (const auto& levels :
       {std::pair<std::size_t, std::size_t>{0, 0}, std::pair<std::size_t, std::size_t>{2, 1}}) {
    // (a structured binding cannot be captured where OpenMP is enabled)
    const std::size_t degree = levels.first;
    const std::size_t refinement = levels.second;
    CAPTURE(degree);
    const auto geometry = surfaceGeometry(degree, refinement);
    Faces faces(8, 7 + static_cast<unsigned>(degree));

    const std::vector<std::string> quantities(MaterialT::Quantities.begin(),
                                              MaterialT::Quantities.end());
    const DerivedProgram derived(
        surfaceProgram(quantities, std::vector<bool>(quantities.size(), true)),
        surfaceSources(),
        SurfacePoints(geometry));
    REQUIRE(derived.pointsPerElement() == geometry.subcells.size() * geometry.dataBase.size());
    Evaluation evaluation(derived, faces.cells, 0, &faces);
    reader::datafield::GridStore store;
    evaluation.run(*kernelFor(evaluation, expr::BackendKind::Interpreter, store));

    const HandWrittenSurface master(geometry, degree);
    const std::size_t pointsPerSubcell = derived.pointsPerSubcell();
    std::vector<double> reference(pointsPerSubcell);
    const auto check = [&](const std::string& output, std::size_t face, std::size_t subcell) {
      for (std::size_t i = 0; i < pointsPerSubcell; ++i) {
        const std::size_t point =
            face * derived.pointsPerElement() + subcell * pointsPerSubcell + i;
        REQUIRE(evaluation.value(output, point) ==
                doctest::Approx(reference[i]).epsilon(Tolerance));
      }
    };

    for (std::size_t face = 0; face < faces.cells.count; ++face) {
      const RealT* dofs = faces.cells.dofs.data() + face * CellStride;
      const RealT* displacement = faces.displacements.data() + face * FaceStride;
      for (std::size_t subcell = 0; subcell < geometry.subcells.size(); ++subcell) {
        for (std::size_t q = 0; q < MaterialT::Quantities.size(); ++q) {
          master.value(reference.data(), dofs + q * QuantityStride, Faces::side(face), subcell);
          check(MaterialT::Quantities[q], face, subcell);
        }
        for (std::size_t component = 0; component < 3; ++component) {
          master.displacement(
              reference.data(), displacement + component * ComponentStride, subcell);
          check("u" + std::to_string(component + 1), face, subcell);
        }
      }
    }
  }
}

TEST_CASE("DerivedOutput: a derivative on a face is the one in the cell at the same points") {
  // with point evaluation on both, a face output point is a cell output point like any other
  auto geometry = surfaceGeometry(2, 1);
  geometry.target = projection::Target::Interpolate;
  Faces faces(8, 11);

  const std::string& quantity = MaterialT::Quantities[MaterialT::VelocityOffset];
  const auto program = expr::compileSderivModule("out def a = dx_" + quantity + " * dz_int_" +
                                                 quantity + " + dy_" + quantity + "\n");
  const DerivedProgram onFaces(program, surfaceSources(), SurfacePoints(geometry));
  Evaluation surface(onFaces, faces.cells, 0, &faces);
  reader::datafield::GridStore store;
  surface.run(*kernelFor(surface, expr::BackendKind::Interpreter, store));

  const auto reference = SurfacePoints(geometry).referencePoints();
  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    DerivedGeometry inCell;
    inCell.order = Cfg::ConvergenceOrder;
    inCell.dataBase.clear();
    for (const auto& point : reference) {
      inCell.dataBase.push_back(faceEmbedding(side)({point[0], point[1]}));
    }
    const DerivedProgram inCells(program, sources(), inCell);
    Evaluation volume(inCells, faces.cells);
    volume.run(*kernelFor(volume, expr::BackendKind::Interpreter, store));
    for (std::size_t face = side; face < faces.cells.count; face += Cell::NumFaces) {
      for (std::size_t p = 0; p < reference.size(); ++p) {
        const std::size_t point = face * reference.size() + p;
        REQUIRE(surface.value("a", point) ==
                doctest::Approx(volume.value("a", point)).epsilon(1e-12));
      }
    }
  }
}

TEST_CASE("DerivedOutput: the coordinates are the affine map of the element at its points") {
  const auto program =
      expr::compileSderivModule("out def px = x\nout def py = y\nout def pz = z\n");
  reader::datafield::GridStore store;

  // the cells of the wave field
  const auto volumeGeometry = refinedGeometry(2);
  Cells cells(5, 8);
  const DerivedProgram inCells(program, sources(), volumeGeometry);
  CHECK(inCells.readsCoordinates());
  CHECK(inCells.program().inputs().empty());
  Evaluation volume(inCells, cells);
  volume.run(*kernelFor(volume, expr::BackendKind::Interpreter, store));
  const auto cellPoints = VolumePoints(volumeGeometry).referencePoints();
  for (std::size_t cell = 0; cell < cells.count; ++cell) {
    for (std::size_t p = 0; p < cellPoints.size(); ++p) {
      const auto expected = cells.shapes[cell].refToSpace(cellPoints[p]);
      const std::size_t point = cell * cellPoints.size() + p;
      REQUIRE(volume.value("px", point) == doctest::Approx(expected[0]).epsilon(1e-14));
      REQUIRE(volume.value("py", point) == doctest::Approx(expected[1]).epsilon(1e-14));
      REQUIRE(volume.value("pz", point) == doctest::Approx(expected[2]).epsilon(1e-14));
    }
  }

  // the faces of the free surface
  const auto surfaceGeometryRefined = surfaceGeometry(2, 1);
  Faces faces(8, 9);
  const DerivedProgram onFaces(program, surfaceSources(), SurfacePoints(surfaceGeometryRefined));
  Evaluation surface(onFaces, faces.cells, 0, &faces);
  surface.run(*kernelFor(surface, expr::BackendKind::Interpreter, store));
  const auto facePoints = SurfacePoints(surfaceGeometryRefined).referencePoints();
  for (std::size_t face = 0; face < faces.cells.count; ++face) {
    for (std::size_t p = 0; p < facePoints.size(); ++p) {
      const auto expected = faces.shapes[face].refToSpace(
          seissol::geometry::FaceTransform::FaceVectorT(facePoints[p][0], facePoints[p][1]));
      const std::size_t point = face * facePoints.size() + p;
      REQUIRE(surface.value("px", point) == doctest::Approx(expected(0)).epsilon(1e-14));
      REQUIRE(surface.value("py", point) == doctest::Approx(expected(1)).epsilon(1e-14));
      REQUIRE(surface.value("pz", point) == doctest::Approx(expected(2)).epsilon(1e-14));
    }
  }
}

TEST_CASE("DerivedOutput: the displacement is read on faces, and without derivatives") {
  const auto reads = [](const std::string& expression) {
    return expr::compileSderivModule("out def a = " + expression + "\n");
  };
  const auto faceGeometry = surfaceGeometry(1, 0);
  // a cell has no face nodes
  CHECK_THROWS_AS(DerivedProgram(reads("u1"), surfaceSources(), refinedGeometry(1)),
                  std::invalid_argument);
  // the face nodes give no derivatives across the face
  CHECK_THROWS_AS(DerivedProgram(reads("dz_u2"), surfaceSources(), SurfacePoints(faceGeometry)),
                  std::invalid_argument);
  CHECK_THROWS_AS(DerivedProgram(reads("u3_r0"), surfaceSources(), SurfacePoints(faceGeometry)),
                  std::invalid_argument);
  // the derivatives of the quantities of the cell are there
  const DerivedProgram derived(
      reads("u1 + dx_" + MaterialT::Quantities[0]), surfaceSources(), SurfacePoints(faceGeometry));
  CHECK(derived.readsJacobian());
  CHECK(derived.matrixBases(3).size() == derived.program().matrices().size());
  CHECK(derived.matrixBases(0) != derived.matrixBases(3));
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
