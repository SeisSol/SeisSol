// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DerivedOutput.h"

#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Expr/Rewrite.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "IO/Instance/Geometry/Points.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "Numerical/Projection.h"
#include "Parallel/OpenMP.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"
#include "SeisSol.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/MultipleSimulations.h"
#include "Solver/Simulator.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <ios>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

namespace {

namespace projection = seissol::numerical::projection;
using expr::Arena;
using expr::Fn;
using expr::NodeId;
using reader::scripting::DataTable;
using reader::scripting::DataType;
using reader::scripting::Direction;

constexpr std::array<char, 3> AxisNames{'x', 'y', 'z'};

/// The prefix of a quantity that names its time integral.
const std::string IntegralPrefix = "int_";

/// The affine map of an element is contracted against the matrix of this name, whose row p is
/// (1, xi_0, xi_1, xi_2) for output point p, and the blocks of this prefix, one per axis.
const std::string GeometryMatrix = "geometry";
const std::string GeometryPrefix = "geometry_";
constexpr std::size_t GeometryColumns = 4;
constexpr std::size_t TransformSize = GeometryColumns * 3;

std::string jacobianName(std::size_t row, std::size_t column) {
  return "jinv" + std::to_string(row) + std::to_string(column);
}

std::string physicalName(std::size_t axis, const std::string& quantity) {
  return std::string("d") + AxisNames[axis] + "_" + quantity;
}

/// What an input channel of a derived program reads.
struct Channel {
  enum class Kind : std::uint8_t {
    Value,
    Reference,
    Physical,
    Jacobian,
    Coordinate,
    Time,
    TimeStep
  };
  Kind kind{Kind::Value};
  std::size_t source{0};
  /// The reference direction (Reference), the axis (Physical, Coordinate), or the entry
  /// row * 3 + column (Jacobian).
  std::size_t direction{0};

  [[nodiscard]] bool readsSource() const {
    return kind == Kind::Value || kind == Kind::Reference || kind == Kind::Physical;
  }
};

std::optional<std::size_t> findSource(const std::vector<DerivedSource>& sources,
                                      const std::string& name) {
  for (std::size_t i = 0; i < sources.size(); ++i) {
    if (sources[i].name == name) {
      return i;
    }
  }
  return std::nullopt;
}

std::optional<std::size_t> axisOf(char letter) {
  for (std::size_t axis = 0; axis < AxisNames.size(); ++axis) {
    if (AxisNames[axis] == letter) {
      return axis;
    }
  }
  return std::nullopt;
}

std::optional<std::size_t> directionOf(char digit) {
  if (digit >= '0' && digit <= '2') {
    return static_cast<std::size_t>(digit - '0');
  }
  return std::nullopt;
}

/// The meaning of the channel `name`, or nothing if it is outside the vocabulary. A source name
/// is looked up first, so a quantity of the material is never shadowed by a derived spelling.
std::optional<Channel> classify(const std::string& name,
                                const std::vector<DerivedSource>& sources) {
  using Kind = Channel::Kind;
  if (const auto source = findSource(sources, name)) {
    return Channel{Kind::Value, *source, 0};
  }
  if (name.size() == 1) {
    if (const auto axis = axisOf(name[0])) {
      return Channel{Kind::Coordinate, 0, *axis};
    }
  }
  if (name == "t") {
    return Channel{Kind::Time, 0, 0};
  }
  if (name == "dt") {
    return Channel{Kind::TimeStep, 0, 0};
  }
  if (name.size() == 6 && name.compare(0, 4, "jinv") == 0) {
    const auto row = directionOf(name[4]);
    const auto column = directionOf(name[5]);
    if (row.has_value() && column.has_value()) {
      return Channel{Kind::Jacobian, 0, *row * 3 + *column};
    }
  }
  if (name.size() > 3 && name.compare(name.size() - 3, 2, "_r") == 0) {
    const auto direction = directionOf(name.back());
    const auto source = findSource(sources, name.substr(0, name.size() - 3));
    if (direction.has_value() && source.has_value()) {
      return Channel{Kind::Reference, *source, *direction};
    }
  }
  if (name.size() > 3 && name[0] == 'd' && name[2] == '_') {
    const auto axis = axisOf(name[1]);
    const auto source = findSource(sources, name.substr(3));
    if (axis.has_value() && source.has_value()) {
      return Channel{Kind::Physical, *source, *axis};
    }
  }
  return std::nullopt;
}

std::string unknownChannel(const std::string& name, const std::vector<DerivedSource>& sources) {
  std::ostringstream message;
  message << "derived output: the program reads '" << name
          << "', which is neither a quantity nor derived from one. Available are the quantities";
  for (const auto& source : sources) {
    message << " " << source.name;
  }
  message << "; for a quantity q of the cell also q_r0, q_r1, q_r2 (reference derivatives) and "
             "dx_q, dy_q, dz_q; and jinv00 ... jinv22, x, y, z, t, dt.";
  return message.str();
}

/// Matrices are keyed by the representation of their source and the derivative: 0 is the value,
/// 1 + k the derivative along reference direction k.
constexpr std::size_t DerivativeKinds = 4;
constexpr std::size_t Representations = 3;
constexpr std::size_t MatrixKinds = Representations * DerivativeKinds;

std::size_t matrixKey(Representation representation, std::optional<std::size_t> derivative) {
  return static_cast<std::size_t>(representation) * DerivativeKinds +
         (derivative.has_value() ? 1 + *derivative : 0);
}

Representation representationOf(std::size_t key) {
  return static_cast<Representation>(key / DerivativeKinds);
}

std::optional<std::size_t> derivativeOf(std::size_t key) {
  if (key % DerivativeKinds == 0) {
    return std::nullopt;
  }
  return key % DerivativeKinds - 1;
}

std::string matrixName(std::size_t key) {
  const std::array<std::string, Representations> names{"modal", "nodal", "face_nodal"};
  const std::string& base = names[key / DerivativeKinds];
  const auto derivative = derivativeOf(key);
  return derivative.has_value() ? base + "_r" + std::to_string(*derivative) : base;
}

projection::Spec specOf(std::size_t order,
                        projection::Source source,
                        projection::Target target,
                        projection::NodalSet nodalSet,
                        std::optional<std::size_t> derivative) {
  projection::Spec spec;
  spec.order = order;
  spec.source = source;
  spec.target = target;
  spec.nodalSet = nodalSet;
  spec.derivative = derivative;
  return spec;
}

/// The projection matrices of all subcells stacked into one: row s * pointsPerSubcell + p is
/// output point p of subcell s, entry (row, m) at row + m * rows.
template <std::size_t From, std::size_t To>
std::vector<double> stackedMatrix(const std::vector<std::array<double, From>>& dataBase,
                                  std::size_t dataOrder,
                                  const std::vector<numerical::AffineMap<From, To>>& subcells,
                                  const projection::Spec& spec,
                                  std::size_t cols) {
  const std::size_t pointsPerSubcell = dataBase.size();
  const std::size_t rows = subcells.size() * pointsPerSubcell;
  std::vector<double> values(rows * cols);
  for (std::size_t subcell = 0; subcell < subcells.size(); ++subcell) {
    const auto local = projection::build<From, To>(dataBase, dataOrder, subcells[subcell], spec);
    if (local.rows() != pointsPerSubcell || local.cols() != cols) {
      throw std::logic_error("derived output: a projection matrix has an unexpected shape");
    }
    for (std::size_t point = 0; point < pointsPerSubcell; ++point) {
      const std::size_t row = subcell * pointsPerSubcell + point;
      for (std::size_t m = 0; m < cols; ++m) {
        values[row + m * rows] = local(point, m);
      }
    }
  }
  return values;
}

/// The affine embedding of the reference triangle into side `side` of the reference tetrahedron.
numerical::AffineMap<2, 3> faceEmbedding(std::size_t side) {
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

bool startsWith(const std::string& text, const std::string& prefix) {
  return text.compare(0, prefix.size(), prefix) == 0;
}

} // namespace

VolumePoints::VolumePoints(DerivedGeometry geometry) : geometry_(std::move(geometry)) {}

std::size_t VolumePoints::pointsPerElement() const {
  return geometry_.subcells.size() * geometry_.dataBase.size();
}

std::size_t VolumePoints::pointsPerSubcell() const { return geometry_.dataBase.size(); }

std::size_t VolumePoints::columns(Representation representation) const {
  switch (representation) {
  case Representation::Modal:
    return projection::modalSize(3, geometry_.order);
  case Representation::Nodal:
    return projection::nodalSize(3, geometry_.order, geometry_.nodalSet);
  case Representation::FaceNodal:
    return 0;
  }
  return 0;
}

bool VolumePoints::derivatives(Representation representation) const {
  return representation != Representation::FaceNodal;
}

std::vector<double> VolumePoints::matrix(std::size_t /*variant*/,
                                         Representation representation,
                                         std::optional<std::size_t> derivative) const {
  if (representation == Representation::FaceNodal) {
    throw std::logic_error("derived output: a cell has no face nodes");
  }
  const auto spec = specOf(geometry_.order,
                           representation == Representation::Nodal ? projection::Source::Nodal
                                                                   : projection::Source::Modal,
                           geometry_.target,
                           geometry_.nodalSet,
                           derivative);
  return stackedMatrix<3, 3>(
      geometry_.dataBase, geometry_.dataOrder, geometry_.subcells, spec, columns(representation));
}

std::vector<std::array<double, 3>> VolumePoints::referencePoints() const {
  std::vector<std::array<double, 3>> points;
  points.reserve(pointsPerElement());
  for (const auto& subcell : geometry_.subcells) {
    for (const auto& point : geometry_.dataBase) {
      points.push_back(subcell(point));
    }
  }
  return points;
}

SurfacePoints::SurfacePoints(DerivedSurfaceGeometry geometry) : geometry_(std::move(geometry)) {}

std::size_t SurfacePoints::pointsPerElement() const {
  return geometry_.subcells.size() * geometry_.dataBase.size();
}

std::size_t SurfacePoints::pointsPerSubcell() const { return geometry_.dataBase.size(); }

std::size_t SurfacePoints::variants() const { return Cell::NumFaces; }

std::size_t SurfacePoints::columns(Representation representation) const {
  switch (representation) {
  case Representation::Modal:
    return projection::modalSize(3, geometry_.order);
  case Representation::Nodal:
    return projection::nodalSize(3, geometry_.order, geometry_.nodalSet);
  case Representation::FaceNodal:
    return projection::nodalSize(2, geometry_.order, projection::NodalSet::WarpBlend);
  }
  return 0;
}

bool SurfacePoints::derivatives(Representation representation) const {
  // the face nodes would only give the derivatives along the face
  return representation != Representation::FaceNodal;
}

std::vector<double> SurfacePoints::matrix(std::size_t variant,
                                          Representation representation,
                                          std::optional<std::size_t> derivative) const {
  if (representation == Representation::FaceNodal) {
    if (derivative.has_value()) {
      throw std::logic_error("derived output: the face nodes have no derivatives across the face");
    }
    // the face displacement is stored at the nodes of the face, which are always warp&blend
    const auto spec = specOf(geometry_.order,
                             projection::Source::Nodal,
                             geometry_.target,
                             projection::NodalSet::WarpBlend,
                             std::nullopt);
    return stackedMatrix<2, 2>(
        geometry_.dataBase, geometry_.dataOrder, geometry_.subcells, spec, columns(representation));
  }
  const auto embedding = faceEmbedding(variant);
  std::vector<numerical::AffineMap<2, 3>> embedded;
  embedded.reserve(geometry_.subcells.size());
  for (const auto& subcell : geometry_.subcells) {
    embedded.push_back(embedding.compose(subcell));
  }
  const auto spec = specOf(geometry_.order,
                           representation == Representation::Nodal ? projection::Source::Nodal
                                                                   : projection::Source::Modal,
                           geometry_.target,
                           geometry_.nodalSet,
                           derivative);
  return stackedMatrix<2, 3>(
      geometry_.dataBase, geometry_.dataOrder, embedded, spec, columns(representation));
}

std::vector<std::array<double, 3>> SurfacePoints::referencePoints() const {
  std::vector<std::array<double, 3>> points;
  points.reserve(pointsPerElement());
  for (const auto& subcell : geometry_.subcells) {
    for (const auto& point : geometry_.dataBase) {
      const auto onFace = subcell(point);
      points.push_back({onFace[0], onFace[1], 0.0});
    }
  }
  return points;
}

DerivedProgram::DerivedProgram(expr::Program program,
                               const std::vector<DerivedSource>& sources,
                               const DerivedGeometry& geometry)
    : DerivedProgram(std::move(program), sources, VolumePoints(geometry)) {}

DerivedProgram::DerivedProgram(expr::Program program,
                               const std::vector<DerivedSource>& sources,
                               const DerivedPoints& points)
    : program_(std::move(program)), pointsPerSubcell_(points.pointsPerSubcell()),
      pointsPerElement_(points.pointsPerElement()) {
  if (pointsPerElement_ == 0) {
    throw std::invalid_argument("derived output: there are no output points");
  }
  for (const auto& output : program_.outputs()) {
    if (output.type != DataType::F64) {
      throw std::invalid_argument("derived output: the output " + output.name +
                                  " is not a 64-bit floating-point value");
    }
  }

  // The channels the program reads, and what they need.
  using Kind = Channel::Kind;
  std::vector<std::pair<std::string, Channel>> channels;
  std::vector<bool> sourceUsed(sources.size(), false);
  std::array<bool, MatrixKinds> matrixUsed{};
  for (const auto& input : program_.inputs()) {
    const auto channel = classify(input.name, sources);
    if (!channel.has_value()) {
      throw std::invalid_argument(unknownChannel(input.name, sources));
    }
    const auto representation =
        channel->readsSource() ? sources[channel->source].representation : Representation::Modal;
    if (channel->readsSource()) {
      const auto& name = sources[channel->source].name;
      if (points.columns(representation) == 0) {
        throw std::invalid_argument("derived output: the program reads '" + input.name + "', but " +
                                    name + " is not available at these output points");
      }
      if (channel->kind != Kind::Value && !points.derivatives(representation)) {
        throw std::invalid_argument("derived output: the program reads '" + input.name + "', but " +
                                    name + " has no derivatives at these output points");
      }
    }
    switch (channel->kind) {
    case Kind::Value:
      sourceUsed[channel->source] = true;
      matrixUsed[matrixKey(representation, std::nullopt)] = true;
      break;
    case Kind::Reference:
      sourceUsed[channel->source] = true;
      matrixUsed[matrixKey(representation, channel->direction)] = true;
      break;
    case Kind::Physical:
      sourceUsed[channel->source] = true;
      for (std::size_t k = 0; k < 3; ++k) {
        matrixUsed[matrixKey(representation, k)] = true;
      }
      readsJacobian_ = true;
      break;
    case Kind::Jacobian:
      readsJacobian_ = true;
      break;
    case Kind::Coordinate:
      readsCoordinates_ = true;
      break;
    case Kind::Time:
      readsTime_ = true;
      break;
    case Kind::TimeStep:
      readsTimeStep_ = true;
      break;
    }
    channels.emplace_back(input.name, *channel);
  }

  // One matrix per representation and derivative, with values per variant; one block per
  // quantity.
  const std::size_t variants = points.variants();
  std::array<expr::MatrixId, MatrixKinds> matrixIds{};
  matrixIds.fill(expr::NoMatrix);
  for (std::size_t key = 0; key < MatrixKinds; ++key) {
    if (!matrixUsed[key]) {
      continue;
    }
    const std::size_t cols = points.columns(representationOf(key));
    matrixIds[key] = program_.internMatrix(
        matrixName(key), expr::MatrixShape{pointsPerElement_, cols, pointsPerElement_});
    Matrix matrix{matrixName(key), {}, cols};
    for (std::size_t variant = 0; variant < variants; ++variant) {
      matrix.values.push_back(points.matrix(variant, representationOf(key), derivativeOf(key)));
    }
    matrices_.push_back(std::move(matrix));
  }
  std::vector<expr::BlockId> blockIds(sources.size(), expr::NoBlock);
  for (std::size_t i = 0; i < sources.size(); ++i) {
    if (sourceUsed[i]) {
      blockIds[i] =
          program_.internBlock(sources[i].name, points.columns(sources[i].representation));
      usedSources_.push_back(i);
    }
  }

  // The coordinates are the affine map of the element at its reference points.
  expr::MatrixId geometryMatrix = expr::NoMatrix;
  std::array<expr::BlockId, 3> geometryBlocks{expr::NoBlock, expr::NoBlock, expr::NoBlock};
  if (readsCoordinates_) {
    geometryMatrix = program_.internMatrix(
        GeometryMatrix, expr::MatrixShape{pointsPerElement_, GeometryColumns, pointsPerElement_});
    const auto reference = points.referencePoints();
    std::vector<double> values(pointsPerElement_ * GeometryColumns);
    for (std::size_t point = 0; point < pointsPerElement_; ++point) {
      values[point] = 1.0;
      for (std::size_t direction = 0; direction < 3; ++direction) {
        values[point + (1 + direction) * pointsPerElement_] = reference[point][direction];
      }
    }
    matrices_.push_back(Matrix{
        GeometryMatrix, std::vector<std::vector<double>>(variants, values), GeometryColumns});
    for (std::size_t axis = 0; axis < 3; ++axis) {
      geometryBlocks[axis] =
          program_.internBlock(GeometryPrefix + AxisNames[axis], GeometryColumns);
    }
  }

  // Values, reference derivatives and coordinates become contractions; a physical derivative
  // becomes the chain rule over the reference derivatives, summed in the order of the reference
  // directions and with the derivative on the left, as the hand-written outputs did.
  std::map<std::string, expr::ChannelBuilder> builders;
  for (const auto& [name, channel] : channels) {
    const auto block = channel.readsSource() ? blockIds[channel.source] : expr::NoBlock;
    const auto representation =
        channel.readsSource() ? sources[channel.source].representation : Representation::Modal;
    switch (channel.kind) {
    case Kind::Value:
    case Kind::Reference: {
      const auto matrix = matrixIds[matrixKey(representation,
                                              channel.kind == Kind::Value
                                                  ? std::nullopt
                                                  : std::optional<std::size_t>(channel.direction))];
      builders.emplace(name,
                       [matrix, block](Arena& arena) { return arena.contract(matrix, block); });
      break;
    }
    case Kind::Physical: {
      const std::array<expr::MatrixId, 3> matrices{matrixIds[matrixKey(representation, 0)],
                                                   matrixIds[matrixKey(representation, 1)],
                                                   matrixIds[matrixKey(representation, 2)]};
      const std::size_t axis = channel.direction;
      builders.emplace(name, [matrices, block, axis](Arena& arena) {
        NodeId sum = expr::NoNode;
        for (std::size_t k = 0; k < 3; ++k) {
          const NodeId reference = arena.contract(matrices[k], block);
          const NodeId jacobian = arena.field(jacobianName(k, axis));
          const NodeId term = arena.pw(Fn::Mul, reference, jacobian);
          sum = sum == expr::NoNode ? term : arena.pw(Fn::Add, sum, term);
        }
        return sum;
      });
      break;
    }
    case Kind::Coordinate: {
      const auto axisBlock = geometryBlocks[channel.direction];
      builders.emplace(name, [geometryMatrix, axisBlock](Arena& arena) {
        return arena.contract(geometryMatrix, axisBlock);
      });
      break;
    }
    default:
      break;
    }
  }
  if (!builders.empty()) {
    expr::substituteChannels(program_, builders);
  }
}

void DerivedProgram::bindMatrices(DataTable& table) const {
  for (const auto& matrix : matrices_) {
    table.bindMatrix<double>(matrix.name,
                             matrix.values.front().data(),
                             pointsPerElement_,
                             matrix.cols,
                             pointsPerElement_);
  }
}

std::vector<const void*> DerivedProgram::matrixBases(std::size_t variant) const {
  std::vector<const void*> bases;
  bases.reserve(program_.matrices().size());
  for (const auto& spec : program_.matrices()) {
    const auto found = std::find_if(matrices_.begin(), matrices_.end(), [&](const Matrix& matrix) {
      return matrix.name == spec.name;
    });
    if (found == matrices_.end()) {
      throw std::logic_error("derived output: the program contracts against an unknown matrix " +
                             spec.name);
    }
    bases.push_back(found->values.at(variant).data());
  }
  return bases;
}

void DerivedProgram::bindGeometry(DataTable& table, const double* transforms) const {
  if (!readsCoordinates_) {
    return;
  }
  for (std::size_t axis = 0; axis < 3; ++axis) {
    table.bindBlock<double>(
        GeometryPrefix + AxisNames[axis], transforms + axis, GeometryColumns, TransformSize, 3);
  }
}

bool readsTimeIntegral(const std::string& name) {
  if (startsWith(name, IntegralPrefix)) {
    return true;
  }
  return name.size() > 3 && name[0] == 'd' && axisOf(name[1]).has_value() && name[2] == '_' &&
         name.compare(3, IntegralPrefix.size(), IntegralPrefix) == 0;
}

expr::Program waveFieldProgram(const WaveFieldSelection& selection) {
  expr::Program program;
  Arena& arena = program.arena();
  const auto output = [&](const std::string& name, NodeId root) {
    program.addOutput(name, DataType::F64, root);
  };
  const auto velocity = [&](std::size_t i) {
    return selection.quantities.at(selection.velocityOffset + i);
  };

  // The inputs in the order of first use, which is what the binding sees.
  std::vector<std::string> inputs;
  const auto read = [&](const std::string& name) {
    if (std::find(inputs.begin(), inputs.end(), name) == inputs.end()) {
      inputs.push_back(name);
    }
    return arena.field(name);
  };

  for (std::size_t q = 0; q < selection.quantities.size(); ++q) {
    if (q < selection.outputMask.size() && selection.outputMask[q]) {
      output(selection.quantities[q], read(selection.quantities[q]));
    }
    if (q < selection.integrationMask.size() && selection.integrationMask[q]) {
      output("int-" + selection.quantities[q], read(IntegralPrefix + selection.quantities[q]));
    }
  }

  using Index = std::tuple<std::string, std::size_t, std::size_t>;
  if (selection.strain) {
    // (d_j v_i + d_i v_j) / 2 of the time integral of the velocity; d_i v_i on the diagonal
    for (const auto& [name, i, j] : {Index{"xx", 0, 0},
                                     Index{"yy", 1, 1},
                                     Index{"zz", 2, 2},
                                     Index{"xy", 0, 1},
                                     Index{"yz", 1, 2},
                                     Index{"xz", 0, 2}}) {
      const NodeId first = read(physicalName(j, IntegralPrefix + velocity(i)));
      if (i == j) {
        output("eps" + name, first);
      } else {
        const NodeId second = read(physicalName(i, IntegralPrefix + velocity(j)));
        output("eps" + name, arena.pw(Fn::Div, arena.pw(Fn::Add, first, second), arena.konst(2.0)));
      }
    }
  }

  if (selection.rotation) {
    // d_j v_i - d_i v_j
    for (const auto& [name, i, j] : {Index{"1", 2, 1}, Index{"2", 0, 2}, Index{"3", 1, 0}}) {
      const NodeId first = read(physicalName(j, velocity(i)));
      const NodeId second = read(physicalName(i, velocity(j)));
      output("rot" + name, arena.pw(Fn::Sub, first, second));
    }
  }

  for (std::size_t p = 0; p < selection.plasticQuantities.size(); ++p) {
    if (p < selection.plasticityMask.size() && selection.plasticityMask[p]) {
      output(selection.plasticQuantities[p], read(selection.plasticQuantities[p]));
    }
  }

  for (const auto& name : inputs) {
    program.addInput(name, DataType::F64);
  }
  expr::validate(program);
  return program;
}

expr::Program surfaceProgram(const std::vector<std::string>& quantities,
                             const std::vector<bool>& outputMask) {
  std::vector<std::string> names;
  for (std::size_t q = 0; q < quantities.size(); ++q) {
    if (q < outputMask.size() && outputMask[q]) {
      names.push_back(quantities[q]);
    }
  }
  for (const auto* component : {"u1", "u2", "u3"}) {
    names.emplace_back(component);
  }

  expr::Program program;
  for (const auto& name : names) {
    program.addOutput(name, DataType::F64, program.arena().field(name));
  }
  for (const auto& name : names) {
    program.addInput(name, DataType::F64);
  }
  expr::validate(program);
  return program;
}

DerivedGeometry waveFieldGeometry(const parameters::WaveFieldOutputParameters& parameters) {
  using parameters::VolumeRefinement;
  namespace geometry = io::instance::geometry;
  DerivedGeometry result;
  result.dataOrder = static_cast<std::size_t>(std::max(0, parameters.vtkorder));
  result.dataBase = geometry::pointsTetrahedron(result.dataOrder);
  result.subcells = geometry::unrefined<3>();
  if (parameters.refinement == VolumeRefinement::Refine4) {
    result.subcells = geometry::subdivideMaps(result.subcells, geometry::TetrahedronRefine4);
  }
  if (parameters.refinement == VolumeRefinement::Refine8) {
    result.subcells = geometry::subdivideMaps(result.subcells, geometry::TetrahedronRefine8);
  }
  if (parameters.refinement == VolumeRefinement::Refine32) {
    // the edge division has to come first; the legacy refinement::DivideTetrahedronBy32
    // subdivided by 8 and then split each of those subcells by its center point, which (as
    // subdivideMaps enumerates input-major) is the ordering 4*i + j the output cells had
    result.subcells = geometry::subdivideMaps(result.subcells, geometry::TetrahedronRefine8);
    result.subcells = geometry::subdivideMaps(result.subcells, geometry::TetrahedronRefine4);
  }
  result.target = parameters.projection == parameters::ProjectionMethod::L2
                      ? projection::Target::Project
                      : projection::Target::Interpolate;
  return result;
}

DerivedSurfaceGeometry surfaceGeometry(const parameters::FreeSurfaceOutputParameters& parameters) {
  namespace geometry = io::instance::geometry;
  DerivedSurfaceGeometry result;
  result.dataOrder = static_cast<std::size_t>(std::max(0, parameters.vtkorder));
  result.dataBase = geometry::pointsTriangle(result.dataOrder);
  result.subcells = geometry::unrefined<2>();
  for (unsigned i = 0; i < parameters.refinement; ++i) {
    result.subcells = geometry::subdivideMaps(result.subcells, geometry::TriangleRefine4);
  }
  result.target = parameters.projection == parameters::ProjectionMethod::L2
                      ? projection::Target::Project
                      : projection::Target::Interpolate;
  return result;
}

std::string DerivedStateLayout::checkpointName() const {
  std::ostringstream name;
  name << "derivedstate-" << std::hex << std::setw(16) << std::setfill('0') << fingerprint;
  return name.str();
}

DerivedStateLayout derivedStateLayout(seissol::SeisSol& seissolInstance, DerivedOutputKind kind) {
  const auto& output = seissolInstance.parameters().output;
  std::string script;
  std::size_t pointsPerElement = 0;
  std::size_t dataOrder = 0;
  std::size_t subcells = 0;
  projection::Target target{};
  if (kind == DerivedOutputKind::WaveField) {
    if (!output.waveFieldParameters.enabled) {
      return {};
    }
    script = output.waveFieldParameters.script;
    const auto geometry = waveFieldGeometry(output.waveFieldParameters);
    pointsPerElement = geometry.subcells.size() * geometry.dataBase.size();
    dataOrder = geometry.dataOrder;
    subcells = geometry.subcells.size();
    target = geometry.target;
  } else {
    if (!output.freeSurfaceParameters.enabled) {
      return {};
    }
    script = output.freeSurfaceParameters.script;
    const auto geometry = surfaceGeometry(output.freeSurfaceParameters);
    pointsPerElement = geometry.subcells.size() * geometry.dataBase.size();
    dataOrder = geometry.dataOrder;
    subcells = geometry.subcells.size();
    target = geometry.target;
  }
  if (script.empty()) {
    return {};
  }
  const auto program = loadDerivedProgram(script);
  if (program.state().empty()) {
    return {};
  }

  DerivedStateLayout layout;
  for (const auto& state : program.state()) {
    layout.states.push_back(state.name);
    layout.initial.push_back(state.initial);
  }
  layout.pointsPerElement = pointsPerElement;
  for (const auto config : seissolInstance.parameters().model.configs()) {
    dispatchConfig(config, [&](auto cfg) {
      layout.simulations = std::max<std::size_t>(layout.simulations, decltype(cfg)::NumSimulations);
    });
  }

  std::uint64_t hash = 0xcbf29ce484222325ULL;
  const auto mix = [&hash](std::uint64_t value) {
    hash ^= value;
    hash *= 0x100000001b3ULL;
  };
  mix(static_cast<std::uint64_t>(kind));
  for (const auto& name : layout.states) {
    for (const char c : name) {
      mix(static_cast<unsigned char>(c));
    }
    mix(0);
  }
  mix(layout.pointsPerElement);
  mix(layout.simulations);
  mix(dataOrder);
  mix(subcells);
  mix(static_cast<std::uint64_t>(target));
  layout.fingerprint = hash;
  return layout;
}

namespace {

/// Sets the state of every element of the storage to its initial values.
template <typename VariableT, typename StorageT>
void initializeState(StorageT& storage, const DerivedStateLayout& layout) {
  const std::size_t size = layout.size();
  for (auto& layer : storage.leaves(LayerMask(Ghost))) {
    double* state = layer.template var<VariableT>();
    if (state == nullptr) {
      continue;
    }
#pragma omp parallel for schedule(static)
    for (std::size_t element = 0; element < layer.size(); ++element) {
      for (std::size_t s = 0; s < layout.states.size(); ++s) {
        for (std::size_t point = 0; point < layout.pointsPerElement; ++point) {
          for (std::size_t sim = 0; sim < layout.simulations; ++sim) {
            state[element * size + layout.offset(s, point, sim)] = layout.initial[s];
          }
        }
      }
    }
  }
}

} // namespace

void initializeDerivedState(seissol::SeisSol& seissolInstance) {
  auto& memoryManager = seissolInstance.memoryManager();
  const auto waveField = derivedStateLayout(seissolInstance, DerivedOutputKind::WaveField);
  if (waveField.size() > 0) {
    initializeState<LTS::DerivedState>(memoryManager.ltsStorage(), waveField);
  }
  const auto surface = derivedStateLayout(seissolInstance, DerivedOutputKind::Surface);
  if (surface.size() > 0) {
    initializeState<SurfaceLTS::DerivedState>(memoryManager.surfaceStorage(), surface);
  }
}

void registerDerivedStateCheckpoints(io::instance::checkpoint::CheckpointManager& checkpoint,
                                     seissol::SeisSol& seissolInstance) {
  auto& memoryManager = seissolInstance.memoryManager();
  const auto waveField = derivedStateLayout(seissolInstance, DerivedOutputKind::WaveField);
  if (waveField.size() > 0) {
    auto& storage = memoryManager.ltsStorage();
    checkpoint.registerOptionalArray(
        waveField.checkpointName(), storage, storage.var<LTS::DerivedState>(), waveField.size());
  }
  const auto surface = derivedStateLayout(seissolInstance, DerivedOutputKind::Surface);
  if (surface.size() > 0) {
    auto& storage = memoryManager.surfaceStorage();
    checkpoint.registerOptionalArray(
        surface.checkpointName(), storage, storage.var<SurfaceLTS::DerivedState>(), surface.size());
  }
}

expr::Program loadDerivedProgram(const std::string& path) {
  std::string reason;
  auto program = reader::scripting::buildProgram(path, &reason);
  if (!program.has_value()) {
    // there is no interpreted fallback for a program that reads contractions
    logError() << "derived output: the program" << path << "cannot be used --" << reason
               << "; a derived output is an sderiv module or a Lua model that traces.";
  }
  return std::move(*program);
}

bool derivedOutputsReadIntegrals(seissol::SeisSol& seissolInstance) {
  const auto& output = seissolInstance.parameters().output;
  const auto scriptReadsIntegrals = [](const std::string& path) {
    const auto program = loadDerivedProgram(path);
    return std::any_of(program.inputs().begin(),
                       program.inputs().end(),
                       [](const expr::VarSpec& input) { return readsTimeIntegral(input.name); });
  };

  const auto& waveField = output.waveFieldParameters;
  const bool masked = std::any_of(waveField.integrationMask.begin(),
                                  waveField.integrationMask.end(),
                                  [](bool integrate) { return integrate; });
  if (masked || waveField.computeStrain) {
    return true;
  }
  if (waveField.enabled && !waveField.script.empty() && scriptReadsIntegrals(waveField.script)) {
    return true;
  }
  const auto& surface = output.freeSurfaceParameters;
  return surface.enabled && !surface.script.empty() && scriptReadsIntegrals(surface.script);
}

namespace {

/// Which array an element reads a source from, and where in it.
enum class Storage : std::uint8_t { Dofs, Integrals, PlasticStrain, FaceDisplacement };

struct SourceLocation {
  Storage storage{Storage::Dofs};
  /// In reals, from the start of the array of the cell (or face).
  std::size_t offset{0};
};

/// A written element of the configuration -- a cell, or a face of one -- where it is stored, and
/// which one it is.
struct WrittenElement {
  std::size_t color{0};
  /// The side of a face; 0 for a cell.
  std::size_t variant{0};
  /// The cell in its layer of the cell storage, and the face in its layer of the surface storage.
  std::size_t cell{0};
  std::size_t face{0};
  std::size_t meshId{0};
  /// Among the written elements.
  std::size_t index{0};
};

/// Elements per work item of an evaluation: enough to amortise a call, few enough to balance.
constexpr std::size_t ChunkElements = 64;

template <typename Cfg>
class DerivedElementOutput final : public DerivedOutput {
  public:
  using RealT = Real<Cfg>;

  /// `faces`: the elements are faces of the free surface, which also offer their displacement.
  DerivedElementOutput(seissol::SeisSol& seissolInstance,
                       std::size_t writtenCount,
                       std::vector<WrittenElement> written,
                       const DerivedPoints& points,
                       const expr::Program& program,
                       bool faces);

  [[nodiscard]] const std::vector<std::string>& names() const override { return names_; }

  [[nodiscard]] bool accumulates() const override { return !derived_.program().state().empty(); }

  void write(double time) override {
    const std::lock_guard lock(mutex_);
    std::vector<std::size_t> pending;
    for (std::size_t i = 0; i < runs_.size(); ++i) {
      if (!accumulates() || !runs_[i].stepped) {
        pending.push_back(i);
      }
      runs_[i].stepped = false;
    }
    evaluate(pending, time);
  }

  void step(std::size_t layerId, double time) override {
    const std::lock_guard lock(mutex_);
    const auto found = runsOfLayer_.find(layerId);
    if (found == runsOfLayer_.end()) {
      return;
    }
    for (const auto run : found->second) {
      runs_[run].stepped = true;
    }
    evaluate(found->second, time);
  }

  void copy(std::size_t output,
            double* target,
            std::size_t index,
            std::size_t subcell) const override {
    const auto& values = simulations_[output / outputCount_].values;
    const std::size_t pointsPerSubcell = derived_.pointsPerSubcell();
    const std::size_t first = (output % outputCount_) * numPoints_ +
                              location_.at(index) * derived_.pointsPerElement() +
                              subcell * pointsPerSubcell;
    std::copy_n(values.data() + first, pointsPerSubcell, target);
  }

  private:
  /// The written elements of one layer and variant: a contiguous range of the point set.
  struct Run {
    std::size_t layer{0};
    std::size_t variant{0};
    std::size_t firstElement{0};
    std::size_t elementCount{0};
    /// Per simulation, the bases of the blocks of the sources in this layer.
    std::vector<std::vector<const void*>> blocks;
    /// The bases of the matrices of the variant.
    std::vector<const void*> matrices;
    /// Per simulation, the bases of the states in this layer, where the storage keeps them.
    std::vector<std::vector<void*>> states;
    /// The time of the previous evaluation, and the time and time step of the current one.
    double lastTime{0};
    double time{0};
    double timeStep{0};
    bool stepped{false};
  };

  struct Simulation {
    std::unique_ptr<DataTable> table;
    std::unique_ptr<expr::Binding> binding;
    /// One per thread: a kernel is not to be run from two threads at once.
    std::vector<std::unique_ptr<expr::Kernel>> kernels;
    /// Output j at point p: values[j * numPoints + p].
    std::vector<double> values;
  };

  /// Every quantity a program of this configuration can read; the time integral and the
  /// plastic strain are checked to be stored only once a program reads them.
  static std::vector<DerivedSource>
      sourcesOf(bool plasticity, bool faces, std::vector<SourceLocation>& locations);

  /// The reals from one cell (or face) to the next in the array of `storage`.
  static std::size_t cellStride(Storage storage);

  const void* blockBase(const SourceLocation& location, std::size_t color, std::size_t simulation);

  void computeGeometry(const std::vector<WrittenElement>& written, bool faces);

  /// The storage keeps the state, unless the program keeps state the configured layout does not
  /// know; the Binding keeps it then, without checkpoints.
  void resolveState(bool faces);

  /// The base of state `state` of simulation `simulation` in the layer `color`.
  void* stateBase(std::size_t state, std::size_t color, std::size_t simulation) const;

  void evaluate(const std::vector<std::size_t>& runs, double time);

  seissol::SeisSol& seissolInstance_;
  std::vector<SourceLocation> locations_;
  DerivedProgram derived_;
  std::vector<std::string> names_;
  std::size_t outputCount_{0};
  std::size_t numPoints_{0};
  /// Written element -> element of the point set; only read for the elements of this
  /// configuration.
  std::vector<std::size_t> location_;
  /// Element of the point set -> its cell in its layer of the cell storage, and its face in its
  /// layer of the surface storage.
  std::vector<std::uint32_t> cellIndex_;
  std::vector<std::uint32_t> faceIndex_;
  /// Per element of the point set: the inverse Jacobian of its cell, row-major, and its affine
  /// map (see DerivedProgram::bindGeometry()).
  std::vector<double> jacobians_;
  std::vector<double> transforms_;
  /// Where the storage keeps the state: of the faces or of the cells, in which layout, and the
  /// state of the layout each state of the program is.
  bool faces_{false};
  DerivedStateLayout stateLayout_;
  std::vector<std::size_t> stateOfLayout_;
  std::vector<Run> runs_;
  std::map<std::size_t, std::vector<std::size_t>> runsOfLayer_;
  std::vector<Simulation> simulations_;
  std::optional<std::size_t> timeInput_;
  std::optional<std::size_t> timeStepInput_;
  // the bases the time columns are bound with; every call moves them to its run
  double boundTime_{0};
  double boundTimeStep_{0};
  std::mutex mutex_;
};

template <typename Cfg>
std::vector<DerivedSource> DerivedElementOutput<Cfg>::sourcesOf(
    bool plasticity, bool faces, std::vector<SourceLocation>& locations) {
  using MaterialT = model::MaterialOf<Cfg>;
  // the stride between two quantities, i.e. the padded coefficients of one
  constexpr auto ModalStride =
      tensor::Q<Cfg>::Size / tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  constexpr auto NodalStride = tensor::QStressNodal<Cfg>::Size /
                               tensor::QStressNodal<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  constexpr auto FaceStride = tensor::faceDisplacement<Cfg>::Size /
                              tensor::faceDisplacement<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];

  std::vector<DerivedSource> sources;
  for (std::size_t q = 0; q < MaterialT::Quantities.size(); ++q) {
    sources.push_back(DerivedSource{MaterialT::Quantities[q], Representation::Modal});
    locations.push_back(SourceLocation{Storage::Dofs, q * ModalStride});
  }
  for (std::size_t q = 0; q < MaterialT::Quantities.size(); ++q) {
    sources.push_back(
        DerivedSource{IntegralPrefix + MaterialT::Quantities[q], Representation::Modal});
    locations.push_back(SourceLocation{Storage::Integrals, q * ModalStride});
  }
  if (plasticity) {
    for (std::size_t p = 0; p < model::PlasticityData<Cfg>::Quantities.size(); ++p) {
      sources.push_back(
          DerivedSource{model::PlasticityData<Cfg>::Quantities[p], Representation::Nodal});
      locations.push_back(SourceLocation{Storage::PlasticStrain, p * NodalStride});
    }
  }
  if (faces) {
    for (std::size_t component = 0; component < 3; ++component) {
      sources.push_back(
          DerivedSource{"u" + std::to_string(component + 1), Representation::FaceNodal});
      locations.push_back(SourceLocation{Storage::FaceDisplacement, component * FaceStride});
    }
  }
  return sources;
}

template <typename Cfg>
std::size_t DerivedElementOutput<Cfg>::cellStride(Storage storage) {
  switch (storage) {
  case Storage::Dofs:
  case Storage::Integrals:
    return tensor::Q<Cfg>::size();
  case Storage::PlasticStrain:
    return tensor::QStressNodal<Cfg>::size() + tensor::QEtaNodal<Cfg>::size();
  case Storage::FaceDisplacement:
    return tensor::faceDisplacement<Cfg>::size();
  }
  return 0;
}

template <typename Cfg>
void DerivedElementOutput<Cfg>::resolveState(bool faces) {
  faces_ = faces;
  const auto& states = derived_.program().state();
  if (states.empty()) {
    return;
  }
  auto layout = derivedStateLayout(
      seissolInstance_, faces ? DerivedOutputKind::Surface : DerivedOutputKind::WaveField);
  std::vector<std::size_t> stateOfLayout;
  for (const auto& state : states) {
    const auto found = std::find(layout.states.begin(), layout.states.end(), state.name);
    if (found == layout.states.end()) {
      logWarning() << "derived output: the state" << state.name
                   << "is not one of the configured program, so it is not written to "
                      "checkpoints.";
      return;
    }
    stateOfLayout.push_back(static_cast<std::size_t>(found - layout.states.begin()));
  }
  if (layout.simulations < Cfg::NumSimulations) {
    logError() << "derived output: the storage keeps the state of" << layout.simulations
               << "simulations, not of" << Cfg::NumSimulations << ".";
  }
  stateLayout_ = std::move(layout);
  stateOfLayout_ = std::move(stateOfLayout);
}

template <typename Cfg>
void* DerivedElementOutput<Cfg>::stateBase(std::size_t state,
                                           std::size_t color,
                                           std::size_t simulation) const {
  double* base = nullptr;
  if (faces_) {
    base = seissolInstance_.memoryManager()
               .surfaceStorage()
               .layer(color)
               .template var<SurfaceLTS::DerivedState>();
  } else {
    base = seissolInstance_.memoryManager()
               .ltsStorage()
               .layer(color)
               .template var<LTS::DerivedState>();
  }
  return base == nullptr ? nullptr
                         : base + stateLayout_.offset(stateOfLayout_[state], 0, simulation);
}

template <typename Cfg>
const void* DerivedElementOutput<Cfg>::blockBase(const SourceLocation& location,
                                                 std::size_t color,
                                                 std::size_t simulation) {
  auto& layer = seissolInstance_.memoryManager().ltsStorage().layer(color);
  const RealT* base = nullptr;
  switch (location.storage) {
  case Storage::Dofs:
    base = reinterpret_cast<const RealT*>(layer.template var<LTS::Dofs>(Cfg()));
    break;
  case Storage::Integrals:
    base = reinterpret_cast<const RealT*>(layer.template var<LTS::Integrals>(Cfg()));
    break;
  case Storage::PlasticStrain:
    base = reinterpret_cast<const RealT*>(layer.template var<LTS::PStrain>(Cfg()));
    break;
  case Storage::FaceDisplacement: {
    // the surface storage has a layer for every layer of the cell storage
    auto& surfaceLayer = seissolInstance_.freeSurfaceIntegrator().surfaceStorage->layer(color);
    base = reinterpret_cast<const RealT*>(
        surfaceLayer.template var<SurfaceLTS::DisplacementDofs>(Cfg()));
    break;
  }
  }
  // the simulation is the fastest dimension of the fused layout
  return base == nullptr ? nullptr : base + location.offset + simulation;
}

template <typename Cfg>
void DerivedElementOutput<Cfg>::computeGeometry(const std::vector<WrittenElement>& written,
                                                bool faces) {
  if (!derived_.readsJacobian() && !derived_.readsCoordinates()) {
    return;
  }
  const auto& meshReader = seissolInstance_.meshReader();
  const auto barycenter =
      seissol::geometry::CellTransform::VectorEigenT(Cell::ReferenceBarycenter.data());
  jacobians_.resize(written.size() * 9);
  transforms_.assign(written.size() * TransformSize, 0.0);
  for (std::size_t i = 0; i < written.size(); ++i) {
    const auto cell =
        seissol::geometry::AffineTransform::fromMeshCell(written[i].meshId, meshReader);
    // the transform is affine, so its inverse Jacobian is the one at the barycenter
    const auto inverse = cell.refToSpaceJacobianInverse(barycenter);
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t column = 0; column < 3; ++column) {
        jacobians_[i * 9 + row * 3 + column] = inverse(row, column);
      }
    }

    // the origin, then the images of the reference unit vectors less the origin; a face has two
    double* transform = transforms_.data() + i * TransformSize;
    if (faces) {
      using FaceVectorT = seissol::geometry::FaceTransform::FaceVectorT;
      const auto face = seissol::geometry::AffineFaceTransform::fromMeshCell(
          written[i].meshId, written[i].variant, meshReader);
      const auto origin = face.refToSpace(FaceVectorT(0.0, 0.0));
      for (std::size_t axis = 0; axis < 3; ++axis) {
        transform[axis] = origin(axis);
      }
      for (std::size_t direction = 0; direction < 2; ++direction) {
        FaceVectorT unit(0.0, 0.0);
        unit(direction) = 1.0;
        const auto image = face.refToSpace(unit);
        for (std::size_t axis = 0; axis < 3; ++axis) {
          transform[3 + direction * 3 + axis] = image(axis) - origin(axis);
        }
      }
    } else {
      const auto origin = cell.refToSpace(std::array<double, 3>{0, 0, 0});
      for (std::size_t axis = 0; axis < 3; ++axis) {
        transform[axis] = origin[axis];
      }
      for (std::size_t direction = 0; direction < 3; ++direction) {
        std::array<double, 3> unit{0, 0, 0};
        unit[direction] = 1;
        const auto image = cell.refToSpace(unit);
        for (std::size_t axis = 0; axis < 3; ++axis) {
          transform[3 + direction * 3 + axis] = image[axis] - origin[axis];
        }
      }
    }
  }
}

template <typename Cfg>
DerivedElementOutput<Cfg>::DerivedElementOutput(seissol::SeisSol& seissolInstance,
                                                std::size_t writtenCount,
                                                std::vector<WrittenElement> written,
                                                const DerivedPoints& points,
                                                const expr::Program& program,
                                                bool faces)
    : seissolInstance_(seissolInstance),
      derived_(program,
               sourcesOf(seissolInstance.parameters().model.plasticity, faces, locations_),
               points) {
  constexpr std::size_t Simulations = Cfg::NumSimulations;
  const std::size_t pointsPerElement = derived_.pointsPerElement();

  for (std::size_t sim = 0; sim < Simulations; ++sim) {
    for (const auto& output : derived_.program().outputs()) {
      names_.push_back(multisim::MultisimHelperWrapper<Cfg>::MultisimEnabled
                           ? output.name + "-" + std::to_string(sim + 1)
                           : output.name);
    }
  }
  outputCount_ = derived_.program().outputs().size();

  // The point set: the written elements of this configuration, layer by layer and variant by
  // variant, in memory order.
  std::sort(written.begin(), written.end(), [](const WrittenElement& a, const WrittenElement& b) {
    return std::tie(a.color, a.variant, a.cell, a.face) <
           std::tie(b.color, b.variant, b.cell, b.face);
  });
  numPoints_ = written.size() * pointsPerElement;
  location_.assign(writtenCount, std::numeric_limits<std::size_t>::max());
  cellIndex_.reserve(written.size());
  faceIndex_.reserve(written.size());
  for (std::size_t i = 0; i < written.size(); ++i) {
    location_.at(written[i].index) = i;
    cellIndex_.push_back(static_cast<std::uint32_t>(written[i].cell));
    faceIndex_.push_back(static_cast<std::uint32_t>(written[i].face));
    if (runs_.empty() || runs_.back().layer != written[i].color ||
        runs_.back().variant != written[i].variant) {
      runsOfLayer_[written[i].color].push_back(runs_.size());
      runs_.emplace_back();
      runs_.back().layer = written[i].color;
      runs_.back().variant = written[i].variant;
      runs_.back().firstElement = i;
    }
    ++runs_.back().elementCount;
  }

  // The bases of the blocks per layer and simulation, and of the matrices per variant.
  const auto& sources = derived_.usedSources();
  for (auto& run : runs_) {
    run.blocks.resize(Simulations);
    for (std::size_t sim = 0; sim < Simulations; ++sim) {
      for (std::size_t b = 0; b < sources.size(); ++b) {
        const void* base = blockBase(locations_[sources[b]], run.layer, sim);
        if (base == nullptr) {
          logError() << "derived output: the program reads"
                     << derived_.program().blocks()[b].name.c_str()
                     << ", which this run does not keep (the time integral of the solution is "
                        "only kept when an output reads it, the plastic strain only with "
                        "plasticity).";
        }
        run.blocks[sim].push_back(base);
      }
    }
    run.matrices = derived_.matrixBases(run.variant);
  }

  resolveState(faces);
  // the time of the previous evaluation: the start of the run, or the time it restarts from
  const double start = seissolInstance.simulator().getCurrentTime();
  for (auto& run : runs_) {
    run.lastTime = start;
    if (!stateOfLayout_.empty()) {
      run.states.resize(Simulations);
      for (std::size_t sim = 0; sim < Simulations; ++sim) {
        for (std::size_t state = 0; state < stateOfLayout_.size(); ++state) {
          run.states[sim].push_back(stateBase(state, run.layer, sim));
        }
      }
    }
  }

  computeGeometry(written, faces);

  const auto& inputs = derived_.program().inputs();
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    if (inputs[i].name == "t") {
      timeInput_ = i;
    }
    if (inputs[i].name == "dt") {
      timeStepInput_ = i;
    }
  }

  // One table, binding and set of kernels per simulation: the state of a point is its own in
  // every simulation.
  const std::size_t threads = std::max<std::size_t>(1, seissol::OpenMP::threadCount());
  auto& grids = reader::datafield::sharedGridStore();
  simulations_.resize(Simulations);
  for (std::size_t sim = 0; sim < Simulations; ++sim) {
    auto& simulation = simulations_[sim];
    simulation.values.assign(outputCount_ * numPoints_, 0.0);
    simulation.table = std::make_unique<DataTable>(numPoints_);
    DataTable& table = *simulation.table;

    derived_.bindMatrices(table);
    for (std::size_t b = 0; b < sources.size(); ++b) {
      const auto& location = locations_[sources[b]];
      reader::scripting::BlockInput block;
      block.base = runs_.front().blocks[sim][b];
      block.cellStride = cellStride(location.storage) * sizeof(RealT);
      block.modeStride = Simulations * sizeof(RealT);
      block.cellIndex =
          location.storage == Storage::FaceDisplacement ? faceIndex_.data() : cellIndex_.data();
      block.type = reader::scripting::DataTypeTraits<RealT>::Type;
      block.length = derived_.program().blocks()[b].length;
      table.bindBlock(derived_.program().blocks()[b].name, block);
    }
    derived_.bindGeometry(table, transforms_.data());
    // the state per element where the storage keeps it, [state][point][simulation] in one
    for (std::size_t state = 0; state < stateOfLayout_.size(); ++state) {
      reader::scripting::StateInput input;
      input.base = runs_.front().states[sim][state];
      input.pointsPerCell = pointsPerElement;
      input.cellStride = stateLayout_.size() * sizeof(double);
      input.pointStride = stateLayout_.simulations * sizeof(double);
      input.cellIndex = faces ? faceIndex_.data() : cellIndex_.data();
      input.type = DataType::F64;
      table.bindState(derived_.program().state()[state].name, input);
    }
    if (derived_.readsJacobian()) {
      for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
          table.bindCellView<double>(
              jacobianName(row, column), jacobians_.data(), pointsPerElement, 9, row * 3 + column);
        }
      }
    }
    if (timeInput_.has_value()) {
      table.bindViewConst<double>("t", Direction::In, &boundTime_, 0);
    }
    if (timeStepInput_.has_value()) {
      table.bindViewConst<double>("dt", Direction::In, &boundTimeStep_, 0);
    }
    for (std::size_t j = 0; j < outputCount_; ++j) {
      table.bindView<double>(derived_.program().outputs()[j].name,
                             Direction::Out,
                             simulation.values.data() + j * numPoints_);
    }

    simulation.binding =
        std::make_unique<expr::Binding>(expr::Binding::bind(derived_.program(), table));

    // The first kernel says which backend it got and how it went; the others, made for the same
    // program, follow it quietly.
    expr::BackendOptions options;
    options.preferred = expr::BackendKind::RtcCpu;
    for (std::size_t thread = 0; thread < threads; ++thread) {
      auto kernel = expr::makeKernel(derived_.program(), *simulation.binding, grids, options);
      kernel->precompute(table);
      options.preferred = kernel->kind();
      options.quiet = true;
      simulation.kernels.push_back(std::move(kernel));
    }
  }
}

template <typename Cfg>
void DerivedElementOutput<Cfg>::evaluate(const std::vector<std::size_t>& runs, double time) {
  struct WorkItem {
    std::size_t run;
    std::size_t simulation;
    std::size_t firstElement;
    std::size_t elementCount;
  };
  std::vector<WorkItem> items;
  for (const auto r : runs) {
    auto& run = runs_[r];
    run.time = time;
    run.timeStep = time - run.lastTime;
    run.lastTime = time;
    for (std::size_t sim = 0; sim < simulations_.size(); ++sim) {
      for (std::size_t element = 0; element < run.elementCount; element += ChunkElements) {
        items.push_back(WorkItem{r,
                                 sim,
                                 run.firstElement + element,
                                 std::min(ChunkElements, run.elementCount - element)});
      }
    }
  }

  const std::size_t pointsPerElement = derived_.pointsPerElement();
  const std::size_t inputCount = derived_.program().inputs().size();
  // at most one thread per kernel
  const int threads = static_cast<int>(simulations_.front().kernels.size());
#pragma omp parallel num_threads(threads)
  {
    std::vector<const void*> inputs(inputCount, nullptr);
    const std::size_t thread = seissol::OpenMP::threadId();
#pragma omp for schedule(dynamic)
    for (std::size_t i = 0; i < items.size(); ++i) {
      const auto& item = items[i];
      const auto& run = runs_[item.run];
      auto& simulation = simulations_[item.simulation];
      if (timeInput_.has_value()) {
        inputs[*timeInput_] = &run.time;
      }
      if (timeStepInput_.has_value()) {
        inputs[*timeStepInput_] = &run.timeStep;
      }
      expr::KernelArgs args;
      args.inputs = inputs.data();
      args.inputCount = inputs.size();
      args.matrices = run.matrices.data();
      args.matrixCount = run.matrices.size();
      // the blocks of the affine maps come after those of the sources, and stay where bound
      args.blocks = run.blocks[item.simulation].data();
      args.blockCount = run.blocks[item.simulation].size();
      if (!run.states.empty()) {
        args.states = run.states[item.simulation].data();
        args.stateCount = run.states[item.simulation].size();
      }
      args.first = item.firstElement * pointsPerElement;
      args.count = item.elementCount * pointsPerElement;
      simulation.kernels[thread]->run(args);
    }
  }
}

template <typename Cfg>
std::shared_ptr<DerivedOutput> makeDerivedOutput(seissol::SeisSol& seissolInstance,
                                                 std::size_t writtenCount,
                                                 std::vector<WrittenElement> written,
                                                 const DerivedPoints& points,
                                                 const expr::Program& program,
                                                 bool faces) {
  if (written.empty()) {
    return nullptr;
  }
  try {
    return std::make_shared<DerivedElementOutput<Cfg>>(
        seissolInstance, writtenCount, std::move(written), points, program, faces);
  } catch (const std::invalid_argument& error) {
    logError() << error.what();
  }
  return nullptr;
}

} // namespace

template <typename Cfg>
std::shared_ptr<DerivedOutput>
    makeDerivedVolumeOutput(seissol::SeisSol& seissolInstance,
                            const std::shared_ptr<const std::vector<std::size_t>>& cells,
                            const DerivedGeometry& geometry,
                            const expr::Program& program) {
  auto& storage = seissolInstance.memoryManager().ltsStorage();
  auto& backmap = seissolInstance.memoryManager().backmap();
  std::vector<WrittenElement> written;
  for (std::size_t index = 0; index < cells->size(); ++index) {
    const auto position = backmap.get(cells->at(index));
    if (storage.layer(position.color).getIdentifier().config == configIdOf<Cfg>()) {
      written.push_back(
          WrittenElement{position.color, 0, position.cell, 0, cells->at(index), index});
    }
  }
  return makeDerivedOutput<Cfg>(
      seissolInstance, cells->size(), std::move(written), VolumePoints(geometry), program, false);
}

template <typename Cfg>
std::shared_ptr<DerivedOutput> makeDerivedSurfaceOutput(seissol::SeisSol& seissolInstance,
                                                        const DerivedSurfaceGeometry& geometry,
                                                        const expr::Program& program) {
  auto& storage = seissolInstance.memoryManager().ltsStorage();
  auto& backmap = seissolInstance.memoryManager().backmap();
  const auto& integrator = seissolInstance.freeSurfaceIntegrator();
  auto& surfaceStorage = *integrator.surfaceStorage;
  const auto* meshIds = surfaceStorage.var<SurfaceLTS::MeshId>();
  const auto* sides = surfaceStorage.var<SurfaceLTS::Side>();

  // The free surface integrator numbers the faces through the surface layers in order; a face
  // lies in the layer of its cell.
  std::vector<std::pair<std::size_t, std::size_t>> layerOfFace;
  for (const auto& layer : surfaceStorage.leaves(LayerMask(Ghost))) {
    for (std::size_t face = 0; face < layer.size(); ++face) {
      layerOfFace.emplace_back(layer.id(), face);
    }
  }

  std::vector<WrittenElement> written;
  for (std::size_t index = 0; index < integrator.backmap.size(); ++index) {
    const auto face = integrator.backmap[index];
    const auto meshId = meshIds[face];
    const auto position = backmap.get(meshId);
    if (storage.layer(position.color).getIdentifier().config != configIdOf<Cfg>()) {
      continue;
    }
    const auto [color, local] = layerOfFace.at(face);
    if (color != position.color) {
      logError() << "derived output: a face of the free surface does not lie in the layer of "
                    "its cell.";
    }
    written.push_back(WrittenElement{position.color,
                                     static_cast<std::size_t>(sides[face]),
                                     position.cell,
                                     local,
                                     meshId,
                                     index});
  }
  return makeDerivedOutput<Cfg>(seissolInstance,
                                integrator.backmap.size(),
                                std::move(written),
                                SurfacePoints(geometry),
                                program,
                                true);
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template std::shared_ptr<DerivedOutput> makeDerivedVolumeOutput<Cfg>(                            \
      seissol::SeisSol&,                                                                           \
      const std::shared_ptr<const std::vector<std::size_t>>&,                                      \
      const DerivedGeometry&,                                                                      \
      const expr::Program&);                                                                       \
  template std::shared_ptr<DerivedOutput> makeDerivedSurfaceOutput<Cfg>(                           \
      seissol::SeisSol&, const DerivedSurfaceGeometry&, const expr::Program&);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::initializer
