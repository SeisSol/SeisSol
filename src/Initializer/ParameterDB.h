// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_

#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Geometry/MeshReader.h"
#include "Geometry/PUMLReader.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "easi/Query.h"
#include "easi/ResultAdapter.h"

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef USE_HDF
#include <PUML/PUML.h>
#endif

#include <Eigen/Dense>

namespace easi {
class Component;
} // namespace easi

namespace seissol::initializer {
constexpr auto NumQuadpoints = ConvergenceOrder * ConvergenceOrder * ConvergenceOrder;

class QueryGenerator;

// temporary struct until we have something like a lazy vector/iterator "map" (as in on-demand,
// element-wise function application)
struct CellToVertexArray {
  using CellToVertexFunction = std::function<std::array<Eigen::Vector3d, 4>(size_t)>;
  using CellToGroupFunction = std::function<int(size_t)>;

  CellToVertexArray(size_t size,
                    const CellToVertexFunction& elementCoordinates,
                    const CellToGroupFunction& elementGroups);

  size_t size;
  CellToVertexFunction elementCoordinates;
  CellToGroupFunction elementGroups;

  static CellToVertexArray fromMeshReader(const seissol::geometry::MeshReader& meshReader);
#ifdef USE_HDF
  /// The cells of a PUML mesh before PUMLReader::getMesh(), with their vertices in the given
  /// (canonical) order rather than in the order of the mesh file. The array refers to
  /// vertexOrders, which therefore has to outlive it.
  static CellToVertexArray
      fromPUML(const seissol::geometry::PumlMesh& mesh,
               const std::vector<seissol::geometry::VertexOrder>& vertexOrders);
  static CellToVertexArray
      fromPUML(const seissol::geometry::PumlMesh& mesh,
               std::vector<seissol::geometry::VertexOrder>&& vertexOrders) = delete;
#endif
  static CellToVertexArray
      fromVectors(const std::vector<std::array<std::array<double, 3>, 4>>& vertices,
                  const std::vector<int>& groups);
  static CellToVertexArray join(std::vector<CellToVertexArray> arrays);
};

std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool useCellHomogenizedMaterial,
                                                      const CellToVertexArray& cellToVertex);

class QueryGenerator {
  public:
  virtual ~QueryGenerator() = default;
  [[nodiscard]] virtual easi::Query generate() const = 0;
  [[nodiscard]] virtual std::size_t outputPerCell() const { return 1; }
};

class ElementBarycenterGenerator : public QueryGenerator {
  public:
  explicit ElementBarycenterGenerator(const CellToVertexArray& cellToVertex)
      : cellToVertex_(cellToVertex) {}
  [[nodiscard]] easi::Query generate() const override;

  private:
  CellToVertexArray cellToVertex_;
};

class ElementAverageGenerator : public QueryGenerator {
  public:
  explicit ElementAverageGenerator(const CellToVertexArray& cellToVertex);
  [[nodiscard]] easi::Query generate() const override;
  [[nodiscard]] const std::array<double, NumQuadpoints>& getQuadratureWeights() const {
    return quadratureWeights_;
  };

  private:
  CellToVertexArray cellToVertex_;
  std::array<double, NumQuadpoints> quadratureWeights_{};
  std::array<std::array<double, 3>, NumQuadpoints> quadraturePoints_{};
};

class PlasticityPointGenerator : public QueryGenerator {
  public:
  explicit PlasticityPointGenerator(const CellToVertexArray& cellToVertex, bool pointwise = true)
      : cellToVertex_(cellToVertex), pointwise_(pointwise) {}
  [[nodiscard]] easi::Query generate() const override;
  [[nodiscard]] std::size_t outputPerCell() const override;

  private:
  CellToVertexArray cellToVertex_;
  bool pointwise_{true};
};

class FaultBarycenterGenerator : public QueryGenerator {
  public:
  FaultBarycenterGenerator(const seissol::geometry::MeshReader& meshReader,
                           std::size_t numberOfPoints)
      : meshReader_(meshReader), numberOfPoints_(numberOfPoints) {}
  [[nodiscard]] easi::Query generate() const override;

  private:
  const seissol::geometry::MeshReader& meshReader_;
  std::size_t numberOfPoints_;
};

/// The quadrature points of the given fault faces, in the quadrature rule of the configuration
/// `Cfg`.
template <typename Cfg>
class FaultGPGenerator : public QueryGenerator {
  public:
  FaultGPGenerator(const seissol::geometry::MeshReader& meshReader,
                   const std::vector<std::size_t>& faceIDs)
      : meshReader_(meshReader), faceIDs_(faceIDs) {}
  [[nodiscard]] easi::Query generate() const override;

  private:
  const seissol::geometry::MeshReader& meshReader_;
  const std::vector<std::size_t>& faceIDs_;
};

class ParameterDB {
  public:
  virtual ~ParameterDB() = default;
  virtual void evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) = 0;
  static easi::Component* loadModel(const std::string& fileName);
};

template <class T>
class MaterialParameterDB : public ParameterDB {
  public:
  void evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) override;
  void setMaterialVector(std::vector<T>* materials) { materials_ = materials; }

  private:
  std::vector<T>* materials_{};
};

/**
 * The parameters of the fault faces of the simulation `simulation` of `numSimulations` fused ones,
 * written into arrays of `T`.
 */
template <typename T>
class FaultParameterDB : public ParameterDB {
  public:
  FaultParameterDB(std::size_t simulation, std::size_t numSimulations)
      : simid_(simulation), numSimulations_(numSimulations) {}
  ~FaultParameterDB() override = default;
  void addParameter(const std::string& parameter, T* memory, unsigned stride = 1) {
    parameters_[parameter] = std::make_pair(memory, stride);
  }
  void evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) override;

  private:
  std::size_t simid_;
  std::size_t numSimulations_;
  std::unordered_map<std::string, std::pair<T*, unsigned>> parameters_;
};

/// The parameters a fault parameter file provides.
std::set<std::string> faultProvides(const std::string& fileName);

/**
 * The frame the affine boundary condition is stated in. Global is the default; face-aligned
 * lets a condition be stated in terms of the face normal, which a condition that mirrors or
 * fixes a direction needs on a boundary that is not axis-aligned.
 */
enum class BoundaryFrame { Global, FaceAligned };

class DirichletCondition {
  public:
  explicit DirichletCondition(const std::string& fileName);

  DirichletCondition() : model_(nullptr) {};
  DirichletCondition(const DirichletCondition&) = delete;
  DirichletCondition& operator=(const DirichletCondition&) = delete;
  DirichletCondition(DirichletCondition&& other) noexcept;
  DirichletCondition& operator=(DirichletCondition&& other) noexcept;

  ~DirichletCondition();

  [[nodiscard]] BoundaryFrame
      query(const double* barycenter, real* mapTermsData, real* constantTermsData) const;

  private:
  easi::Component* model_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_
