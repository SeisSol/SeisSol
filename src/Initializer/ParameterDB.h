// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_

#include <memory>
#include <set>
#include <string>
#include <unordered_map>

#include "Geometry/MeshReader.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"

#include "easi/Query.h"
#include "easi/ResultAdapter.h"
#include "generated_code/init.h"

#include "Equations/Datastructures.h"

#ifdef USE_HDF
#include "PUML/PUML.h"
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
  static CellToVertexArray fromPUML(const PUML::TETPUML& mesh);
#endif
  static CellToVertexArray
      fromVectors(const std::vector<std::array<std::array<double, 3>, 4>>& vertices,
                  const std::vector<int>& groups);
};

easi::Component* loadEasiModel(const std::string& fileName);
std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool plasticity,
                                                      bool useCellHomogenizedMaterial,
                                                      const CellToVertexArray& cellToVertex);

class QueryGenerator {
  public:
  virtual ~QueryGenerator() = default;
  [[nodiscard]] virtual easi::Query generate() const = 0;
};

class ElementBarycenterGenerator : public QueryGenerator {
  public:
  explicit ElementBarycenterGenerator(const CellToVertexArray& cellToVertex)
      : m_cellToVertex(cellToVertex) {}
  [[nodiscard]] easi::Query generate() const override;

  private:
  CellToVertexArray m_cellToVertex;
};

class ElementAverageGenerator : public QueryGenerator {
  public:
  explicit ElementAverageGenerator(const CellToVertexArray& cellToVertex);
  [[nodiscard]] easi::Query generate() const override;
  [[nodiscard]] const std::array<double, NumQuadpoints>& getQuadratureWeights() const {
    return m_quadratureWeights;
  };

  private:
  CellToVertexArray m_cellToVertex;
  std::array<double, NumQuadpoints> m_quadratureWeights{};
  std::array<std::array<double, 3>, NumQuadpoints> m_quadraturePoints{};
};

class FaultBarycenterGenerator : public QueryGenerator {
  public:
  FaultBarycenterGenerator(const seissol::geometry::MeshReader& meshReader, unsigned numberOfPoints)
      : m_meshReader(meshReader), m_numberOfPoints(numberOfPoints) {}
  [[nodiscard]] easi::Query generate() const override;

  private:
  const seissol::geometry::MeshReader& m_meshReader;
  unsigned m_numberOfPoints;
};

class FaultGPGenerator : public QueryGenerator {
  public:
  FaultGPGenerator(const seissol::geometry::MeshReader& meshReader,
                   const std::vector<unsigned>& faceIDs)
      : m_meshReader(meshReader), m_faceIDs(faceIDs) {}
  [[nodiscard]] easi::Query generate() const override;

  private:
  const seissol::geometry::MeshReader& m_meshReader;
  const std::vector<unsigned>& m_faceIDs;
};

class ParameterDB {
  public:
  virtual ~ParameterDB() = default;
  virtual void evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) = 0;
  static easi::Component* loadModel(const std::string& fileName);
};

template <class T>
class MaterialParameterDB : ParameterDB {
  public:
  T computeAveragedMaterial(unsigned elementIdx,
                            const std::array<double, NumQuadpoints>& quadratureWeights,
                            const std::vector<T>& materialsFromQuery);
  void evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) override;
  void setMaterialVector(std::vector<T>* materials) { m_materials = materials; }
  void addBindingPoints(easi::ArrayOfStructsAdapter<T>& adapter) {};

  private:
  std::vector<T>* m_materials{};
};

class FaultParameterDB : ParameterDB {
  public:
  ~FaultParameterDB() override = default;
  void addParameter(const std::string& parameter, real* memory, unsigned stride = 1) {
    m_parameters[parameter] = std::make_pair(memory, stride);
  }
  void evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) override;
  static std::set<std::string> faultProvides(const std::string& fileName);

  private:
  std::unordered_map<std::string, std::pair<real*, unsigned>> m_parameters;
};

class EasiBoundary {
  public:
  explicit EasiBoundary(const std::string& fileName);

  EasiBoundary() : model(nullptr) {};
  EasiBoundary(const EasiBoundary&) = delete;
  EasiBoundary& operator=(const EasiBoundary&) = delete;
  EasiBoundary(EasiBoundary&& other) noexcept;
  EasiBoundary& operator=(EasiBoundary&& other) noexcept;

  ~EasiBoundary();

  void query(const real* nodes, real* mapTermsData, real* constantTermsData) const;

  private:
  easi::Component* model;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_
