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

#include <Common/Templating.h>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>

#include "Geometry/MeshReader.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"

#include "GeneratedCode/init.h"
#include "easi/Query.h"
#include "easi/ResultAdapter.h"

#include "Equations/Datastructures.h"

#ifdef USE_HDF
#include "PUML/PUML.h"
#endif

#include <Eigen/Dense>

namespace easi {
class Component;
} // namespace easi

namespace seissol::initializer {
constexpr std::size_t AveragingOrder = 4;
constexpr auto NumQuadpoints = AveragingOrder * AveragingOrder * AveragingOrder;

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
  static CellToVertexArray join(std::vector<CellToVertexArray> arrays);

  [[nodiscard]] CellToVertexArray filter(const std::vector<bool>& keep) const;
};

easi::Component* loadEasiModel(const std::string& fileName);

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

class FaultGPGenerator : public QueryGenerator {
  public:
  FaultGPGenerator(const seissol::geometry::MeshReader& meshReader,
                   const std::vector<std::size_t>& faceIDs,
                   std::size_t configId)
      : m_meshReader(meshReader), m_faceIDs(faceIDs), configId(configId) {}
  [[nodiscard]] easi::Query generate() const override;

  private:
  const seissol::geometry::MeshReader& m_meshReader;
  const std::vector<std::size_t>& m_faceIDs;
  std::size_t configId;
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

template <typename T>
class FaultParameterDB : ParameterDB {
  public:
  explicit FaultParameterDB(std::size_t simulation, std::size_t simCount)
      : simid(simulation), simCount(simCount) {}
  ~FaultParameterDB() override = default;
  void addParameter(const std::string& parameter, T* memory, unsigned stride = 1) {
    m_parameters[parameter] = std::make_pair(memory, stride);
  }
  void evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) override;
  static std::set<std::string> faultProvides(const std::string& fileName);

  private:
  std::size_t simid;
  std::size_t simCount;
  std::unordered_map<std::string, std::pair<T*, unsigned>> m_parameters;
};

template <typename Cfg>
class EasiBoundary {
  public:
  using real = Real<Cfg>;
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

using EasiBoundaryT =
    TransformVariadicT<std::optional, TransformVariadicT<EasiBoundary, ConfigVariant>>;

template <typename MaterialT>
std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool plasticity,
                                                      bool useCellHomogenizedMaterial,
                                                      const CellToVertexArray& cellToVertex) {
  std::shared_ptr<QueryGenerator> queryGen;
  if (!useCellHomogenizedMaterial) {
    queryGen = std::make_shared<ElementBarycenterGenerator>(cellToVertex);
  } else {
    if (MaterialT::Type != model::MaterialType::Viscoelastic &&
        MaterialT::Type != model::MaterialType::Elastic) {
      logWarning() << "Material Averaging is not implemented for " << MaterialT::Text
                   << " materials. Falling back to "
                      "material properties sampled from the element barycenters instead.";
      queryGen = std::make_shared<ElementBarycenterGenerator>(cellToVertex);
    } else if (plasticity) {
      logWarning()
          << "Material Averaging is not implemented for plastic materials. Falling back to "
             "material properties sampled from the element barycenters instead.";
      queryGen = std::make_shared<ElementBarycenterGenerator>(cellToVertex);
    } else {
      queryGen = std::make_shared<ElementAverageGenerator>(cellToVertex);
    }
  }
  return queryGen;
}

using MaterialVariant =
    RemoveDuplicateVariadicT<TransformVariadicT<model::MaterialTT, ConfigVariant>>;

template <typename T>
std::vector<T> queryDB(const std::shared_ptr<seissol::initializer::QueryGenerator>& queryGen,
                       const std::string& fileName,
                       size_t size) {
  std::vector<T> vectorDB(size);
  seissol::initializer::MaterialParameterDB<T> parameterDB;
  parameterDB.setMaterialVector(&vectorDB);
  parameterDB.evaluateModel(fileName, *queryGen);
  return vectorDB;
}

std::vector<MaterialVariant>
    queryMaterials(const parameters::ModelParameters& params,
                   const seissol::initializer::CellToVertexArray& ctvArray);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_
