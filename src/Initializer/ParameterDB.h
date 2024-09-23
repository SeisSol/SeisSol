// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 */

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

namespace seissol {
namespace initializer {
constexpr auto NumQuadpoints = ConvergenceOrder * ConvergenceOrder * ConvergenceOrder;

class QueryGenerator;
class ElementBarycentreGenerator;
class ElementAverageGenerator;
class ElementBarycentreGeneratorPUML;
class FaultBarycentreGenerator;
class FaultGPGenerator;
class ParameterDB;
template <class T>
class MaterialParameterDB;
class FaultParameterDB;
class EasiBoundary;

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
QueryGenerator* getBestQueryGenerator(bool anelasticity,
                                      bool plasticity,
                                      bool anisotropy,
                                      bool poroelasticity,
                                      bool useCellHomogenizedMaterial,
                                      const CellToVertexArray& cellToVertex);
} // namespace initializer
} // namespace seissol

class seissol::initializer::QueryGenerator {
  public:
  virtual ~QueryGenerator() = default;
  virtual easi::Query generate() const = 0;
};

class seissol::initializer::ElementBarycentreGenerator
    : public seissol::initializer::QueryGenerator {
  public:
  explicit ElementBarycentreGenerator(const CellToVertexArray& cellToVertex)
      : m_cellToVertex(cellToVertex) {}
  easi::Query generate() const override;

  private:
  CellToVertexArray m_cellToVertex;
};

class seissol::initializer::ElementAverageGenerator : public seissol::initializer::QueryGenerator {
  public:
  explicit ElementAverageGenerator(const CellToVertexArray& cellToVertex);
  easi::Query generate() const override;
  const std::array<double, NumQuadpoints>& getQuadratureWeights() const {
    return m_quadratureWeights;
  };

  private:
  CellToVertexArray m_cellToVertex;
  std::array<double, NumQuadpoints> m_quadratureWeights;
  std::array<std::array<double, 3>, NumQuadpoints> m_quadraturePoints;
};

class seissol::initializer::FaultBarycentreGenerator : public seissol::initializer::QueryGenerator {
  public:
  FaultBarycentreGenerator(const seissol::geometry::MeshReader& meshReader, unsigned numberOfPoints)
      : m_meshReader(meshReader), m_numberOfPoints(numberOfPoints) {}
  easi::Query generate() const override;

  private:
  const seissol::geometry::MeshReader& m_meshReader;
  unsigned m_numberOfPoints;
};

class seissol::initializer::FaultGPGenerator : public seissol::initializer::QueryGenerator {
  public:
  FaultGPGenerator(const seissol::geometry::MeshReader& meshReader,
                   const std::vector<unsigned>& faceIDs)
      : m_meshReader(meshReader), m_faceIDs(faceIDs) {}
  easi::Query generate() const override;

  private:
  const seissol::geometry::MeshReader& m_meshReader;
  const std::vector<unsigned>& m_faceIDs;
};

class seissol::initializer::ParameterDB {
  public:
  virtual void evaluateModel(const std::string& fileName, const QueryGenerator* const queryGen) = 0;
  static easi::Component* loadModel(const std::string& fileName);
};

template <class T>
class seissol::initializer::MaterialParameterDB : seissol::initializer::ParameterDB {
  public:
  T computeAveragedMaterial(unsigned elementIdx,
                            const std::array<double, NumQuadpoints>& quadratureWeights,
                            const std::vector<T>& materialsFromQuery);
  void evaluateModel(const std::string& fileName, const QueryGenerator* const queryGen) override;
  void setMaterialVector(std::vector<T>* materials) { m_materials = materials; }
  void addBindingPoints(easi::ArrayOfStructsAdapter<T>& adapter) {};

  private:
  std::vector<T>* m_materials;
};

class seissol::initializer::FaultParameterDB : seissol::initializer::ParameterDB {
  public:
  void addParameter(const std::string& parameter, real* memory, unsigned stride = 1) {
    m_parameters[parameter] = std::make_pair(memory, stride);
  }
  void evaluateModel(const std::string& fileName, const QueryGenerator* const queryGen) override;
  static std::set<std::string> faultProvides(const std::string& fileName);

  private:
  std::unordered_map<std::string, std::pair<real*, unsigned>> m_parameters;
};

class seissol::initializer::EasiBoundary {
  public:
  explicit EasiBoundary(const std::string& fileName);

  EasiBoundary() : model(nullptr) {};
  EasiBoundary(const EasiBoundary&) = delete;
  EasiBoundary& operator=(const EasiBoundary&) = delete;
  EasiBoundary(EasiBoundary&& other);
  EasiBoundary& operator=(EasiBoundary&& other);

  ~EasiBoundary();

  void query(const real* nodes, real* mapTermsData, real* constantTermsData) const;

  private:
  easi::Component* model;
};

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERDB_H_
