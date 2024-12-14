/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017 - 2020, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Setup of SeisSol's cell local matrices.
 **/

#ifndef INITIALIZER_PARAMETERDB_H_
#define INITIALIZER_PARAMETERDB_H_

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
QueryGenerator* getBestQueryGenerator(bool anelasticity,
                                      bool plasticity,
                                      bool anisotropy,
                                      bool poroelasticity,
                                      bool useCellHomogenizedMaterial,
                                      const CellToVertexArray& cellToVertex);

class QueryGenerator {
  public:
  virtual ~QueryGenerator() = default;
  [[nodiscard]] virtual easi::Query generate() const = 0;
};

class ElementBarycentreGenerator : public QueryGenerator {
  public:
  explicit ElementBarycentreGenerator(const CellToVertexArray& cellToVertex)
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

class FaultBarycentreGenerator : public QueryGenerator {
  public:
  FaultBarycentreGenerator(const seissol::geometry::MeshReader& meshReader, unsigned numberOfPoints)
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
  virtual void evaluateModel(const std::string& fileName, const QueryGenerator* queryGen) = 0;
  static easi::Component* loadModel(const std::string& fileName);
};

template <class T>
class MaterialParameterDB : ParameterDB {
  public:
  T computeAveragedMaterial(unsigned elementIdx,
                            const std::array<double, NumQuadpoints>& quadratureWeights,
                            const std::vector<T>& materialsFromQuery);
  void evaluateModel(const std::string& fileName, const QueryGenerator* queryGen) override;
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
  void evaluateModel(const std::string& fileName, const QueryGenerator* queryGen) override;
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

#endif
