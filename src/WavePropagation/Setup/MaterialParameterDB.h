// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SETUP_MATERIALPARAMETERDB_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SETUP_MATERIALPARAMETERDB_H_

#include "Core/Data/DataReader.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Geometry/QueryGenerator.h"

#include <array>
#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace seissol::initializer {

/**
 * The query generator for the material `MaterialT`: if useCellHomogenizedMaterial and the material
 * can be averaged, the average over each cell with the quadrature for the given convergence order;
 * otherwise, the barycenter of each cell.
 */
template <typename MaterialT>
std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool useCellHomogenizedMaterial,
                                                      const CellToVertexArray& cellToVertex,
                                                      std::size_t convergenceOrder);

class ElementBarycenterGenerator : public QueryGenerator {
  public:
  explicit ElementBarycenterGenerator(CellToVertexArray cellToVertex)
      : cellToVertex_(std::move(cellToVertex)) {}
  [[nodiscard]] data::DataTable generate() const override;

  private:
  CellToVertexArray cellToVertex_;
};

/// Queries the points of the quadrature for the given convergence order in each cell.
class ElementAverageGenerator : public QueryGenerator {
  public:
  ElementAverageGenerator(CellToVertexArray cellToVertex, std::size_t convergenceOrder);
  [[nodiscard]] data::DataTable generate() const override;
  [[nodiscard]] const std::vector<double>& getQuadratureWeights() const {
    return quadratureWeights_;
  };

  private:
  CellToVertexArray cellToVertex_;
  std::vector<double> quadratureWeights_;
  std::vector<std::array<double, Cell::Dim>> quadraturePoints_;
};

template <class T>
class MaterialParameterDB {
  public:
  //! @brief Evaluates the material model @p model at the points of @p queryGen .
  void evaluate(data::DataReader& model, const QueryGenerator& queryGen);
  void setMaterialVector(std::vector<T>* materials) { materials_ = materials; }

  private:
  std::vector<T>* materials_{};
};

/**
 * @brief The parameters that a model which supplies @p supplied lacks for materials of the type
 * @p T , as MaterialParameterDB evaluates it: none if it supplies those of T, or those of an
 * elastic or acoustic material that T can be made of; else the ones of T it does not supply.
 */
template <class T>
std::vector<std::string> missingMaterialParameters(const std::set<std::string>& supplied);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_WAVEPROPAGATION_SETUP_MATERIALPARAMETERDB_H_
