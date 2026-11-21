// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FAULTPARAMETERDB_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FAULTPARAMETERDB_H_

#include "Core/Data/DataReader.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Geometry/QueryGenerator.h"

#include <cstddef>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::initializer {

/// The quadrature points of the given fault faces, in the quadrature rule of the configuration
/// `Cfg`.
template <typename Cfg>
class FaultGPGenerator : public QueryGenerator {
  public:
  FaultGPGenerator(const seissol::geometry::MeshReader& meshReader,
                   const std::vector<std::size_t>& faceIDs)
      : meshReader_(meshReader), faceIDs_(faceIDs) {}
  [[nodiscard]] data::DataTable generate() const override;

  private:
  const seissol::geometry::MeshReader& meshReader_;
  const std::vector<std::size_t>& faceIDs_;
};

/**
 * The parameters of the fault faces of the simulation `simulation` of `numSimulations` fused ones,
 * written into arrays of `T`.
 */
template <typename T>
class FaultParameterDB {
  public:
  FaultParameterDB(std::size_t simulation, std::size_t numSimulations)
      : simid_(simulation), numSimulations_(numSimulations) {}
  void addParameter(const std::string& parameter, T* memory, unsigned stride = 1) {
    parameters_[parameter] = std::make_pair(memory, stride);
  }
  /**
   * Writes the queried faces to the positions `faces` of the arrays, with `pointsPerFace` points
   * per face and fused simulation; by default, the queried faces are all faces of the arrays, in
   * order.
   */
  void setFaces(std::vector<std::size_t> faces, std::size_t pointsPerFace) {
    faces_ = std::move(faces);
    pointsPerFace_ = pointsPerFace;
  }
  /// Evaluates the fault model `model` at the points of `queryGen`.
  void evaluate(data::DataReader& model, const QueryGenerator& queryGen);

  private:
  std::size_t simid_;
  std::size_t numSimulations_;
  std::unordered_map<std::string, std::pair<T*, unsigned>> parameters_;
  std::optional<std::vector<std::size_t>> faces_;
  std::size_t pointsPerFace_{1};
};

/// The parameters that the fault parameter file `fileName` provides, which `models` opens; none if
/// there is no file.
std::set<std::string> faultProvides(const data::ModelLoader& models, const std::string& fileName);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FAULTPARAMETERDB_H_
