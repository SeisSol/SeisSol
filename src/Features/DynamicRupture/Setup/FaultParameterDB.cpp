// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "FaultParameterDB.h"

#include "Config.h"
#include "Core/Common/Constants.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Data/DataReader.h"
#include "Core/Data/DataTable.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Geometry/CellTransform.h"
#include "Core/Geometry/FaceTransform.h"
#include "Core/Geometry/MeshDefinition.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Geometry/QueryGenerator.h"
#include "Features/DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace seissol::initializer {

template <typename Cfg>
data::DataTable FaultGPGenerator<Cfg>::generate() const {
  constexpr size_t NumPoints = dr::misc::NumPaddedPointsSingleSim<Cfg>;

  // element, side and the face transform of each fault face managed by this generator (we have
  // one generator per LTS layer), set up once rather than per point and column
  struct FaultFace {
    std::size_t element{0};
    std::int8_t side{0};
    seissol::geometry::AffineFaceTransform transform;
  };
  const std::vector<Fault>& fault = meshReader_.getFault();
  const auto cellToVertex = CellToVertexArray::fromMeshReader(meshReader_);
  auto faces = std::make_shared<std::vector<FaultFace>>();
  faces->reserve(faceIDs_.size());
  for (const auto& faultId : faceIDs_) {
    const Fault& f = fault.at(faultId);
    std::size_t element = 0;
    std::int8_t side = 0;
    auto sideOrientation = seissol::geometry::FaceOrientation::Local;
    if (f.element.hasValue()) {
      element = f.element.value();
      side = f.side;
    } else {
      assert(f.neighborElement.hasValue());

      element = f.neighborElement.value();
      side = f.neighborSide;
      // the canonical vertex numbering pins the face orientation index to zero
      sideOrientation = seissol::geometry::FaceOrientation::Rotate0;
    }
    faces->push_back(
        FaultFace{.element = element,
                  .side = side,
                  .transform = seissol::geometry::AffineFaceTransform(
                      seissol::geometry::AffineTransform(cellToVertex.elementCoordinates(element)),
                      seissol::geometry::ReferenceFaceMap(side, sideOrientation))});
  }

  data::DataTable table(NumPoints * faceIDs_.size());
  bindPointSet(
      table,
      [faces](std::size_t first, std::size_t count, double* xyz) {
        const auto pointsView = init::quadpoints<Cfg>::view::create(init::quadpoints<Cfg>::Values);
        for (std::size_t i = 0; i < count; ++i) {
          const auto q = first + i;
          const auto n = q % NumPoints;
          auto localPoints = seissol::geometry::FaceTransform::FaceVectorT(
              seissol::multisim::multisimTranspose<Cfg>(pointsView, n, 0),
              seissol::multisim::multisimTranspose<Cfg>(pointsView, n, 1));
          // padded points are in the middle of the tetrahedron
          if (n >= dr::misc::NumBoundaryGaussPoints<Cfg>) {
            localPoints =
                seissol::geometry::FaceTransform::FaceVectorT(Face::ReferenceBarycenter.data());
          }
          const auto point = (*faces)[q / NumPoints].transform.refToSpace(localPoints);
          for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
            xyz[i * Cell::Dim + dim] = point(dim);
          }
        }
      },
      [faces, &elements = meshReader_.getElements()](
          std::size_t first, std::size_t count, std::int32_t* out) {
        for (std::size_t i = 0; i < count; ++i) {
          const auto& face = (*faces)[(first + i) / NumPoints];
          out[i] = elements[face.element].faultTags[face.side];
        }
      });
  return table;
}

#define SEISSOL_INSTANTIATE(Cfg) template class FaultGPGenerator<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

template <typename T>
void FaultParameterDB<T>::evaluate(data::DataReader& model, const QueryGenerator& queryGen) {
  auto table = queryGen.generate();

  // the point q of the query lies on the face faces[q / pointsPerFace] of the arrays, if the faces
  // are given; i.e. the arrays may hold the points of faces which the query does not contain, in
  // between the ones it does
  std::vector<std::uint32_t> faceIndex;
  if (faces_.has_value()) {
    faceIndex.assign(faces_->begin(), faces_->end());
  }
  for (const auto& [name, binding] : parameters_) {
    const auto [memory, stride] = binding;
    const auto pointStride = static_cast<std::size_t>(stride) * numSimulations_;
    if (faces_.has_value()) {
      table.bindCellPointView(name,
                              data::Direction::Out,
                              memory,
                              pointsPerFace_,
                              pointsPerFace_ * pointStride,
                              pointStride,
                              simid_,
                              faceIndex.data());
    } else {
      table.bindView(name, data::Direction::Out, memory, pointStride, simid_);
    }
  }

  data::evaluate(model, table, "fault material");
}

template class FaultParameterDB<float>;
template class FaultParameterDB<double>;

std::set<std::string> faultProvides(const data::ModelLoader& models, const std::string& fileName) {
  if (fileName.empty()) {
    return {};
  }

  const auto model = models(fileName, {"x", "y", "z"});
  return data::suppliedParameters(*model);
}

} // namespace seissol::initializer
