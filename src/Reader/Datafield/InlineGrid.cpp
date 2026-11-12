// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Datafield/InlineGrid.h"

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace seissol::reader::datafield {

std::string InlineGrid::problemOf(const GridData& data, std::size_t components) {
  if (data.knots.empty() || data.knots.size() > MaxGridDimensions) {
    return "a grid has 1 to " + std::to_string(MaxGridDimensions) + " axes, not " +
           std::to_string(data.knots.size());
  }
  if (components == 0) {
    return "a grid has components";
  }
  std::size_t nodes = 1;
  for (std::size_t d = 0; d < data.knots.size(); ++d) {
    const auto& knots = data.knots[d];
    if (knots.empty()) {
      return "axis " + std::to_string(d) + " has no knots";
    }
    for (std::size_t i = 0; i < knots.size(); ++i) {
      if (!std::isfinite(knots[i])) {
        return "knot " + std::to_string(i) + " of axis " + std::to_string(d) + " is not finite";
      }
      if (i > 0 && !(knots[i - 1] < knots[i])) {
        return "the knots of axis " + std::to_string(d) + " do not increase strictly, at knot " +
               std::to_string(i);
      }
    }
    nodes *= knots.size();
  }
  if (data.values.size() != nodes * components) {
    return std::to_string(data.values.size()) + " values for " + std::to_string(nodes) +
           " nodes of " + std::to_string(components) + " components";
  }
  return {};
}

InlineGrid::InlineGrid(std::shared_ptr<const GridData> data, std::size_t components)
    : data_(std::move(data)), components_(components) {
  if (data_ == nullptr) {
    throw std::invalid_argument("datafield: an inline grid without data");
  }
  if (const auto problem = problemOf(*data_, components_); !problem.empty()) {
    throw std::invalid_argument("datafield: " + problem);
  }
  view_.data = data_->values.data();
  view_.type = ElementType::Double;
  view_.dimensions = static_cast<unsigned>(data_->knots.size());
  view_.components = static_cast<unsigned>(components_);
  view_.componentStride = 1;
  std::size_t stride = components_;
  inverseSpacings_.resize(data_->knots.size());
  for (std::size_t d = 0; d < data_->knots.size(); ++d) {
    const auto& knots = data_->knots[d];
    for (std::size_t i = 0; i + 1 < knots.size(); ++i) {
      inverseSpacings_[d].push_back(1.0 / (knots[i + 1] - knots[i]));
    }
    // a lookup reads the spacing of its cell, the first of a single knot as well
    if (inverseSpacings_[d].empty()) {
      inverseSpacings_[d].push_back(0.0);
    }
    view_.min[d] = knots.front();
    view_.num[d] = static_cast<unsigned>(knots.size());
    view_.stride[d] = stride;
    view_.knots[d] = knots.data();
    view_.inverseSpacings[d] = inverseSpacings_[d].data();
    stride *= knots.size();
  }
}

void InlineGrid::geometry(GridGeometry& out) const {
  out = GridGeometry{};
  out.dimensions = view_.dimensions;
  for (std::size_t d = 0; d < data_->knots.size(); ++d) {
    const auto& knots = data_->knots[d];
    out.min[d] = knots.front();
    // knots have no one spacing: the mean, which spans the axis as the knots do
    out.delta[d] = knots.size() > 1
                       ? (knots.back() - knots.front()) / static_cast<double>(knots.size() - 1)
                       : 0.0;
    out.num[d] = static_cast<unsigned>(knots.size());
  }
}

void InlineGrid::sample(const double* /*x*/, double* /*out*/) const {
  throw std::invalid_argument(
      "datafield: an inline grid is sampled through the written-out lookup only");
}

void InlineGrid::sampleBatch(const double* const* /*x*/,
                             std::size_t /*count*/,
                             double* const* /*out*/) const {
  throw std::invalid_argument(
      "datafield: an inline grid is sampled through the written-out lookup only");
}

void InlineGrid::sampleBatch(const double* const* /*x*/,
                             std::size_t /*count*/,
                             float* const* /*out*/) const {
  throw std::invalid_argument(
      "datafield: an inline grid is sampled through the written-out lookup only");
}

} // namespace seissol::reader::datafield
