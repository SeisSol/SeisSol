// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Datafield/ScecGrid.h"

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
#include "utils/logger.h"

#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace seissol::reader::datafield {

namespace {

/// The numbers of a text, one after another, as `>>` reads them: separated by space, line breaks
/// included; a failure names the line.
class Numbers {
  public:
  Numbers(std::string path, std::string text) : path_(std::move(path)), text_(std::move(text)) {}

  /// The next number, an integer from 0 to `limit`.
  std::uint64_t index(const char* what, std::uint64_t limit) {
    const std::string token = next(what);
    std::uint64_t value = 0;
    const auto* end = token.data() + token.size();
    const auto [last, error] = std::from_chars(token.data(), end, value);
    if (error != std::errc() || last != end || value > limit) {
      fail(std::string(what) + " is an integer from 0 to " + std::to_string(limit) + ", not `" +
           token + "`");
    }
    return value;
  }

  /// The next number, an integer.
  std::int64_t integer(const char* what) {
    std::string token = next(what);
    std::int64_t value = 0;
    const char* begin = token.data() + (token.size() > 1 && token[0] == '+' ? 1 : 0);
    const auto* end = token.data() + token.size();
    const auto [last, error] = std::from_chars(begin, end, value);
    if (error != std::errc() || last != end) {
      fail(std::string(what) + " is an integer, not `" + token + "`");
    }
    return value;
  }

  /// The next number, rounded to the nearest double as strtod does.
  double real(const char* what) {
    const std::string token = next(what);
    double value = 0.0;
    const char* begin = token.data() + (token.size() > 1 && token[0] == '+' ? 1 : 0);
    const auto* end = token.data() + token.size();
    const auto [last, error] = std::from_chars(begin, end, value, std::chars_format::general);
    if (error != std::errc() || last != end) {
      fail(std::string(what) + " is a number, not `" + token + "`");
    }
    return value;
  }

  [[noreturn]] void fail(const std::string& what) const {
    throw std::runtime_error(path_ + ":" + std::to_string(line_) + ": " + what);
  }

  private:
  std::string next(const char* what) {
    while (at_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[at_])) != 0) {
      if (text_[at_] == '\n') {
        ++line_;
      }
      ++at_;
    }
    if (at_ == text_.size()) {
      fail(std::string("the file ends where it has ") + what);
    }
    const std::size_t start = at_;
    while (at_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[at_])) == 0) {
      ++at_;
    }
    return text_.substr(start, at_ - start);
  }

  std::string path_;
  std::string text_;
  std::size_t at_{0};
  int line_{1};
};

std::string textOf(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream.good()) {
    throw std::runtime_error("Could not open " + path + " for reading.");
  }
  return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

/// The header: the seed, the cells and lengths, the hypocenter.
ScecData header(Numbers& numbers) {
  ScecData data;
  static_cast<void>(numbers.integer("the first integer of the seed"));
  static_cast<void>(numbers.integer("the second integer of the seed"));
  // a node more than cells, and as many as an index can address
  constexpr std::uint64_t MostCells = (std::uint64_t{1} << 31U) - 2;
  for (std::size_t d = 0; d < 2; ++d) {
    data.nodes[d] = static_cast<unsigned>(
        numbers.index(d == 0 ? "the cells along strike" : "the cells along dip", MostCells) + 1);
  }
  data.lengths[0] = numbers.real("the length along strike");
  data.lengths[1] = numbers.real("the length along dip");
  // easi's message, where it samples the grid
  for (std::size_t d = 0; d < 2; ++d) {
    if (!(data.spacing(d) > 0.0)) {
      numbers.fail("Degenerate grid in dimension " + std::to_string(d) + " (spacing " +
                   std::to_string(data.spacing(d)) + ", " + std::to_string(data.nodes[d]) +
                   " points).");
    }
  }
  static_cast<void>(numbers.integer("the first index of the hypocenter"));
  static_cast<void>(numbers.integer("the second index of the hypocenter"));
  for (const char* what : {"the first coordinate of the hypocenter",
                           "the second coordinate of the hypocenter",
                           "the third number of the hypocenter",
                           "the fourth number of the hypocenter"}) {
    static_cast<void>(numbers.real(what));
  }
  return data;
}

} // namespace

const std::vector<std::string>& scecComponentNames() {
  static const std::vector<std::string> Names = {"s_normal",
                                                 "s_strike",
                                                 "s_dip",
                                                 "s_strike_over_s_normal",
                                                 "s_dip_over_s_normal",
                                                 "mu_s",
                                                 "mu_d",
                                                 "d_c",
                                                 "cohesion",
                                                 "forced_rupture_time"};
  return Names;
}

ScecData readScecHeader(const std::string& path) {
  Numbers numbers(path, textOf(path));
  return header(numbers);
}

ScecData readScecFile(const std::string& path) {
  Numbers numbers(path, textOf(path));
  ScecData data = header(numbers);
  const std::size_t nodes = std::size_t{data.nodes[0]} * data.nodes[1];
  data.values.assign(nodes * ScecData::Components, 0.0);
  std::vector<bool> given(nodes, false);
  for (std::size_t line = 0; line < nodes; ++line) {
    const auto i = numbers.index("the index along strike", data.nodes[0] - 1);
    const auto j = numbers.index("the index along dip", data.nodes[1] - 1);
    const std::size_t node = i + j * data.nodes[0];
    if (given[node]) {
      numbers.fail("the node (" + std::to_string(i) + ", " + std::to_string(j) +
                   ") is given twice");
    }
    given[node] = true;
    static_cast<void>(numbers.real("the coordinate along strike"));
    static_cast<void>(numbers.real("the coordinate along dip"));
    for (std::size_t c = 0; c < ScecData::Components; ++c) {
      data.values[node * ScecData::Components + c] = numbers.real("a value of the node");
    }
  }
  return data;
}

ScecGrid::ScecGrid(ScecData data,
                   std::vector<std::size_t> columns,
                   std::string path,
                   Interpolation interpolation,
                   BoundaryMode boundary)
    : data_(std::move(data)), columns_(std::move(columns)), path_(std::move(path)),
      interpolation_(interpolation), boundary_(boundary) {
  const std::size_t nodes = std::size_t{data_.nodes[0]} * data_.nodes[1];
  if (data_.values.size() != nodes * ScecData::Components) {
    throw std::invalid_argument("datafield: " + std::to_string(data_.values.size()) +
                                " values for the " + std::to_string(nodes) +
                                " nodes of the SCEC file " + path_);
  }
  if (columns_.empty()) {
    throw std::invalid_argument("datafield: no component of the SCEC file " + path_);
  }
  for (const auto column : columns_) {
    if (column >= ScecData::Components) {
      throw std::invalid_argument("datafield: a SCEC file has " +
                                  std::to_string(ScecData::Components) + " values per node, not " +
                                  std::to_string(column + 1));
    }
  }
  // the components the grid has, of every node
  values_.resize(nodes * columns_.size());
  for (std::size_t node = 0; node < nodes; ++node) {
    for (std::size_t c = 0; c < columns_.size(); ++c) {
      values_[node * columns_.size() + c] = data_.values[node * ScecData::Components + columns_[c]];
    }
  }
  data_.values = std::vector<double>();
  view_.data = values_.data();
  view_.type = ElementType::Double;
  view_.dimensions = 2;
  view_.components = static_cast<unsigned>(columns_.size());
  view_.componentStride = 1;
  std::size_t stride = columns_.size();
  for (std::size_t d = 0; d < 2; ++d) {
    view_.min[d] = 0.0;
    view_.deltaInv[d] = 1.0 / data_.spacing(d);
    view_.num[d] = data_.nodes[d];
    view_.stride[d] = stride;
    stride *= data_.nodes[d];
  }
}

void ScecGrid::geometry(GridGeometry& out) const {
  out = GridGeometry{};
  out.dimensions = 2;
  for (std::size_t d = 0; d < 2; ++d) {
    out.min[d] = 0.0;
    out.delta[d] = data_.spacing(d);
    out.num[d] = data_.nodes[d];
  }
}

void ScecGrid::sample(const double* x, double* out) const {
  reportClamped(datafield::sample(view_, interpolation_, x, out), 1);
}

void ScecGrid::sampleBatch(const double* const* x, std::size_t count, double* const* out) const {
  reportClamped(datafield::sampleBatch(view_, interpolation_, x, count, out, scratch_), count);
}

void ScecGrid::sampleBatch(const double* const* x, std::size_t count, float* const* out) const {
  reportClamped(datafield::sampleBatch(view_, interpolation_, x, count, out, scratch_), count);
}

void ScecGrid::reportClamped(std::size_t clamped, std::size_t count) const {
  if (clamped != 0 && boundary_ == BoundaryMode::Fail) {
    logError() << "datafield: the SCEC file" << path_ << "was sampled outside its fault for"
               << clamped << "of" << count << "points.";
  }
}

} // namespace seissol::reader::datafield
