// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Datafield/TableGrid.h"

#include "Expr/Core.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::reader::datafield {

namespace {

/// The fields of a line: separated by spaces or commas.
std::vector<std::string> fieldsOf(const std::string& line) {
  std::string spaced = line;
  std::replace(spaced.begin(), spaced.end(), ',', ' ');
  std::istringstream stream(spaced);
  std::vector<std::string> fields;
  std::string field;
  while (stream >> field) {
    fields.push_back(field);
  }
  return fields;
}

/// The number `spelling`: one of sderiv, nan or inf, with a sign. Throws std::invalid_argument.
double numberOf(const std::string& spelling) {
  const bool negative = !spelling.empty() && spelling[0] == '-';
  const std::string body =
      negative || (!spelling.empty() && spelling[0] == '+') ? spelling.substr(1) : spelling;
  double value = 0.0;
  if (body == "nan") {
    value = std::numeric_limits<double>::quiet_NaN();
  } else if (body == "inf") {
    value = std::numeric_limits<double>::infinity();
  } else {
    value = expr::numberOf(body);
  }
  return negative ? -value : value;
}

/// `value` as briefly as it reads back.
std::string spelled(double value) {
  std::array<char, 32> text{};
  const auto result = std::to_chars(text.data(), text.data() + text.size(), value);
  return {text.data(), result.ptr};
}

/// The names `names` as a list for a message.
std::string listed(const std::vector<std::string>& names) {
  std::string list;
  for (const auto& name : names) {
    list += (list.empty() ? "" : ", ") + name;
  }
  return list;
}

} // namespace

const std::vector<double>* TextTable::column(const std::string& name) const {
  const auto found = std::find(names.begin(), names.end(), name);
  return found == names.end() ? nullptr : &columns[static_cast<std::size_t>(found - names.begin())];
}

TextTable parseTextTable(const std::string& text, const std::string& path) {
  std::istringstream lines(text);
  TextTable table;
  std::string line;
  std::size_t number = 0;
  const auto fail = [&](const std::string& what) {
    throw std::runtime_error(path + ":" + std::to_string(number) + ": " + what);
  };
  while (std::getline(lines, line)) {
    ++number;
    const auto hash = line.find('#');
    const auto fields = fieldsOf(hash == std::string::npos ? line : line.substr(0, hash));
    if (fields.empty()) {
      continue;
    }
    if (table.names.empty()) {
      for (std::size_t c = 0; c < fields.size(); ++c) {
        if (std::find(fields.begin(), fields.begin() + static_cast<std::ptrdiff_t>(c), fields[c]) !=
            fields.begin() + static_cast<std::ptrdiff_t>(c)) {
          fail("the column `" + fields[c] + "` is named twice");
        }
      }
      table.names = fields;
      table.columns.resize(fields.size());
      continue;
    }
    if (fields.size() != table.names.size()) {
      fail(std::to_string(fields.size()) + " values for " + std::to_string(table.names.size()) +
           " columns");
    }
    for (std::size_t c = 0; c < fields.size(); ++c) {
      try {
        table.columns[c].push_back(numberOf(fields[c]));
      } catch (const std::exception&) {
        fail("`" + fields[c] + "` is no number");
      }
    }
    table.lines.push_back(number);
  }
  if (table.names.empty()) {
    throw std::runtime_error(path + ": no line names the columns of the table");
  }
  return table;
}

TextTable readTextTable(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream.good()) {
    throw std::runtime_error("Could not open " + path + " for reading.");
  }
  return parseTextTable(
      std::string(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()), path);
}

GridData tableGridData(const TextTable& table,
                       const std::vector<std::string>& coordinates,
                       const std::vector<std::string>& components,
                       const std::string& path) {
  const auto columnOf = [&](const std::string& name) {
    const auto* column = table.column(name);
    if (column == nullptr) {
      throw std::runtime_error("the table " + path + " has no column `" + name + "`; it has " +
                               listed(table.names));
    }
    return column;
  };
  if (coordinates.empty() || coordinates.size() > MaxGridDimensions) {
    throw std::runtime_error("a table has 1 to " + std::to_string(MaxGridDimensions) +
                             " coordinates, not " + std::to_string(coordinates.size()));
  }
  std::vector<const std::vector<double>*> axes;
  for (const auto& name : coordinates) {
    axes.push_back(columnOf(name));
  }
  std::vector<const std::vector<double>*> values;
  for (const auto& name : components) {
    values.push_back(columnOf(name));
  }

  // the knots of every axis: the distinct values of its column
  GridData data;
  std::size_t nodes = 1;
  for (std::size_t d = 0; d < axes.size(); ++d) {
    for (std::size_t row = 0; row < table.rows(); ++row) {
      if (!std::isfinite((*axes[d])[row])) {
        throw std::runtime_error(path + ":" + std::to_string(table.lines[row]) +
                                 ": the coordinate `" + coordinates[d] + "` is " +
                                 spelled((*axes[d])[row]) + ", not finite");
      }
    }
    std::vector<double> knots = *axes[d];
    std::sort(knots.begin(), knots.end());
    knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
    if (knots.empty()) {
      throw std::runtime_error("the table " + path + " has no rows");
    }
    if (knots.size() > MaxKnotsOfFile) {
      throw std::runtime_error("the coordinate `" + coordinates[d] + "` of the table " + path +
                               " has " + std::to_string(knots.size()) + " values, more than " +
                               std::to_string(MaxKnotsOfFile));
    }
    nodes *= knots.size();
    data.knots.push_back(std::move(knots));
  }

  // every row a node, every node once
  const auto where = [&](const std::vector<std::size_t>& index) {
    std::string text = "(";
    for (std::size_t d = 0; d < index.size(); ++d) {
      text += (d == 0 ? "" : ", ") + coordinates[d] + " = " + spelled(data.knots[d][index[d]]);
    }
    return text + ")";
  };
  constexpr std::size_t None = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> rowOf(nodes, None);
  std::vector<std::size_t> index(axes.size());
  for (std::size_t row = 0; row < table.rows(); ++row) {
    std::size_t node = 0;
    for (std::size_t d = axes.size(); d-- > 0;) {
      const auto& knots = data.knots[d];
      index[d] = static_cast<std::size_t>(
          std::lower_bound(knots.begin(), knots.end(), (*axes[d])[row]) - knots.begin());
      node = node * knots.size() + index[d];
    }
    if (rowOf[node] != None) {
      throw std::runtime_error(path + ":" + std::to_string(table.lines[row]) + ": the node " +
                               where(index) + " is given twice, before in line " +
                               std::to_string(table.lines[rowOf[node]]));
    }
    rowOf[node] = row;
  }
  data.values.resize(nodes * components.size());
  for (std::size_t node = 0; node < nodes; ++node) {
    if (rowOf[node] == None) {
      std::size_t rest = node;
      for (std::size_t d = 0; d < axes.size(); ++d) {
        index[d] = rest % data.knots[d].size();
        rest /= data.knots[d].size();
      }
      throw std::runtime_error("the table " + path + " gives no node " + where(index) +
                               ": its rows are the nodes of a grid, each once");
    }
    for (std::size_t c = 0; c < components.size(); ++c) {
      data.values[node * components.size() + c] = (*values[c])[rowOf[node]];
    }
  }
  return data;
}

} // namespace seissol::reader::datafield
