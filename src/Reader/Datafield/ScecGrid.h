// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_SCECGRID_H_
#define SEISSOL_SRC_READER_DATAFIELD_SCECGRID_H_

// A grid of a SCEC file of fault parameters (GridKind::Scec): the layout of the input files of the
// SCEC benchmarks TPV16 and TPV17, as easi's !SCECFile reads it (easi's SCECFile.cpp). Numbers
// separated by space: a header of three lines -- two integers (a seed); the numbers of cells along
// strike and along dip, and the lengths of the fault along them; the hypocenter, two integers and
// four numbers -- then a line per node: its indices along strike and along dip, its coordinates,
// and ten values. The components of the grid are those ten values, in the order of the file: the
// normal stress, the shear stresses along strike and along dip, their ratios to the normal stress,
// the static and the dynamic coefficient of friction, the critical slip distance, the cohesion,
// and the time of forced rupture.
//
// As in easi, a node is where its indices put it -- its coordinates are not read -- on a uniform
// grid from 0 to the lengths, with a node more than cells along each axis, whose spacing is the
// length over the cells; axis 0 is along strike. Unlike easi, every node is to be given exactly
// once, and the lines after the last node are not read.

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::reader::datafield {

/// What a SCEC file holds: the numbers of nodes and the lengths along strike and dip, and the
/// values of the nodes (empty where only the header is read), the ten of a node together, the
/// nodes along strike first.
struct ScecData {
  static constexpr std::size_t Components = 10;
  std::array<unsigned, 2> nodes{};
  std::array<double, 2> lengths{};
  std::vector<double> values;

  /// The spacing along `axis`: the length over the cells, as easi's regular grid has it (1 for an
  /// axis of one node, whose spacing no sample reads).
  [[nodiscard]] double spacing(std::size_t axis) const {
    return nodes[axis] > 1 ? (lengths[axis] - 0.0) / static_cast<double>(nodes[axis] - 1) : 1.0;
  }
};

/// The names of the ten components of a SCEC file, in their order.
[[nodiscard]] const std::vector<std::string>& scecComponentNames();

/// The header of the SCEC file `path`. Throws std::runtime_error, with the line, where the file
/// cannot be read or has no header of the layout.
[[nodiscard]] ScecData readScecHeader(const std::string& path);

/// The SCEC file `path`, its values included. Throws as readScecHeader(), and where a node is
/// missing, given twice, or outside the grid.
[[nodiscard]] ScecData readScecFile(const std::string& path);

class ScecGrid final : public Grid {
  public:
  /// The grid of `data`, read from `path`, which sample() and sampleBatch() interpolate with
  /// `interpolation`.
  ScecGrid(ScecData data, std::string path, Interpolation interpolation, BoundaryMode boundary);

  [[nodiscard]] std::size_t dimensions() const override { return 2; }
  [[nodiscard]] std::size_t components() const override { return ScecData::Components; }
  void geometry(GridGeometry& out) const override;
  void sample(const double* x, double* out) const override;
  void sampleBatch(const double* const* x, std::size_t count, double* const* out) const override;
  void sampleBatch(const double* const* x, std::size_t count, float* const* out) const override;
  [[nodiscard]] const ArrayView* arrayView() const override { return &view_; }

  private:
  void reportClamped(std::size_t clamped, std::size_t count) const;

  ScecData data_;
  std::string path_;
  Interpolation interpolation_;
  BoundaryMode boundary_;
  ArrayView view_;
  mutable SampleScratch scratch_;
};

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_SCECGRID_H_
