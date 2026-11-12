// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_INLINEGRID_H_
#define SEISSOL_SRC_READER_DATAFIELD_INLINEGRID_H_

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace seissol::reader::datafield {

/// A grid whose values and knots its program gives (GridKind::Inline): the layers of a layered
/// model, say. Every axis is knots; the values lie in memory as GridData has them. Only the
/// written-out lookup samples it.
class InlineGrid final : public Grid {
  public:
  /// Throws std::invalid_argument for data that does not form a grid: an axis without knots,
  /// knots that do not strictly increase or are not finite, or values of another number than the
  /// nodes times `components`.
  InlineGrid(std::shared_ptr<const GridData> data, std::size_t components);

  [[nodiscard]] std::size_t dimensions() const override { return data_->knots.size(); }
  [[nodiscard]] std::size_t components() const override { return components_; }
  void geometry(GridGeometry& out) const override;
  void sample(const double* x, double* out) const override;
  void sampleBatch(const double* const* x, std::size_t count, double* const* out) const override;
  void sampleBatch(const double* const* x, std::size_t count, float* const* out) const override;
  [[nodiscard]] const ArrayView* arrayView() const override { return &view_; }

  /// What is wrong with `data` as the data of a grid of `components` components, or nothing.
  [[nodiscard]] static std::string problemOf(const GridData& data, std::size_t components);

  private:
  std::shared_ptr<const GridData> data_;
  std::size_t components_;
  std::vector<std::vector<double>> inverseSpacings_;
  ArrayView view_;
};

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_INLINEGRID_H_
