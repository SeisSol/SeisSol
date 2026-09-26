// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_POINT_HDF5TABLE_H_
#define SEISSOL_SRC_IO_INSTANCE_POINT_HDF5TABLE_H_

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Instance/Point/Grouping.h"
#include "IO/Writer/Writer.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

namespace seissol::io::instance::point {

/**
 * @brief The sampled values of a set of points, as one HDF5 dataset per quantity set.
 *
 * A point records the quantities its material has, so not all of them record the same ones. The
 * points are therefore gathered into groups that share a quantity set (see Grouping), and every
 * group becomes a dataset of its own rather than one wide table with holes in it.
 *
 * A dataset is indexed by sample and by point, and one element holds a whole sample of one point
 * as a compound of its quantities. That is the same memory as a (sample, point, quantity) array
 * of numbers, with the quantity axis named instead of numbered: a reader gets the names and the
 * types out of the file, and a group whose quantities are not all of one type still fits.
 *
 * The sample axis grows with every write and the point axis is split across the ranks, so all the
 * samples of one point lie together -- which is how the files are read afterwards, one point at a
 * time.
 */
class Hdf5Table {
  public:
  /**
   * @brief Sets up the tables for the points this rank holds.
   *
   * @p pointQuantities is what each local point records, in the caller's own point order.
   * Collective on @p comm : which groups exist has to be agreed on, since a rank takes part in
   * declaring a dataset even when it holds no point of that group.
   *
   * @p sampleChunk is how far a storage chunk reaches along the sample axis. Zero lets the
   * writer decide. The chunking is fixed when the dataset is created, so a run whose writes
   * carry a varying number of samples wants this set rather than derived from whichever number
   * the first write happened to have.
   */
  Hdf5Table(std::string name,
            const std::vector<std::vector<TableQuantity>>& pointQuantities,
            MPI_Comm comm,
            std::size_t sampleChunk = 0);

  //! @brief How the points were gathered and renumbered.
  [[nodiscard]] const Grouping& grouping() const;

  //! @brief Bytes one sample of one point of @p group takes.
  [[nodiscard]] std::size_t sampleSize(std::size_t group) const;

  //! @brief Points of @p group this rank holds.
  [[nodiscard]] std::size_t localPointCount(std::size_t group) const;

  //! @brief Where local point @p point sits in the block of its group that this rank holds.
  [[nodiscard]] std::size_t localRow(std::size_t point) const;

  /**
   * @brief Adds a value per point that does not change, written once next to the tables.
   *
   * In the caller's point order, like the point map, and split across the ranks the same way, so
   * that a reader lines the two up row by row. @p shape is what the value of one point holds.
   */
  template <typename T>
  void addPointData(const std::string& name,
                    const std::vector<std::size_t>& shape,
                    const std::vector<T>& values) {
    PointData entry;
    entry.name = name;
    entry.shape = shape;
    entry.datatype = datatype::inferDatatype<T>();
    entry.bytes.resize(values.size() * sizeof(T));
    std::memcpy(entry.bytes.data(), values.data(), entry.bytes.size());
    pointData_.emplace_back(std::move(entry));
  }

  /**
   * @brief Storage for the samples of @p group that the next write is to carry.
   *
   * Laid out (sample, point) with the points in the order the grouping renumbered them into, and
   * one sample of one point taking sampleSize() bytes. The caller fills it and it keeps its
   * contents until the next call for that group, which is what the write needs: the plan is
   * handed on and read later.
   *
   * Every rank has to announce the same @p samples for a group, since that is how far the
   * dataset grows.
   */
  [[nodiscard]] char* prepare(std::size_t group, std::size_t samples);

  //! @brief Drops the samples of every group, after a write has been planned.
  void clear();

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter();

  private:
  //! A value per point that is written once, alongside the tables.
  struct PointData {
    std::string name;
    std::vector<std::size_t> shape;
    std::shared_ptr<datatype::Datatype> datatype;
    std::vector<char> bytes;
  };

  std::string name_;
  Grouping grouping_;
  std::size_t sampleChunk_;
  //! Per group, the samples handed over for the next write, laid out (sample, point).
  std::vector<std::vector<char>> storage_;
  //! Per group, how many samples that is.
  std::vector<std::size_t> samples_;
  //! Per group, how many points this rank holds.
  std::vector<std::size_t> localPoints_;
  //! Per local point, where it sits in the block of its group that this rank holds.
  std::vector<std::size_t> localRow_;
  //! Per group, whether its dataset still lacks the attributes that describe it. They go with
  //! the first write that carries samples, since that is the one that creates the dataset.
  std::vector<bool> undescribed_;
  //! Whether this run has written the table before.
  bool started_{false};
  //! The point map, kept alive for as long as the write that carries it.
  std::vector<std::uint64_t> index_;
  std::vector<PointData> pointData_;
};

} // namespace seissol::io::instance::point

#endif // SEISSOL_SRC_IO_INSTANCE_POINT_HDF5TABLE_H_
