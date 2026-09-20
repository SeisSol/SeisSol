// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_

#include "Common/Real.h"
#include "Deduplicate.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Instance/Mesh/Xdmf.h"

#include <algorithm>
#include <optional>
#include <string>
#include <utils/env.h>
#include <variant>
#include <vector>

namespace seissol::io::instance::geometry {

/**
 * @brief Whether coinciding points of an order 0 output are merged into one.
 *
 * On by default. Merging costs a pass over the local points at setup and saves about a factor of
 * twenty in the point array of a tetrahedral mesh; turning it off is for the case where a reader
 * wants the points of a cell to belong to that cell alone.
 *
 * SEISSOL_IO_VERTEXFILTER is the name to use; SEISSOL_VERTEXFILTER is what the previous writer
 * read and is accepted so that existing job scripts keep working.
 */
inline bool vertexFilterEnabled() {
  auto value = utils::Env("SEISSOL_IO_").getOptional<bool>("VERTEXFILTER");
  if (!value.has_value()) {
    value = utils::Env("SEISSOL_").getOptional<bool>("VERTEXFILTER");
  }
  return value.value_or(true);
}

enum class WriterFormat : int32_t { Xdmf, Vtk };

enum class WriterBackend : int32_t { Binary, Hdf5 };

enum class WriterGroup : int32_t {
  //! Every output step is a file of its own, and repeats the data that does not change.
  FullSnapshot,
  //! The data that does not change is written once, and the snapshots link to it.
  IncrementalSnapshot,
  //! One file holding all output steps, as a VTKHDF time series.
  Monolith
};

struct WriterConfig {
  uint32_t order{0};
  WriterFormat format{WriterFormat::Vtk};
  WriterBackend backend{WriterBackend::Hdf5};
  WriterGroup time{WriterGroup::FullSnapshot};
  uint32_t compress{0};
};

class GeometryWriter {
  private:
  template <typename F>
  static std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter>
      getUnderlyingWriter(const std::string& name,
                          std::size_t localElementCount,
                          Shape shape,
                          const WriterConfig& config,
                          std::size_t subdivide,
                          F projector) {
    if (config.time == WriterGroup::Monolith && config.format == WriterFormat::Xdmf) {
      logError() << "A monolithic time series output is only available for the VTKHDF format.";
    }
    if (config.time == WriterGroup::IncrementalSnapshot && config.format == WriterFormat::Xdmf) {
      // Xdmf keeps its payload in one file per output anyway and references the unchanging parts
      // from every time step, so there is nothing for a separate const file to save here.
      logError() << "An incremental snapshot output is only available for the VTKHDF format.";
    }
    // the projector is handed one cell at a time, subcells included
    const auto cellProjector = [projector, subdivide](double* data, std::size_t index) {
      std::invoke(projector, data, index / subdivide, index % subdivide);
    };

    std::optional<mesh::VertexMap> vertexMap;
    std::vector<double> uniquePoints;
    if (config.order == 0 && vertexFilterEnabled()) {
      // At degree 0 the points of a cell are its corners, and a corner belongs to every cell
      // around it -- about twenty of them in a tetrahedral mesh. Writing each of them once makes
      // the point array that much smaller. From degree 1 on the points are Lagrange nodes, which
      // neighbouring cells deliberately do not share, since the solution is discontinuous there.
      const auto pointsPerCell = numPoints(1, shape);
      std::vector<double> coordinates(localElementCount * pointsPerCell * 3);
      for (std::size_t cell = 0; cell < localElementCount; ++cell) {
        cellProjector(coordinates.data() + cell * pointsPerCell * 3, cell);
      }

      auto merged = deduplicatePoints(coordinates);
      // read the count before the points are moved out from under it
      const auto uniqueCount = merged.pointCount();
      uniquePoints = std::move(merged.points);
      vertexMap = mesh::VertexMap{uniqueCount, std::move(merged.indices)};
    }

    const auto pointProjector = [&](auto& writer) {
      if (vertexMap.has_value()) {
        writer.addPointProjector([points = uniquePoints](double* target, std::size_t index) {
          std::copy_n(points.data() + index * 3, 3, target);
        });
      } else {
        writer.addPointProjector(cellProjector);
      }
    };

    if (config.format == WriterFormat::Xdmf) {
      auto writer = mesh::XdmfWriter(name,
                                     localElementCount,
                                     shape,
                                     config.order,
                                     config.backend == WriterBackend::Binary,
                                     config.compress,
                                     vertexMap);
      pointProjector(writer);
      return writer;
    }

    auto writer = mesh::VtkHdfWriter(name,
                                     localElementCount,
                                     shape,
                                     config.order,
                                     config.time == WriterGroup::Monolith,
                                     config.compress,
                                     config.time == WriterGroup::IncrementalSnapshot,
                                     vertexMap);
    pointProjector(writer);
    return writer;
  }

  public:
  /**
   * @brief Creates the writer.
   *
   * @p projector is handed a cell and a subcell and fills the coordinates of that subcell's
   * points. It is taken here rather than added afterwards because the number of points a file
   * holds is only known once they have been looked at: at degree 0 the coinciding ones are merged.
   */
  template <typename F>
  GeometryWriter(const std::string& name,
                 std::size_t localElementCount,
                 Shape shape,
                 const WriterConfig& config,
                 std::size_t subdivide,
                 F projector)
      : config_(config), subdivide_(subdivide),
        underlying_(getUnderlyingWriter(
            name, subdivide * localElementCount, shape, config, subdivide, projector)) {}

  template <typename T, typename F>
  void addGeometryOutput(const std::string& name,
                         const std::vector<std::size_t>& dimensions,
                         bool isConst,
                         F writerFunction) {

    const auto subdivide = this->subdivide_;
    if (config_.order == 0) {
      // cell output
      std::visit(
          [&](auto& writer) {
            writer.template addCellData<T>(
                name,
                dimensions,
                isConst,
                [writerFunction, subdivide](auto* data, std::size_t index) {
                  std::invoke(writerFunction, data, index / subdivide, index % subdivide);
                });
          },
          underlying_);
    } else {
      // point output
      std::visit(
          [&](auto& writer) {
            writer.template addPointData<T>(
                name,
                dimensions,
                isConst,
                [writerFunction, subdivide](auto* data, std::size_t index) {
                  std::invoke(writerFunction, data, index / subdivide, index % subdivide);
                });
          },
          underlying_);
    }
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   bool isConst,
                   F writerFunction) {

    const auto subdivide = this->subdivide_;
    std::visit(
        [&](auto& writer) {
          writer.template addCellData<T>(
              name,
              dimensions,
              isConst,
              [writerFunction, subdivide](auto* data, std::size_t index) {
                std::invoke(writerFunction, data, index / subdivide, index % subdivide);
              });
        },
        underlying_);
  }

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() {
    return std::visit([&](auto& writer) { return writer.makeWriter(); }, underlying_);
  }

  void addHook(const std::function<void(std::size_t, double)>& hook) {
    std::visit([&](auto& writer) { return writer.addHook(hook); }, underlying_);
  }

  protected:
  WriterConfig config_;
  std::size_t subdivide_;
  std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter> underlying_;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
