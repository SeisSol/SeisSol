// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_

#include "Common/Real.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Instance/Mesh/Xdmf.h"

#include <variant>

namespace seissol::io::instance::geometry {

enum class WriterFormat : int32_t { Xdmf, Vtk };

enum class WriterBackend : int32_t { Binary, Hdf5 };

enum class WriterGroup : int32_t {
  FullSnapshot,
  // IncrementalSnapshot, ?
  Monolith
};

struct WriterConfig {
  uint32_t order{0};
  WriterFormat format{WriterFormat::Vtk};
  WriterBackend backend{WriterBackend::Hdf5};
  WriterGroup time{WriterGroup::FullSnapshot};
  int32_t compress{0};
};

class GeometryWriter {
  private:
  static std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter>
      getUnderlyingWriter(const std::string& name,
                          std::size_t localElementCount,
                          Shape shape,
                          const WriterConfig& config) {
    if (config.format == WriterFormat::Xdmf) {
      return mesh::XdmfWriter(name,
                              localElementCount,
                              shape,
                              config.order,
                              config.backend == WriterBackend::Binary,
                              config.compress);
    } else {
      return mesh::VtkHdfWriter(name,
                                localElementCount,
                                shape,
                                config.order,
                                config.time == WriterGroup::Monolith,
                                config.compress);
    }
  }

  public:
  GeometryWriter(const std::string& name,
                 std::size_t localElementCount,
                 Shape shape,
                 const WriterConfig& config,
                 std::size_t subdivide = 1)
      : config(config), subdivide(subdivide),
        underlying(getUnderlyingWriter(name, subdivide * localElementCount, shape, config)) {}

  template <typename F>
  void addPointProjector(F projector) {
    const auto subdivide = this->subdivide;
    std::visit(
        [&](auto& writer) {
          writer.addPointProjector([projector, subdivide](auto* data, std::size_t index) {
            std::invoke(projector, data, index / subdivide, index % subdivide);
          });
        },
        underlying);
  }

  template <typename T, typename F>
  void addGeometryOutput(const std::string& name,
                         const std::vector<std::size_t>& dimensions,
                         bool isConst,
                         F writerFunction) {

    const auto subdivide = this->subdivide;
    if (config.order == 0) {
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
          underlying);
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
          underlying);
    }
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   bool isConst,
                   F writerFunction) {

    const auto subdivide = this->subdivide;
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
        underlying);
  }

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() {
    return std::visit([&](auto& writer) { return writer.makeWriter(); }, underlying);
  }

  void addHook(const std::function<void(std::size_t, double)>& hook) {
    std::visit([&](auto& writer) { return writer.addHook(hook); }, underlying);
  }

  protected:
  WriterConfig config;
  std::size_t subdivide;
  std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter> underlying;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
