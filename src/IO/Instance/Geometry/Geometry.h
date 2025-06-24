// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_

#include <IO/Instance/Mesh/VtkHdf.h>
#include <IO/Instance/Mesh/Xdmf.h>
#include <variant>

namespace seissol::io::instance::geometry {

class GeometryWriter {
  private:
  static std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter>
      getUnderlyingWriter(const std::string& name,
                          std::size_t localElementCount,
                          std::size_t dimension,
                          int targetDegree) {
    if (targetDegree < 0) {
      return mesh::XdmfWriter();
    } else {
      return mesh::VtkHdfWriter(name, localElementCount, dimension, targetDegree);
    }
  }

  public:
  GeometryWriter(const std::string& name,
                 std::size_t localElementCount,
                 std::size_t dimension,
                 int targetDegree)
      : underlying(getUnderlyingWriter(name, localElementCount, dimension, targetDegree)) {}

  template <typename F>
  void addPointProjector(F&& projector) {
    std::visit([&](auto& writer) { writer.addPointProjector(std::forward<F>(projector)); },
               underlying);
  }

  template <typename T, typename F>
  void addGeometryOutput(const std::string& name,
                         const std::vector<std::size_t>& dimensions,
                         F&& writerFunction) {
    if (order == 0) {
      // cell output
      std::visit(
          [&](auto& writer) {
            writer.template addCellData<T>(name, dimensions, std::forward<F>(writerFunction));
          },
          underlying);
    } else {
      // point output
      std::visit(
          [&](auto& writer) {
            using WriterT = std::decay_t<decltype(writer)>;
            if constexpr (std::is_same_v<WriterT, mesh::VtkHdfWriter>) {
              writer.template addPointData<T>(name, dimensions, std::forward<F>(writerFunction));
            } else {
              logError() << "High-order (point) output is not supported by the Xdmf writer.";
            }
          },
          underlying);
    }
  }

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() {
    return std::visit([&](auto& writer) { return writer.makeWriter(); }, underlying);
  }

  void addHook(const std::function<void(std::size_t, double)>& hook) {
    std::visit([&](auto& writer) { return writer.addHook(hook); }, underlying);
  }

  protected:
  int order;
  int refinement;
  std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter> underlying;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
