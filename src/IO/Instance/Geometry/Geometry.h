// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_

#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Instance/Mesh/Xdmf.h"

#include <variant>

namespace seissol::io::instance::geometry {

class GeometryWriter {
  private:
  static std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter> getUnderlyingWriter(
      const std::string& name, std::size_t localElementCount, Shape shape, int targetDegree) {
    if (targetDegree < 0) {
      return mesh::XdmfWriter(name, localElementCount, shape, 0);
    } else {
      return mesh::VtkHdfWriter(name, localElementCount, shape, targetDegree);
    }
  }

  public:
  GeometryWriter(const std::string& name,
                 std::size_t localElementCount,
                 Shape shape,
                 int targetDegree)
      : underlying(getUnderlyingWriter(name, localElementCount, shape, targetDegree)) {}

  template <typename F>
  void addPointProjector(F projector) {
    std::visit([&](auto& writer) { writer.addPointProjector(std::forward<F>(projector)); },
               underlying);
  }

  template <typename T, typename F>
  void addGeometryOutput(const std::string& name,
                         const std::vector<std::size_t>& dimensions,
                         bool isConst,
                         F writerFunction) {
    if (order <= 0) {
      // cell output
      std::visit(
          [&](auto& writer) {
            writer.template addCellData<T>(
                name, dimensions, isConst, std::forward<F>(writerFunction));
          },
          underlying);
    } else {
      // point output
      std::visit(
          [&](auto& writer) {
            writer.template addPointData<T>(
                name, dimensions, isConst, std::forward<F>(writerFunction));
          },
          underlying);
    }
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   bool isConst,
                   F writerFunction) {
    std::visit(
        [&](auto& writer) {
          writer.template addCellData<T>(
              name, dimensions, isConst, std::forward<F>(writerFunction));
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
  int order{};
  int refinement{};
  std::variant<mesh::VtkHdfWriter, mesh::XdmfWriter> underlying;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_GEOMETRY_H_
