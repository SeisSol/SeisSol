// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_
#define SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_

#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Instruction.h"
#include "IO/Writer/Writer.h"

#include <functional>
#include <memory>
namespace seissol::io::instance::mesh {

class XdmfWriter {
  public:
  struct WriteResult {
    std::string name;
    std::string type;
    std::shared_ptr<datatype::Datatype> datatype;
    std::string format;
    std::vector<std::size_t> offset;
    std::vector<std::size_t> dimensions;
    std::string location;
    std::shared_ptr<writer::instructions::WriteInstruction> instruction;
  };

  XdmfWriter(const std::string& name,
             std::size_t localElementCount,
             geometry::Shape shape,
             std::size_t targetDegree);

  void addData(const std::string& name,
               const std::string& type,
               bool isConst,
               std::size_t localCount,
               const std::shared_ptr<writer::DataSource>& data);

  template <typename F>
  void addPointProjector(F&& projector) {
    const auto data =
        writer::GeneratedBuffer::createElementwise<double>(localElementCount,
                                                           pointsPerElement,
                                                           std::vector<std::size_t>{3},
                                                           std::forward<F>(projector));

    addData("XYZ", "Geometry", true, localPointCount, data);
  }

  template <typename T, typename F>
  void addPointData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    bool isConst,
                    F&& pointMapper) {
    const auto data = writer::GeneratedBuffer::createElementwise<T>(
        localElementCount, pointsPerElement, dimensions, std::forward<F>(pointMapper));
    addData(name, "AttributeNode", isConst, localPointCount, data);
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   bool isConst,
                   F&& cellMapper) {
    const auto data = writer::GeneratedBuffer::createElementwise<T>(
        localElementCount, 1, dimensions, std::forward<F>(cellMapper));
    addData(name, "AttributeCell", isConst, localElementCount, data);
  }

  void addHook(const std::function<void(std::size_t, double)>& hook);

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter();

  private:
  std::string name;
  std::string type;
  bool binary{false};
  std::size_t localElementCount;
  std::size_t globalElementCount;
  std::size_t elementOffset{0};
  std::size_t localPointCount;
  std::size_t globalPointCount;
  std::size_t pointOffset;
  std::size_t pointsPerElement;
  std::size_t datasetCount{0};
  std::vector<std::function<void(std::size_t, double)>> hooks;
  std::vector<std::function<WriteResult(const std::string&, std::size_t)>> instructionsConst;
  std::vector<std::function<WriteResult(const std::string&, std::size_t)>> instructions;
};

} // namespace seissol::io::instance::mesh

#endif // SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_
