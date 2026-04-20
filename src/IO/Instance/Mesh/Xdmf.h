// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_
#define SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_

#include "IO/Instance/Geometry/Typedefs.h"
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
             std::size_t targetDegree,
             bool binary,
             int32_t compress);

  void addData(const std::string& name,
               const std::string& type,
               bool isConst,
               std::size_t localCount,
               const std::shared_ptr<writer::DataSource>& data);

  template <typename F>
  void addPointProjector(F&& projector) {
    const auto data =
        writer::GeneratedBuffer::createElementwise<double>(localElementCount_,
                                                           pointsPerElement_,
                                                           std::vector<std::size_t>{3},
                                                           std::forward<F>(projector));

    addData("XYZ", "Geometry", true, localPointCount_, data);
  }

  template <typename T, typename F>
  void addPointData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    bool isConst,
                    F&& pointMapper) {
    const auto data = writer::GeneratedBuffer::createElementwise<T>(
        localElementCount_, pointsPerElement_, dimensions, std::forward<F>(pointMapper));
    addData(name, "AttributeNode", isConst, localPointCount_, data);
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   bool isConst,
                   F&& cellMapper) {
    const auto data = writer::GeneratedBuffer::createElementwise<T>(
        localElementCount_, 1, dimensions, std::forward<F>(cellMapper));
    addData(name, "AttributeCell", isConst, localElementCount_, data);
  }

  void addHook(const std::function<void(std::size_t, double)>& hook);

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter();

  private:
  std::string name_;
  std::string type_;
  bool binary_{false};
  int32_t compress_{0};
  std::size_t localElementCount_;
  std::size_t globalElementCount_;
  std::size_t elementOffset_{0};
  std::size_t localPointCount_;
  std::size_t globalPointCount_;
  std::size_t pointOffset_;
  std::size_t pointsPerElement_;
  std::size_t datasetCount_{0};
  std::vector<std::function<void(std::size_t, double)>> hooks_;
  std::vector<std::function<WriteResult(const std::string&, std::size_t)>> instructionsConst_;
  std::vector<std::function<WriteResult(const std::string&, std::size_t)>> instructions_;
};

} // namespace seissol::io::instance::mesh

#endif // SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_
