// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_MESH_VTKHDF_H_
#define SEISSOL_SRC_IO_INSTANCE_MESH_VTKHDF_H_

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Datatype/MPIType.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Instructions/Instruction.h"
#include "IO/Writer/Writer.h"
#include "Initializer/MemoryManager.h"
#include "utils/logger.h"

#include <IO/Datatype/Datatype.h>
#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <IO/Instance/Geometry/Typedefs.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <IO/Writer/Instructions/Instruction.h>
#include <IO/Writer/Writer.h>
#include <Initializer/MemoryManager.h>
#include <functional>
#include <memory>
#include <utils/logger.h>

namespace seissol::io::instance::mesh {
class VtkHdfWriter {
  public:
  VtkHdfWriter(const std::string& name,
               std::size_t localElementCount,
               geometry::Shape shape,
               std::size_t targetDegree);

  void addData(const std::string& name,
               const std::optional<std::string>& group,
               bool isConst,
               const std::shared_ptr<writer::DataSource>& data);

  template <typename F>
  void addPointProjector(F&& projector) {
    const auto data =
        writer::GeneratedBuffer::createElementwise<double>(localElementCount,
                                                           pointsPerElement,
                                                           std::vector<std::size_t>{3},
                                                           std::forward<F>(projector));

    addData("Points", std::optional<std::string>(), true, data);
  }

  template <typename T, typename F>
  void addPointData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    bool isConst,
                    F&& pointMapper) {
    const auto data = writer::GeneratedBuffer::createElementwise<T>(
        localElementCount, pointsPerElement, dimensions, std::forward<F>(pointMapper));
    addData(name, PointDataName, isConst, data);
    if (temporal) {
      const auto offset = isConst ? 0 : globalPointCount;
      addData(name,
              PointDataName + "Offsets",
              isConst,
              writer::WriteInline::createArray<uint64_t>({1}, {offset}));
    }
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   bool isConst,
                   F&& cellMapper) {
    const auto data = writer::GeneratedBuffer::createElementwise<T>(
        localElementCount, 1, dimensions, std::forward<F>(cellMapper));
    addData(name, CellDataName, isConst, data);
    if (temporal) {
      const auto offset = isConst ? 0 : globalElementCount;
      addData(name,
              CellDataName + "Offsets",
              false,
              writer::WriteInline::createArray<uint64_t>({1}, {offset}));
    }
  }

  template <typename T>
  void addFieldData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    bool isConst,
                    const std::vector<T>& data) {
    const auto datasource = writer::WriteInline::createArray(dimensions, data);
    addData(name, FieldDataName, isConst, datasource);
    if (temporal) {
      // TODO:
      const auto offset = isConst ? 0UL : 1UL;
      addData(name,
              FieldDataName + "Offsets",
              false,
              writer::WriteInline::createArray<uint64_t>({1}, {offset}));
    }
  }

  void addHook(const std::function<void(std::size_t, double)>& hook);

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter();

  private:
  std::string name;
  std::size_t localElementCount;
  std::size_t globalElementCount;
  std::size_t elementOffset{0};
  std::size_t localPointCount;
  std::size_t globalPointCount;
  std::size_t pointOffset;
  std::size_t pointsPerElement;
  std::vector<std::function<void(std::size_t, double)>> hooks;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, double)>>
      instructionsConst;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, const std::string&)>>
      instructionsConstLink;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, double)>>
      instructions;
  std::size_t type;
  std::size_t targetDegree;
  bool temporal{false};
  const static inline std::string GroupName = "VTKHDF";
  const static inline std::string FieldDataName = "FieldData";
  const static inline std::string CellDataName = "CellData";
  const static inline std::string PointDataName = "PointData";
};
} // namespace seissol::io::instance::mesh

#endif // SEISSOL_SRC_IO_INSTANCE_MESH_VTKHDF_H_
