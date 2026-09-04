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
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Instructions/Instruction.h"
#include "IO/Writer/Writer.h"
#include "utils/logger.h"

#include <functional>
#include <memory>

namespace seissol::io::instance::mesh {
class VtkHdfWriter {
  public:
  VtkHdfWriter(const std::string& name,
               std::size_t localElementCount,
               geometry::Shape shape,
               std::size_t targetDegree,
               bool temporal,
               std::int32_t compress,
               bool constFile = false);

  void addData(const std::string& name,
               const std::optional<std::string>& group,
               bool isConst,
               const std::shared_ptr<writer::DataSource>& data,
               bool attribute = false);

  template <typename F>
  void addPointProjector(F&& projector) {
    const auto data =
        writer::GeneratedBuffer::createElementwise<double>(localElementCount_,
                                                           pointsPerElement_,
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
        localElementCount_, pointsPerElement_, dimensions, std::forward<F>(pointMapper));
    addData(name, PointDataName, isConst, data);
    if (temporal_) {
      addStepOffset(name, PointDataName + "Offsets", isConst ? 0 : globalPointCount_);
    }
  }

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   bool isConst,
                   F&& cellMapper) {
    const auto data = writer::GeneratedBuffer::createElementwise<T>(
        localElementCount_, 1, dimensions, std::forward<F>(cellMapper));
    addData(name, CellDataName, isConst, data);
    if (temporal_) {
      addStepOffset(name, CellDataName + "Offsets", isConst ? 0 : globalElementCount_);
    }
  }

  template <typename T>
  void addFieldData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    bool isConst,
                    const std::vector<T>& data) {
    const auto datasource = writer::WriteInline::createArray(dimensions, data);
    addData(name, FieldDataName, isConst, datasource);
    if (temporal_) {
      // one tuple per step, which is what a reader assumes in the absence of FieldDataSizes
      addStepOffset(name, FieldDataName + "Offsets", isConst ? 0 : 1);
    }
  }

  /**
   * @brief Adds one of the offsets of the Steps group: the position in the flattened array at
   * which the data of a step begins.
   *
   * @p perStep is how far the offset advances from one step to the next; zero for data that is
   * written once and read again by every step.
   */
  void addStepOffset(const std::string& name,
                     const std::optional<std::string>& group,
                     std::size_t perStep);

  void addHook(const std::function<void(std::size_t, double)>& hook);

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter();

  private:
  std::string name_;
  std::size_t localElementCount_;
  std::size_t globalElementCount_;
  std::size_t elementOffset_{0};
  std::size_t localPointCount_;
  std::size_t globalPointCount_;
  std::size_t pointOffset_;
  std::size_t pointsPerElement_;
  std::vector<std::function<void(std::size_t, double)>> hooks_;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, std::size_t, double)>>
      instructionsConst_;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, const std::string&)>>
      instructionsConstLink_;
  std::vector<std::function<std::shared_ptr<writer::instructions::WriteInstruction>(
      const std::string&, std::size_t, double)>>
      instructions_;
  std::size_t type_;
  std::size_t targetDegree_;
  //! Write the data that does not change into a file of its own, and link to it from every
  //! snapshot instead of repeating it.
  bool constFile_{false};
  bool temporal_{false};
  int32_t compress_{0};
  const static inline std::string GroupName = "VTKHDF";
  const static inline std::string StepsName = "Steps";
  const static inline std::string FieldDataName = "FieldData";
  const static inline std::string CellDataName = "CellData";
  const static inline std::string PointDataName = "PointData";
};
} // namespace seissol::io::instance::mesh

#endif // SEISSOL_SRC_IO_INSTANCE_MESH_VTKHDF_H_
