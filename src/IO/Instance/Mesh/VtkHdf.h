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

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace seissol::io::instance::mesh {

/**
 * @brief Which point each corner of each cell refers to.
 *
 * Without one, every cell writes its own copy of its corners. With one, the points are written
 * once and the cells index into them, which is what makes an order 0 output as large as the mesh
 * rather than as large as the mesh times the number of cells a vertex touches.
 */
struct VertexMap {
  std::size_t localPointCount{0};
  //! localElementCount * pointsPerElement entries, each below localPointCount
  std::vector<std::size_t> connectivity;
};
class VtkHdfWriter {
  public:
  VtkHdfWriter(const std::string& name,
               std::size_t localElementCount,
               geometry::Shape shape,
               std::size_t targetDegree,
               bool temporal,
               std::int32_t compress,
               bool constFile = false,
               std::optional<VertexMap> vertexMap = {});

  void addData(const std::string& name,
               const std::optional<std::string>& group,
               bool isConst,
               const std::shared_ptr<writer::DataSource>& data,
               bool attribute = false);

  /**
   * @brief The shape of bulk data in this file.
   *
   * In a time series the steps are concatenated along the same dimension the ranks are split
   * along -- one flat array that the step offsets slice -- so that dimension both moves between
   * ranks and grows. Data that does not change is written once and read again by every step.
   */
  [[nodiscard]] std::vector<writer::Dimension> bulkDimensions(const std::vector<std::size_t>& shape,
                                                              bool isConst) const {
    std::vector<writer::Dimension> result;
    result.push_back((!isConst && temporal_) ? writer::Dimension::distributedAppended()
                                             : writer::Dimension::distributed());
    for (const auto size : shape) {
      result.push_back(writer::Dimension::replicated(size));
    }
    return result;
  }

  /**
   * @brief Installs the source of the point coordinates.
   *
   * Without a vertex map the projector is called once per cell and fills its corners; with one it
   * is called once per point and fills that point.
   */
  template <typename F>
  void addPointProjector(F&& projector) {
    const auto data = writer::GeneratedBuffer::createElementwiseShaped<double>(
        pointSourceCount_, pointsPerSource_, bulkDimensions({3}, true), std::forward<F>(projector));

    addData("Points", std::optional<std::string>(), true, data);
  }

  template <typename T, typename F>
  void addPointData(const std::string& name,
                    const std::vector<std::size_t>& dimensions,
                    bool isConst,
                    F&& pointMapper) {
    const auto data =
        writer::GeneratedBuffer::createElementwiseShaped<T>(localElementCount_,
                                                            pointsPerElement_,
                                                            bulkDimensions(dimensions, isConst),
                                                            std::forward<F>(pointMapper));
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
    const auto data = writer::GeneratedBuffer::createElementwiseShaped<T>(
        localElementCount_, 1, bulkDimensions(dimensions, isConst), std::forward<F>(cellMapper));
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
    // Field data is not split across the ranks, so a time series grows it along a dimension of
    // its own rather than along a distributed one.
    std::vector<writer::Dimension> shape;
    if (!isConst && temporal_) {
      shape.push_back(writer::Dimension::appended(dimensions.empty() ? 1 : dimensions.front()));
      for (std::size_t i = 1; i < dimensions.size(); ++i) {
        shape.push_back(writer::Dimension::replicated(dimensions[i]));
      }
    } else {
      for (const auto size : dimensions) {
        shape.push_back(writer::Dimension::replicated(size));
      }
    }
    const auto datasource = writer::WriteInline::createShaped(shape, data);
    addData(name, FieldDataName, isConst, datasource);
    if (temporal_) {
      const auto tuples = dimensions.empty() ? 1 : dimensions.front();
      std::size_t components = 1;
      for (std::size_t i = 1; i < dimensions.size(); ++i) {
        components *= dimensions[i];
      }
      addStepOffset(name, FieldDataName + "Offsets", isConst ? 0 : tuples);
      addStepFieldDataSize(name, components, tuples);
    }
  }

  /**
   * @brief Adds field data whose values are produced anew for every step.
   *
   * @p provider is handed the number of the step within the file and its time, and returns the
   * tuples of that step laid out tuple-major; how many it returns may differ from one step to the
   * next. @p components is what a single tuple holds, and stays the same for the run.
   *
   * The three datasets that describe a step -- its values, where they start, and how many there
   * are -- share the count the provider just produced. They are built in one pass over the
   * instructions and in this order, so the two that only describe the values see what the first
   * one wrote.
   */
  template <typename T, typename F>
  void addStepFieldData(const std::string& name,
                        const std::vector<std::size_t>& components,
                        F&& provider) {
    std::size_t perTuple = 1;
    for (const auto size : components) {
      perTuple *= size;
    }

    //! What the step being built holds, and where it starts.
    struct StepExtent {
      std::size_t tuples{0};
      std::size_t start{0};
    };
    const auto extent = std::make_shared<StepExtent>();
    const bool temporal = temporal_;

    const std::vector<std::string> dataGroups{GroupName, FieldDataName};
    instructions_.emplace_back([=, provider = std::forward<F>(provider)](
                                   const std::string& filename, std::size_t step, double time) {
      const std::vector<T> values = std::invoke(provider, step, time);
      if (perTuple == 0 || values.size() % perTuple != 0) {
        logError() << "The field data" << name << "produced" << values.size()
                   << "values, which is not a whole number of tuples of" << perTuple
                   << "components.";
      }
      extent->tuples = values.size() / perTuple;

      std::vector<writer::Dimension> shape;
      shape.push_back(temporal ? writer::Dimension::appended(extent->tuples)
                               : writer::Dimension::replicated(extent->tuples));
      for (const auto size : components) {
        shape.push_back(writer::Dimension::replicated(size));
      }

      const auto data = writer::WriteInline::createShaped<T>(shape, values);
      return std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, dataGroups), name, data, data->datatype());
    });

    if (temporal_) {
      const std::vector<std::string> offsetGroups{GroupName, StepsName, FieldDataName + "Offsets"};
      instructions_.emplace_back(
          [=](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
            const auto start = extent->start;
            extent->start += extent->tuples;
            const auto data = writer::WriteInline::createShaped<std::uint64_t>(
                {writer::Dimension::appended(1)}, {static_cast<std::uint64_t>(start)});
            return std::make_shared<writer::instructions::Hdf5DataWrite>(
                writer::instructions::Hdf5Location(filename, offsetGroups),
                name,
                data,
                data->datatype());
          });

      const std::vector<std::string> sizeGroups{GroupName, StepsName, FieldDataName + "Sizes"};
      instructions_.emplace_back([=](const std::string& filename,
                                     std::size_t /*counter*/,
                                     double /*time*/) {
        const auto data = writer::WriteInline::createShaped<std::int64_t>(
            {writer::Dimension::appended(1), writer::Dimension::replicated(2)},
            {static_cast<std::int64_t>(perTuple), static_cast<std::int64_t>(extent->tuples)});
        return std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, sizeGroups), name, data, data->datatype());
      });
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

  /**
   * @brief Records the component and tuple count of a field data array for every step.
   *
   * Without it a reader assumes one tuple per step and the largest component count it finds,
   * which is only right for a scalar.
   */
  void addStepFieldDataSize(const std::string& name, std::size_t components, std::size_t tuples);

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
  //! Where this rank's entries start in the concatenated connectivity, which is not the point
  //! offset once the points are shared.
  std::size_t connectivityOffset_{0};
  //! How the point projector is called; see addPointProjector.
  std::size_t pointSourceCount_{0};
  std::size_t pointsPerSource_{0};
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
