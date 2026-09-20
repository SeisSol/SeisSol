// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_HDF5_H_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_HDF5_H_

#include "Data.h"
#include "Instruction.h"

#include <memory>
#include <optional>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {

/**
 * @brief How a dataset grows from one write to the next.
 *
 * Which of the two growing modes is wanted depends on the reader, not on the data: a VTKHDF time
 * series wants one flat array that it slices with its step offsets, while a table indexed by both
 * time and entry wants the two as separate dimensions.
 */
enum class Append : std::uint8_t {
  //! Written once. Writing the same dataset a second time is an error.
  None,
  //! A leading dimension is added, and every write contributes one entry along it.
  Steps,
  //! The distributed dimension itself grows, so the data stays one flat array.
  Flat
};

//! @brief The name Append::value is serialised under, and back.
std::string appendName(Append append);
Append appendFromName(const std::string& name);
class Hdf5Location {
  public:
  explicit Hdf5Location(const std::string& longstring);

  Hdf5Location(const std::string& file,
               const std::vector<std::string>& groups,
               const std::optional<std::string>& dataset = std::optional<std::string>());

  explicit Hdf5Location(YAML::Node node);

  [[nodiscard]] std::string file() const;
  [[nodiscard]] std::vector<std::string> groups() const;
  [[nodiscard]] std::optional<std::string> dataset() const;
  [[nodiscard]] std::string infilePath() const;

  [[nodiscard]] std::optional<Hdf5Location> commonLocation(const Hdf5Location& other) const;

  YAML::Node serialize();

  private:
  std::string fileP_;
  std::vector<std::string> groupsP_;
  std::optional<std::string> datasetP_;
};

struct Hdf5AttributeWrite : public WriteInstruction {
  ~Hdf5AttributeWrite() override = default;
  Hdf5Location location;
  std::string name;
  std::shared_ptr<writer::DataSource> dataSource;

  YAML::Node serialize() override;

  Hdf5AttributeWrite(const Hdf5Location& location,
                     const std::string& name,
                     std::shared_ptr<writer::DataSource> dataSource);

  explicit Hdf5AttributeWrite(YAML::Node node);

  std::vector<std::shared_ptr<DataSource>> dataSources() override;
};

struct Hdf5DataWrite : public WriteInstruction {
  ~Hdf5DataWrite() override = default;
  Hdf5Location location;
  std::string name;
  std::shared_ptr<writer::DataSource> dataSource;
  std::shared_ptr<datatype::Datatype> targetType;
  Append append;
  int compress;

  Hdf5DataWrite(const Hdf5Location& location,
                const std::string& name,
                std::shared_ptr<writer::DataSource> dataSource,
                std::shared_ptr<datatype::Datatype> targetType,
                Append append = Append::None,
                int compress = 0);

  YAML::Node serialize() override;

  explicit Hdf5DataWrite(YAML::Node node);

  std::vector<std::shared_ptr<DataSource>> dataSources() override;
};

struct Hdf5LinkExternalWrite : public WriteInstruction {
  ~Hdf5LinkExternalWrite() override = default;
  Hdf5Location location;
  std::string name;
  Hdf5Location remote;

  YAML::Node serialize() override;

  Hdf5LinkExternalWrite(const Hdf5Location& location,
                        const std::string& name,
                        const Hdf5Location& remote);

  explicit Hdf5LinkExternalWrite(YAML::Node node);

  std::vector<std::shared_ptr<DataSource>> dataSources() override;
};
} // namespace seissol::io::writer::instructions

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_HDF5_H_
