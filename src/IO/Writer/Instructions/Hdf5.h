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
class Hdf5Location {
  public:
  Hdf5Location(const std::string& longstring);

  Hdf5Location(const std::string& file,
               const std::vector<std::string>& groups,
               const std::optional<std::string>& dataset = std::optional<std::string>());

  explicit Hdf5Location(YAML::Node node);

  [[nodiscard]] std::string file() const;
  [[nodiscard]] std::vector<std::string> groups() const;
  [[nodiscard]] std::optional<std::string> dataset() const;

  [[nodiscard]] std::optional<Hdf5Location> commonLocation(const Hdf5Location& other) const;

  YAML::Node serialize();

  private:
  std::string fileP;
  std::vector<std::string> groupsP;
  std::optional<std::string> datasetP;
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
  int compress;

  Hdf5DataWrite(const Hdf5Location& location,
                const std::string& name,
                std::shared_ptr<writer::DataSource> dataSource,
                std::shared_ptr<datatype::Datatype> targetType,
                int compress = 0);

  YAML::Node serialize() override;

  explicit Hdf5DataWrite(YAML::Node node);

  std::vector<std::shared_ptr<DataSource>> dataSources() override;
};
} // namespace seissol::io::writer::instructions

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_HDF5_H_
