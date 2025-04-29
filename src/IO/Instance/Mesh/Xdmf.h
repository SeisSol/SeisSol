// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_
#define SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_

#include <IO/Writer/Writer.h>
#include <functional>
namespace seissol::io::instance::mesh {

class XdmfWriter {
  public:
  template <typename F>
  void addPointProjector(F projector) {}

  template <typename T, typename F>
  void addCellData(const std::string& name,
                   const std::vector<std::size_t>& dimensions,
                   F&& writerFunction) {}

  void addHook(const std::function<void(std::size_t, double)>& hook) {}

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() {
    return [](const std::string&, std::size_t, double) { return writer::Writer(); };
  }

  private:
};

} // namespace seissol::io::instance::mesh

#endif // SEISSOL_SRC_IO_INSTANCE_MESH_XDMF_H_
