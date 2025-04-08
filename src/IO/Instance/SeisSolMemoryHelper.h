// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_SEISSOLMEMORYHELPER_H_
#define SEISSOL_SRC_IO_INSTANCE_SEISSOLMEMORYHELPER_H_

#include <IO/Datatype/Inference.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <IO/Writer/Writer.h>
#include <Initializer/MemoryManager.h>
#include <Memory/Tree/LTSTree.h>
#include <type_traits>
namespace seissol::io::instance {

template <typename T>
std::shared_ptr<writer::DataSource> variableBuffer(initializer::LTSTree& tree,
                                                   initializer::Variable<T>& variable,
                                                   bool unwrapArray = false) {
  if (unwrapArray) {
    auto dimensions = datatype::arrayTotalExtent<T>();
    size_t subcount = 1;
    for (size_t dim : dimensions) {
      subcount *= dim;
    }
    return writer::WriteBuffer::create(tree.var(variable),
                                       subcount * tree.getNumberOfCells(variable.mask),
                                       datatype::inferDatatype<std::remove_all_extents_t<T>>());
  } else {
    return writer::WriteBuffer::create(tree.var(variable), tree.getNumberOfCells(variable.mask));
  }
}

// (bucketBuffer has, as seen so far, not been needed; hence it has not been implemented at this
// point of time)

} // namespace seissol::io::instance

#endif // SEISSOL_SRC_IO_INSTANCE_SEISSOLMEMORYHELPER_H_
