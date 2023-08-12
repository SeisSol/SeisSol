#pragma once

#include "Common/cellconfigconv.hpp"
#include "Initializer/ConfigFile.hpp"
#include "Initializer/ParameterDB.h"

#include <utility>
#include <vector>

namespace seissol::model {

// retrieves different materials
::std::pair<::std::vector<SupportedMaterials>, ::std::vector<::seissol::model::Plasticity>>
    queryMaterial(const initializer::CellConfigInfoMap& infoMap,
                  const initializers::CellToVertexArray& ctov,
                  bool fullCompute);

} // namespace seissol::model
