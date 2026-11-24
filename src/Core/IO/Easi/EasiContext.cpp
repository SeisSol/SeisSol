// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EasiContext.h"

#include <memory>

#ifdef USE_ASAGI
#include "Core/IO/Easi/AsagiReader.h"

#include <easi/util/AsagiReader.h>
#else
// Without ASAGI, the reader of ASAGI files is never made; the empty class completes its type, as
// the destructor of the pointer that holds it needs (see EasiReader.cpp).
namespace easi {
class AsagiReader {};
} // namespace easi
#endif

namespace seissol::data {

std::unique_ptr<easi::AsagiReader> makeAsagiReader([[maybe_unused]] const EasiContext& context) {
#ifdef USE_ASAGI
  if (context.asagi != nullptr) {
    return std::make_unique<seissol::asagi::AsagiReader>(*context.asagi);
  }
#endif
  return nullptr;
}

} // namespace seissol::data
