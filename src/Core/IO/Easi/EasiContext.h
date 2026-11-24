// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_IO_EASI_EASICONTEXT_H_
#define SEISSOL_SRC_CORE_IO_EASI_EASICONTEXT_H_

#include <memory>

namespace easi {
class AsagiReader;
} // namespace easi

namespace seissol::asagi {
class AsagiModule;
} // namespace seissol::asagi

namespace seissol::data {

/**
 * How a run reads its easi models: the files of ASAGI with the settings of ASAGI of the process
 * (see asagi::AsagiModule), if the run has them; otherwise (e.g. in tests, or without ASAGI), as
 * easi opens them by itself.
 */
struct EasiContext {
  asagi::AsagiModule* asagi{nullptr};
};

/// The reader of the files of ASAGI for a parser of easi, as `context` says; none if easi opens
/// them by itself.
std::unique_ptr<easi::AsagiReader> makeAsagiReader(const EasiContext& context);

} // namespace seissol::data

#endif // SEISSOL_SRC_CORE_IO_EASI_EASICONTEXT_H_
