// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_RUNNER_H_
#define SEISSOL_SRC_PROXY_RUNNER_H_

#include "Common.h"
namespace seissol::proxy {

auto runProxy(const ProxyConfig& config) -> ProxyOutput;

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_RUNNER_H_
