// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_TOOLS_H_
#define SEISSOL_SRC_PROXY_TOOLS_H_

#include <cstdint>
#include <ctime>

auto getCycles() -> std::uint64_t;

auto deriveCyclesFromTime(double time) -> double;

void printHostname();

auto sec(struct timeval start, struct timeval end) -> double;

#endif // SEISSOL_SRC_PROXY_TOOLS_H_
