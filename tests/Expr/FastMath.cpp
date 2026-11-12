// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Linked with an interpreter compiled with -ffast-math (tests/CMakeLists.txt): its self-test has to
// find that it does not compute as IEEE 754 says (the specification of sderiv, "Conformance").

#include "Expr/Interp.h"

#include <cstdio>
#include <string>

int main() {
  const std::string fault = seissol::expr::interpreterFault();
  if (fault.empty()) {
    std::printf("the self-test passed an interpreter built with -ffast-math\n");
    return 1;
  }
  std::printf("the self-test rejects an interpreter built with -ffast-math: %s\n", fault.c_str());
  return 0;
}
