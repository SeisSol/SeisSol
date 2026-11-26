// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// The core on its own (see seissol-core-link in tests/CMakeLists.txt): an executable of all the
// objects of the core and of the generated code, but of nothing above them, which links only if
// none of the objects uses the wave propagation, the features or the composition.
int main() { return 0; }
