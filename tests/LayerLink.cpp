// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// A layer on its own (see seissol-core-link and the other seissol-*-link in tests/CMakeLists.txt):
// an executable of all the objects of the layer, of the layers below it and of the generated code,
// but of nothing else, which links only if none of the objects uses a layer above or beside them.
int main() { return 0; }
