// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_SDERIVLIBRARY_H_
#define SEISSOL_SRC_EXPR_SDERIVLIBRARY_H_

// The library of the full language of sderiv: modules a module imports by their names, `import
// stress`, whose definitions it reads as `stress.andersonian(...)`. Each is a module of the full
// language, of definitions and types only, and part of SeisSol, so that it is the same in every
// run of one version:
//
//   stress    the stress of a fault optimally oriented under friction, as easi's special
//             components compute it, to their bits: str_dip_slip_am, optimal, andersonian
//   rock      gassmann: the elastic material of a poroelastic one (De la Puente et al. 2008, (44))
//   rotation  matrix: the rotation of SeisSol's point sources by strike, dip and rake; rotate
//   motion    angles and rotated: the horizontal motion along every angle, for RotD100

#include <map>
#include <string>

namespace seissol::expr {

/// The modules of the library, by name: their sources.
[[nodiscard]] const std::map<std::string, std::string>& libraryModules();

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_SDERIVLIBRARY_H_
