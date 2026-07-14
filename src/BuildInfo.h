// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_BUILDINFO_H_
#define SEISSOL_SRC_BUILDINFO_H_

#include <string>
namespace seissol {

struct BuildInfo {
  const static std::string VersionString;
  const static std::string CommitHash;
  const static std::string CommitYear;
  const static std::string CommitTimestamp;
  const static std::string SeisSolHostArch;
  const static std::string SeisSolDeviceArch;
  const static std::string SeisSolDeviceBackend;
};

} // namespace seissol

#endif // SEISSOL_SRC_BUILDINFO_H_
