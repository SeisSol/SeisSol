// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_RESULTWRITER_CITATION_H_
#define SEISSOL_SRC_RESULTWRITER_CITATION_H_

#include <string>
#include <vector>

namespace seissol {

enum class CitationFeature { Core };

class CitationHandler {
  public:
  void enable(seissol::CitationFeature feature);
  void print(const std::string& prefix);

  private:
  std::vector<CitationFeature> features{CitationFeature::Core};
};

} // namespace seissol
#endif // SEISSOL_SRC_RESULTWRITER_CITATION_H_
