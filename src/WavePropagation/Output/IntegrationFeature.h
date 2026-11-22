// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_INTEGRATIONFEATURE_H_
#define SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_INTEGRATIONFEATURE_H_

#include "Core/Memory/Descriptor/LTS.h"
#include "Core/TimeStepping/Feature.h"

namespace seissol::kernels {

/**
 * The time integral of the unknowns of the cells, for the outputs of the wave field that write it,
 * as a feature of a run: a cell stage per cluster that integrates the cells after each correction,
 * if `integrate` (their storage holds the integrals then, see LTS::addTo()).
 */
class IntegrationFeature : public seissol::solver::Feature {
  public:
  explicit IntegrationFeature(bool integrate) : integrate_(integrate) {}

  [[nodiscard]] const char* name() const override { return "integrated output"; }

  void addStages(seissol::solver::TimeManager& timeManager, LTS::Storage& cells) override;

  private:
  bool integrate_;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_INTEGRATIONFEATURE_H_
