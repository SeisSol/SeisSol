// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_DOFSYNC_H_
#define SEISSOL_SRC_RESULTWRITER_DOFSYNC_H_

#include "Initializer/TimeStepping/Halo.h"
#include "Modules/Module.h"

namespace seissol::writer {

class DofSync : public seissol::Module {
  public:
  void setup(const initializer::MeshLayout& layout, LTS::Storage* storage);
  void syncDofs(double time);

  private:
  initializer::MeshLayout layout_;
  LTS::Storage* storage_{nullptr};
  double time_{-1};
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_DOFSYNC_H_
