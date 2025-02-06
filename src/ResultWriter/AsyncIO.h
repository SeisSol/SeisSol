// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_ASYNCIO_H_
#define SEISSOL_SRC_RESULTWRITER_ASYNCIO_H_

#include "async/Dispatcher.h"

namespace seissol::io {

class AsyncIO : public async::Dispatcher {
  public:
  /**
   * @return False if this rank is an MPI executor that does not contribute to the
   *  computation.
   */
  bool init();
  void finalize();
};

} // namespace seissol::io

#endif // SEISSOL_SRC_RESULTWRITER_ASYNCIO_H_
