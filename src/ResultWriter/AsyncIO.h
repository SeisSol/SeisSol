// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

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
