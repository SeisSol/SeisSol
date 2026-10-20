// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_ASYNCIO_H_
#define SEISSOL_SRC_RESULTWRITER_ASYNCIO_H_

#include <async/Dispatcher.h>
#include <string>
#include <utils/env.h>

namespace seissol::parallel {
class Pinning;
} // namespace seissol::parallel

namespace seissol::io {

/**
 * The dispatcher of the asynchronous output.
 *
 * Its settings are read from SEISSOL_ASYNC_* (generally: the prefix of the given environment,
 * followed by ASYNC_), and from the plain ASYNC_* names, which are deprecated.
 */
class AsyncIO : public async::Dispatcher {
  public:
  explicit AsyncIO(utils::Env& env);

  /**
   * Settles the mode AUTO: SYNC unless every rank leaves CPUs free for the threads of the output;
   * otherwise, ASYNC chooses between THREAD and NODE from how the ranks sit on the nodes.
   *
   * @return False if this rank is an MPI executor that does not contribute to the
   *  computation.
   */
  bool initDispatcher(const parallel::Pinning& pinning);
  void finalizeDispatcher();

  //! Whether the mode was AUTO, which initDispatcher has settled
  [[nodiscard]] bool automatic() const { return automatic_; }

  private:
  //! The prefix of the settings, for the messages
  std::string prefix_;

  bool automatic_{false};
};

} // namespace seissol::io

#endif // SEISSOL_SRC_RESULTWRITER_ASYNCIO_H_
