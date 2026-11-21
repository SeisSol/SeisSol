// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_RECEIVERPARAMETERS_H_
#define SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_RECEIVERPARAMETERS_H_

#include "Core/IO/WriterParameters.h"

#include <cstddef>
#include <string>

namespace seissol::initializer::parameters {

enum class ReceiverOutputFormat { Csv, Hdf5 };

struct ReceiverOutputParameters {
  ReceiverOutputFormat format{ReceiverOutputFormat::Csv};
  bool computeRotation{false};
  bool computeStrain{false};
  double writeInterval{VeryLongTime};
  double samplingInterval{0};
  std::string fileName;
  bool collectiveio{false};
  /**
   * @brief How far a storage chunk of the HDF5 receiver output reaches along the sample axis.
   *
   * Zero lets the writer decide. The chunking is settled when the file is created, while the
   * number of samples a write carries varies as soon as the sampling interval does not divide
   * the write interval, so which value pays off is a property of how the file is read afterwards.
   */
  std::size_t samplechunk{0};
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_RECEIVERPARAMETERS_H_
