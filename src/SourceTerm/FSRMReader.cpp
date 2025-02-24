// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FSRMReader.h"
#include <Kernels/Precision.h>
#include <cstddef>
#include <utils/logger.h>

#include <cassert>
#include <fstream>
#include <stdexcept>
#include <string>

// this code replicates the behavior of the corresponding FORTRAN code for legacy reasons. In
// particular, this reader is not programmed to be very fail-safe...

namespace {

template <size_t N>
void readArrayOrZero(std::ifstream& filestream,
                     std::string& header,
                     const std::string& keyword,
                     real* data) {
  if (header.find(keyword) != std::string::npos) {
    for (size_t i = 0; i < N; ++i) {
      filestream >> data[i];
    }
    std::getline(filestream, header); // end of line
    std::getline(filestream, header);
  } else {
    for (size_t i = 0; i < N; ++i) {
      data[i] = 0;
    }
  }
}

} // namespace

void seissol::sourceterm::FSRMSource::read(const std::string& filename) {
  logInfo() << "Reading FSRM point sources from file " << filename;

  std::ifstream filestream(filename);

  if (!filestream.good()) {
    throw std::runtime_error("FSRM point source file not found.");
  }

  std::string lineval;

  // comment
  std::getline(filestream, lineval);
  filestream >> this->momentTensor[0][0];
  filestream >> this->momentTensor[0][1];
  filestream >> this->momentTensor[0][2];
  filestream >> this->momentTensor[1][0];
  filestream >> this->momentTensor[1][1];
  filestream >> this->momentTensor[1][2];
  filestream >> this->momentTensor[2][0];
  filestream >> this->momentTensor[2][1];
  filestream >> this->momentTensor[2][2];
  std::getline(filestream, lineval); // end of line

  // comment
  std::getline(filestream, lineval);
  readArrayOrZero<3>(filestream, lineval, "velocity", this->solidVelocityComponent);
  readArrayOrZero<1>(filestream, lineval, "pressure", &this->pressureComponent);
  readArrayOrZero<3>(filestream, lineval, "fluid", this->fluidVelocityComponent);
  // (we've last read a header/comment line at this point here)

  // read faults
  filestream >> this->numberOfSources;
  this->centers.resize(this->numberOfSources);
  this->strikes.resize(this->numberOfSources);
  this->dips.resize(this->numberOfSources);
  this->rakes.resize(this->numberOfSources);
  this->onsets.resize(this->numberOfSources);
  this->areas.resize(this->numberOfSources);

  this->timeHistories.resize(this->numberOfSources);

  std::getline(filestream, lineval); // end of line

  // (read the comment before the actual data starts)
  std::getline(filestream, lineval);

  for (size_t i = 0; i < this->numberOfSources; ++i) {
    filestream >> this->centers[i](0);
    filestream >> this->centers[i](1);
    filestream >> this->centers[i](2);
    filestream >> this->strikes[i];
    filestream >> this->dips[i];
    filestream >> this->rakes[i];
    filestream >> this->areas[i];
    filestream >> this->onsets[i];
  }

  std::getline(filestream, lineval); // end of line

  // read samples
  std::getline(filestream, lineval);
  filestream >> this->timestep;
  filestream >> this->numberOfSamples;

  std::getline(filestream, lineval); // end of line

  // (the last comment)
  std::getline(filestream, lineval);

  for (size_t i = 0; i < this->numberOfSources; ++i) {
    this->timeHistories[i].resize(this->numberOfSamples);
    for (size_t j = 0; j < this->numberOfSamples; ++j) {
      filestream >> this->timeHistories[i][j];
    }
  }

  // and that's it!
}
