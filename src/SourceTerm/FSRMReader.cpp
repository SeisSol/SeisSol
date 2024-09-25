/**
 * @file
 * This file is part of SeisSol.
 *
 * @section LICENSE
 * Copyright (c) 2023, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

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

template <size_t N, typename T>
static void readArrayOrZero(std::ifstream& filestream,
                            std::string& header,
                            const std::string& keyword,
                            T* data) {
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
  readArrayOrZero<3, size_t>(filestream, lineval, "dirNum", this->numInVelComponents);
  readArrayOrZero<3, real>(filestream, lineval, "velocity", this->solidVelocityComponent);
  readArrayOrZero<1, real>(filestream, lineval, "pressure", &this->pressureComponent);
  readArrayOrZero<3, real>(filestream, lineval, "fluid", this->fluidVelocityComponent);
  // (we've last read a header/comment line at this point here)

  // read faults
  filestream >> this->numberOfSources;
  this->centers.resize(this->numberOfSources);
  this->strikes.resize(this->numberOfSources);
  this->dips.resize(this->numberOfSources);
  this->rakes.resize(this->numberOfSources);
  this->onsets.resize(this->numberOfSources);
  this->areas.resize(this->numberOfSources);
  this->solidVel3cVector.resize(this->numberOfSources);
  for (auto& arr : this->solidVel3cVector) {
    arr.fill(0.0);
  }

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
    if (i < this->numInVelComponents[0]) {
      this->solidVel3cVector[i][0] = 1.0;
    } else if (i < this->numInVelComponents[0] + this->numInVelComponents[1]) {
      this->solidVel3cVector[i][1] = 1.0;
    } else if (i < this->numInVelComponents[0] + this->numInVelComponents[1] +
                       this->numInVelComponents[2]) {
      this->solidVel3cVector[i][2] = 1.0;
    } else if (this->numInVelComponents[0] + this->numInVelComponents[1] +
                   this->numInVelComponents[2] <
               1e-5)
      logInfo() << "Multi-components point source is nor specified in the FSRM file.\n";
    else
      throw std::runtime_error(
          "Total number of FSRM point sources does not match the sum of each component.");
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
