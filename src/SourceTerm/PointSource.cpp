// SPDX-FileCopyrightText: 2015 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "PointSource.h"
#include <Kernels/Precision.h>
#include <Memory/MemoryAllocator.h>
#include <SourceTerm/Typedefs.h>
#include <algorithm>
#include <cmath>
#include <generated_code/tensor.h>

void seissol::sourceterm::transformMomentTensor(
    const double localMomentTensor[3][3],
    const double localSolidVelocityComponent[3],
    double localPressureComponent,
    const double localFluidVelocityComponent[3],
    double strike,
    double dip,
    double rake,
    seissol::memory::AlignedArray<real, PointSources::TensorSize>& forceComponents) {
  const double cstrike = std::cos(strike);
  const double sstrike = std::sin(strike);
  const double cdip = std::cos(dip);
  const double sdip = std::sin(dip);
  const double crake = std::cos(rake);
  const double srake = std::sin(rake);

  // Note, that R[j][i] = R_{ij} here.
  const double r[3][3] = {{crake * cstrike + cdip * srake * sstrike,
                           cdip * crake * sstrike - cstrike * srake,
                           sdip * sstrike},
                          {cdip * cstrike * srake - crake * sstrike,
                           srake * sstrike + cdip * crake * cstrike,
                           cstrike * sdip},
                          {-sdip * srake, -crake * sdip, cdip}};

  double m[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  // Calculate M_{ij} = R_{ki} * LM_{kl} * R_{lj}.
  // Note, again, that X[j][i] = X_{ij} here.
  // As M is symmetric, it is sufficient to calculate
  // (i,j) = (0,0), (1,0), (2,0), (1,1), (2,1), (2,2)
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          m[j][i] += r[i][k] * localMomentTensor[l][k] * r[j][l];
        }
      }
    }
  }
  double f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int j = 0; j < 3; ++j) {
    for (int k = 0; k < 3; ++k) {
      f[k] += r[k][j] * localSolidVelocityComponent[j];
      f[k + 3] += r[k][j] * localFluidVelocityComponent[j];
    }
  }

  std::fill(forceComponents.data(), forceComponents.data() + forceComponents.size(), 0);
  // Save in order (\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{xz}, u,
  // v, w, p, u_f, v_f, w_f)

  // TODO: prettify the code
  forceComponents[0] = m[0][0];
  if constexpr (tensor::momentFSRM::Size == 4) {
    forceComponents[1] = f[0];
    forceComponents[2] = f[1];
    forceComponents[3] = f[2];

    // unused for FSRM (but needed for NRF)
    forceComponents[4] = 0;
    forceComponents[5] = 0;
    forceComponents[6] = 0;
    forceComponents[7] = 0;
    forceComponents[8] = 0;
  } else {
    forceComponents[1] = m[1][1];
    forceComponents[2] = m[2][2];
    forceComponents[3] = m[0][1];
    forceComponents[4] = m[1][2];
    forceComponents[5] = m[0][2];
    forceComponents[6] = f[0];
    forceComponents[7] = f[1];
    forceComponents[8] = f[2];
    if constexpr (tensor::momentFSRM::Size >= 13) {
      forceComponents[9] = localPressureComponent;
      forceComponents[10] = f[3];
      forceComponents[11] = f[4];
      forceComponents[12] = f[5];
    }
  }
}
