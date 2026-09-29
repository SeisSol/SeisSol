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

#include "Equations/Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Model/Quantities.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

void seissol::sourceterm::transformMomentTensor(const double localMomentTensor[3][3],
                                                const double localSolidVelocityComponent[3],
                                                double localPressureComponent,
                                                const double localFluidVelocityComponent[3],
                                                double strike,
                                                double dip,
                                                double rake,
                                                real* forceComponents) {
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

  // The source acts on the primary quantities only. They lead the quantity axis in every layout;
  // the memory variables of a fused anelastic layout follow them and receive no source.
  constexpr auto Groups = model::MaterialT::PrimaryGroups;
  static_assert(model::totalExtent(Groups) <= tensor::update::Size,
                "The primary quantities have to fit into the point source update.");
  constexpr auto StressKind = model::roleKind(Groups, model::FaceRole::Traction);
  static_assert(model::roleExtent(Groups, model::FaceRole::Traction) > 0 &&
                    (StressKind == model::QuantityKind::SymTensor2 ||
                     StressKind == model::QuantityKind::Scalar),
                "The moment tensor acts on a stress tensor or a scalar stress only.");

  std::fill(forceComponents, forceComponents + tensor::update::Size, 0);

  constexpr auto StressOffset = model::roleOffset(Groups, model::FaceRole::Traction);
  if constexpr (StressKind == model::QuantityKind::SymTensor2) {
    // Voigt order (xx, yy, zz, xy, yz, xz)
    forceComponents[StressOffset + 0] = m[0][0];
    forceComponents[StressOffset + 1] = m[1][1];
    forceComponents[StressOffset + 2] = m[2][2];
    forceComponents[StressOffset + 3] = m[0][1];
    forceComponents[StressOffset + 4] = m[1][2];
    forceComponents[StressOffset + 5] = m[0][2];
  } else {
    // a scalar stress, i.e. the acoustic pressure, takes the first diagonal entry
    forceComponents[StressOffset] = m[0][0];
  }

  for (std::size_t i = 0; i < 3; ++i) {
    forceComponents[model::MaterialT::VelocityOffset + i] = f[i];
  }

  // the poroelastic material adds the fluid pressure and the fluid velocity
  if constexpr (model::roleExtent(Groups, model::FaceRole::ExtraTraction) > 0) {
    forceComponents[model::roleOffset(Groups, model::FaceRole::ExtraTraction)] =
        localPressureComponent;
  }
  if constexpr (model::roleExtent(Groups, model::FaceRole::ExtraVelocity) > 0) {
    constexpr auto FluidOffset = model::roleOffset(Groups, model::FaceRole::ExtraVelocity);
    for (std::size_t i = 0; i < 3; ++i) {
      forceComponents[FluidOffset + i] = f[i + 3];
    }
  }
}
