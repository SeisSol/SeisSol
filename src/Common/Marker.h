// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_MARKER_H_
#define SEISSOL_SRC_COMMON_MARKER_H_

// host/device function markers

#if defined(__CUDACC__) || defined(__HIP__)
#define SEISSOL_DEVICE __device__
#define SEISSOL_HOST __host__
#else
#define SEISSOL_DEVICE
#define SEISSOL_HOST
#endif
#define SEISSOL_HOSTDEVICE SEISSOL_DEVICE SEISSOL_HOST

#endif // SEISSOL_SRC_COMMON_MARKER_H_
