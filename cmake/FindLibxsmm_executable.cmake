# SPDX-FileCopyrightText: 2020-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#  libxsmm - Library targeting Intel Architecture for specialized 
#  dense and sparse matrix operations, and deep learning primitives
#  source code: https://github.com/hfp/libxsmm
#
#  NOTE: it is not an official cmake search file
#  author: Ravil Dorozhinskii
#  email: ravil.dorozhinskii@tum.de 
#
#  Libxsmm_executable_FOUND - system has Libxsmm_executable
#
#  Additional env. variables that may be used by this module. 
#  They can change the default behaviour and
#  could to be set before calling find_package:
#
#  LIBXSMM_DIR - the root directory of the BLIS installation


include(FindPackageHandleStandardArgs)

find_program(Libxsmm_executable_PROGRAM libxsmm_gemm_generator
  HINTS ENV LIBXSMM_DIR
  PATH_SUFFIXES bin
  DOC "Directory where the libxsmm binary file is located"
)

find_package_handle_standard_args(Libxsmm_executable
  REQUIRED_VARS Libxsmm_executable_PROGRAM)