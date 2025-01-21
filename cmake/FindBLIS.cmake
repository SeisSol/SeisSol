# SPDX-FileCopyrightText: 2020-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#  BLIS - BLAS-like Library Instantiation Software Framework
#  source code: https://github.com/flame/blis
#
#  NOTE: it is not an official cmake search file
#  author: Ravil Dorozhinskii
#  email: ravil.dorozhinskii@tum.de 
#
#  BLIS_FOUND        - system has Blis
#  BLIS_INCLUDE_DIRS - include directories for BLIS
#  BLIS_LIBRARIES    - libraries for BLIS
#
#  Additional env. variables that may be used by this module. 
#  They can change the default behaviour and
#  could to be set before calling find_package:
#
#  BLIS_DIR          - the root directory of the BLIS installation


include(FindPackageHandleStandardArgs)

find_path(BLIS_INCLUDE_DIRS blis.h
  HINTS ENV BLIS_DIR
  PATH_SUFFIXES include/blis
  DOC "Directory where the BLIS header files are located"
)

find_library(BLIS_LIBRARIES
  NAMES libblis.a libblis.so
  HINTS ENV BLIS_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the BLIS library is located"
)

if (BLIS_INCLUDE_DIRS AND BLIS_LIBRARIES)
  
  set(CMAKE_REQUIRED_INCLUDES ${BLIS_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${BLIS_LIBRARIES})
  set(CMAKE_REQUIRED_LINK_OPTIONS "-fopenmp")

  #Build and run test program
  include(CheckCXXSourceRuns)
  check_cxx_source_runs("
  #include <blis.h>

  int main()
  {
    dim_t m = 4, n = 5, k = 3;
    inc_t rsa = 1, csa = m;
    inc_t rsb = 1, csb = k;
    inc_t rsc = 1, csc= m;

    double c[m * n];
    double a[m * k];
    double b[k * n];

    double alpha = 1.0;
	  double beta = 1.0;

    bli_dgemm(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
	            m, n, k, 
              &alpha, a, rsa, csa, 
              b, rsb, csb, 
              &beta, c, rsc, csc);
  return 0;
  }
  " BLIS_TEST_RUNS)
endif()


find_package_handle_standard_args(BLIS
  BLIS_INCLUDE_DIRS 
  BLIS_LIBRARIES
  BLIS_TEST_RUNS)
