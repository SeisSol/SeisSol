# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Cross-compiles with the GNU compilers of a target, named by its triplet in the environment
# variable SEISSOL_CROSS_TRIPLET (e.g. riscv64-linux-gnu), against the libraries that
# .ci/emulated/dependencies.sh builds into the prefix in SEISSOL_CROSS_PREFIX. The compilers are
# <triplet>-gcc, <triplet>-g++ and <triplet>-gfortran, unless CC, CXX and FC name others.
#
# CMake reads this file again for each of its checks, which see the same environment.

if ("$ENV{SEISSOL_CROSS_TRIPLET}" STREQUAL "" OR "$ENV{SEISSOL_CROSS_PREFIX}" STREQUAL "")
  message(FATAL_ERROR "Set SEISSOL_CROSS_TRIPLET and SEISSOL_CROSS_PREFIX to cross-compile.")
endif()

set(_triplet "$ENV{SEISSOL_CROSS_TRIPLET}")

set(CMAKE_SYSTEM_NAME Linux)
# as uname -m names the processor on the target, which is what the build compares it with
string(REGEX REPLACE "-.*" "" CMAKE_SYSTEM_PROCESSOR "${_triplet}")
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "powerpc64le")
  set(CMAKE_SYSTEM_PROCESSOR "ppc64le")
endif()

if (NOT "$ENV{CC}" STREQUAL "")
  set(CMAKE_C_COMPILER "$ENV{CC}")
else()
  set(CMAKE_C_COMPILER "${_triplet}-gcc")
endif()
if (NOT "$ENV{CXX}" STREQUAL "")
  set(CMAKE_CXX_COMPILER "$ENV{CXX}")
else()
  set(CMAKE_CXX_COMPILER "${_triplet}-g++")
endif()
if (NOT "$ENV{FC}" STREQUAL "")
  set(CMAKE_Fortran_COMPILER "$ENV{FC}")
else()
  set(CMAKE_Fortran_COMPILER "${_triplet}-gfortran")
endif()

# Libraries, headers and packages come from the prefix and from the libraries of the cross
# compiler; programs, such as the code generators, from the build machine.
set(CMAKE_FIND_ROOT_PATH "$ENV{SEISSOL_CROSS_PREFIX}" "/usr/${_triplet}")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# A check that runs what it compiles (try_run) does so with QEMU; the programs that the build and
# the tests run, binfmt_misc hands to QEMU on its own.
find_program(SEISSOL_CROSS_QEMU NAMES "qemu-${CMAKE_SYSTEM_PROCESSOR}-static"
                                      "qemu-${CMAKE_SYSTEM_PROCESSOR}")
if (SEISSOL_CROSS_QEMU)
  set(CMAKE_CROSSCOMPILING_EMULATOR "${SEISSOL_CROSS_QEMU}")
endif()
