##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2014-2015, SeisSol Group
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Builds the unit tests.

cmake_minimum_required( VERSION 2.6 )

option( ALIGNMENT         "Alignment in byte" )
option( GENERATED_KERNELS "Architecture and precision (snoarch, dnoarch, swsm, dwsm, ssnb, dsnb, sknc, dknc, shsw, dhsw)" )
option( CONVERGENCE_ORDER "Convergence order of the ADER-DG scheme" )

# enable advanced defintions
add_definitions( -DGENERATEDKERNELS )

# enable MPI
add_definitions( -DUSE_MPI )

# set alignemnt
add_definitions( -DALIGNMENT=${ALIGNMENT} )

# define the architecture
string( TOUPPER ${GENERATED_KERNELS} GENERATED_KERNELS )
add_definitions( -D${GENERATED_KERNELS} )

# set convergence order
add_definitions( -DCONVERGENCE_ORDER=${CONVERGENCE_ORDER} )

# set number of quantities
add_definitions( -DNUMBER_OF_QUANTITIES=9 )

# enable C++11
add_definitions( -std=c++11 -g )

# require mpi
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# add submodules
include_directories( ../
                     ../generated_code
                     ../../submodules
                   )

# add initializers
add_library( internal_state
             ../Initializer/InternalState.cpp )

add_library( mpi_layout
             ../Initializer/MpiLayout.cpp )

# add unit tests
find_package( CxxTest )
if(CXXTEST_FOUND)
  set( CXXTEST_USE_PYTHON TRUE )
  include_directories( ${CXXTEST_INCLUDE_DIR} )
  enable_testing()

  # unit tests for the lts setup
  CXXTEST_ADD_TEST( lts_setup
                    ltsSetup.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/Initializer/time_stepping/commonTestSuite.t.h )

  # unit tests for internal state
  CXXTEST_ADD_TEST( internal_state_unit_tests
                    internalState.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/Initializer/InternalStateTestSuite.t.h
                  )
  target_link_libraries( internal_state_unit_tests internal_state ${MPI_LIBRARIES} )

  # unit tests for the MPI layout
  CXXTEST_ADD_TEST( mpi_layout_unit_tests
                    mpiLayout.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/Initializer/MpiLayoutTestSuite.t.h
                  )
  target_link_libraries( mpi_layout_unit_tests mpi_layout internal_state ${MPI_LIBRARIES} )
endif()
