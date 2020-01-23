# NOTE: it is not an official cmake search file
# source: https://gist.github.com/scivision/5108cf6ab1515f581a84cd9ad1ef72aa
#
################################################################################
# NOTE: this file was originally developer by:
#
# \file      cmake/FindMKL.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Math Kernel Library from Intel
# \date      Thu 26 Jan 2017 02:05:50 PM MST
#
################################################################################

# Find the Math Kernel Library from Intel
#
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - MKL include files directories
#  MKL_LIBRARIES - The MKL libraries
#  MKL_INTERFACE_LIBRARY - MKL interface library
#  MKL_SEQUENTIAL_LAYER_LIBRARY - MKL sequential layer library
#  MKL_CORE_LIBRARY - MKL core library
#  MKL_COMPILER_DEFINITIONS - required compiler definitions for MKL
#
#  The environment variables MKLROOT and INTEL are used to find the library.
#  Everything else is ignored. If MKL is found "-DMKL_ILP64" is added to
#  MKL_COMPILER_DEFINITIONS.
#
#  Example usage:
#
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_include_directories(TARGET PUBLIC ${MKL_INCLUDE_DIRS})
#    target_link_libraries(TARGET PUBLIC ${MKL_LIBRARIES})
#    target_compile_definitions(TARGET PUBLIC ${MKL_COMPILER_DEFINITIONS})
#  endif()

# If already in cache, be silent
if (MKL_INCLUDE_DIRS AND MKL_LIBRARIES AND MKL_INTERFACE_LIBRARY AND
    MKL_SEQUENTIAL_LAYER_LIBRARY AND MKL_CORE_LIBRARY AND MKL_COMPILER_DEFINITIONS)
  set (MKL_FIND_QUIETLY TRUE)
endif()

if(NOT BUILD_SHARED_LIBS)
  set(INT_LIB "libmkl_intel_ilp64.a")
  set(SEQ_LIB "libmkl_sequential.a")
  set(THR_LIB "libmkl_intel_thread.a")
  set(COR_LIB "libmkl_core.a")
else()
  set(INT_LIB "mkl_intel_ilp64")
  set(SEQ_LIB "mkl_sequential")
  set(THR_LIB "mkl_intel_thread")
  set(COR_LIB "mkl_core")
endif()

find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS $ENV{MKLROOT}/include)

find_library(MKL_INTERFACE_LIBRARY
             NAMES ${INT_LIB}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_SEQUENTIAL_LAYER_LIBRARY
             NAMES ${SEQ_LIB}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_CORE_LIBRARY
             NAMES ${COR_LIB}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_SEQUENTIAL_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})

if (MKL_INCLUDE_DIR AND
    MKL_INTERFACE_LIBRARY AND
    MKL_SEQUENTIAL_LAYER_LIBRARY AND
    MKL_CORE_LIBRARY)

    if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
        NOT DEFINED ENV{CRAY_PRGENVGNU} AND
        NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
        NOT DEFINED ENV{CRAY_PRGENVINTEL})
      set(ABI "-m64")
    endif()

    set(MKL_COMPILER_DEFINITIONS "-DMKL_ILP64 ${ABI}")

else()

  set(MKL_INCLUDE_DIRS "")
  set(MKL_LIBRARIES "")
  set(MKL_INTERFACE_LIBRARY "")
  set(MKL_SEQUENTIAL_LAYER_LIBRARY "")
  set(MKL_CORE_LIBRARY "")
  set(MKL_COMPILER_DEFINITIONS "")

endif()

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES 
                                                  MKL_INCLUDE_DIRS 
                                                  MKL_INTERFACE_LIBRARY 
                                                  MKL_SEQUENTIAL_LAYER_LIBRARY 
                                                  MKL_CORE_LIBRARY
                                                  MKL_COMPILER_DEFINITIONS)

MARK_AS_ADVANCED(MKL_INCLUDE_DIRS 
                 MKL_LIBRARIES 
                 MKL_INTERFACE_LIBRARY 
                 MKL_SEQUENTIAL_LAYER_LIBRARY 
                 MKL_CORE_LIBRARY 
                 MKL_COMPILER_DEFINITIONS)