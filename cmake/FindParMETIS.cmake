# FROM DUNE


# .. cmake_module::
#
#    Module that checks whether ParMETIS is available.
#
#    You may set the following variables to configure this modules behaviour:
#
#    :ref:`PARMETIS_ROOT`
#       Prefix where ParMETIS is installed.
#
#    :ref:`METIS_LIB_NAME`
#       Name of the METIS library (default: metis).
#
#    :ref:`PARMETIS_LIB_NAME`
#       Name of the ParMETIS library (default: parmetis).
#
#    :ref:`METIS_LIBRARY`
#       Full path of the METIS library.
#
#    :ref:`PARMETIS_LIBRARY`
#       Full path of the ParMETIS library
#
#    Sets the following variables:
#
#    :code:`PARMETIS_FOUND`
#       True if ParMETIS was found.
#
#    :code:`METIS_LIBRARY`
#       Full path of the METIS library.
#
#    :code:`PARMETIS_LIBRARY`
#       Full path of the ParMETIS library.
#
#    :code:`PARMETIS_LIBRARIES`
#       List of all libraries needed for linking with ParMETIS,
#
# .. cmake_variable:: PARMETIS_ROOT
#
#    You may set this variable to have :ref:`FindParMETIS` look
#    for the ParMETIS library and includes in the given path
#    before inspecting default system paths.
#
# .. cmake_variable:: PARMETIS_LIB_NAME
#
#    You may set this variable to specify the name of the ParMETIS
#    library that :ref:`FindParMETIS` looks for.
#
# .. cmake_variable:: PARMETIS_LIBRARY
#
#    You may set this variable to specify the full path to the ParMETIS
#    library, that should be used by :ref:`FindParMETIS`.
#


find_path(PARMETIS_INCLUDE_DIR parmetis.h
          PATHS ${PARMETIS_DIR} ${PARMETIS_ROOT}
          PATH_SUFFIXES include parmetis
          NO_DEFAULT_PATH
          DOC "Include directory of ParMETIS")
find_path(PARMETIS_INCLUDE_DIR parmetis.h
          PATH_SUFFIXES include parmetis)

set(METIS_LIB_NAME metis
    CACHE STRING "Name of the METIS library (default: metis).")
set(PARMETIS_LIB_NAME parmetis
    CACHE STRING "Name of the ParMETIS library (default: parmetis).")
set(METIS_LIBRARY METIS_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the METIS library")
set(PARMETIS_LIBRARY ParMETIS_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the ParMETIS library")

# check METIS and ParMETIS headers
include(CMakePushCheckState)
cmake_push_check_state() # Save variables
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_DUNE_INCLUDE_PATH} ${PARMETIS_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${MPI_DUNE_COMPILE_FLAGS}")
check_include_file(metis.h METIS_FOUND)
check_include_file(parmetis.h PARMETIS_FOUND)

# check whether metis.h is available
# to work around installation bug in ParMETIS 4.0.3
if(NOT METIS_FOUND)
  message(WARNING "metis.h is missing, you have to copy it manually next to parmetis.h")
endif()

if(PARMETIS_FOUND)
  set(ParMETIS_INCLUDE_PATH ${CMAKE_REQUIRED_INCLUDES})
  set(ParMETIS_COMPILE_FLAGS "${CMAKE_REQUIRED_FLAGS} -DENABLE_PARMETIS=1")

  # search METIS library
  find_library(METIS_LIBRARY metis
               PATHS ${PARMETIS_DIR} ${PARMETIS_ROOT}
               PATH_SUFFIXES lib
               NO_DEFAULT_PATH)
  find_library(METIS_LIBRARY metis)

  # search ParMETIS library
  find_library(PARMETIS_LIBRARY parmetis
               PATHS ${PARMETIS_DIR} ${PARMETIS_ROOT}
               PATH_SUFFIXES lib
               NO_DEFAULT_PATH)
  find_library(PARMETIS_LIBRARY parmetis)

  set(_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}") # do a backup
  # check ParMETIS library
  if(PARMETIS_LIBRARY)
    set(_PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARIES} ${MPI_DUNE_LIBRARIES})
    set(CMAKE_REQUIRED_LIBRARIES ${_PARMETIS_LIBRARIES} ${_CMAKE_REQUIRED_LIBRARIES})
    include(CheckFunctionExists)
    check_function_exists(ParMETIS_V3_PartKway HAVE_PARMETIS)
    if(NOT HAVE_PARMETIS)
      # Maybe we are using static scotch libraries. In this case we need to link
      # the other scotch libraries too. Let's make a best effort.
      # Get the path where ParMETIS_LIBRARY resides
      get_filename_component(_lib_root ${METIS_LIBRARY} DIRECTORY)
      # Search for additional libs only in this directory.
      # Otherwise we might find incompatible ones, e.g. for int instead of long
      find_library(PTSCOTCH_LIBRARY ptscotch PATHS ${_lib_root} "The PT-Scotch library."
        NO_DEFAULT_PATH)
      find_library(PTSCOTCHERR_LIBRARY ptscotcherr PATHS ${_lib_root} "The Scotch error library."
        NO_DEFAULT_PATH)
      if(PTSCOTCH_LIBRARY AND PTSCOTCHERR_LIBRARY)
        set(_PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${PTSCOTCH_LIBRARY}
          ${PTSCOTCHERR_LIBRARY} ${METIS_LIBRARIES} ${MPI_DUNE_LIBRARIES})
        set(CMAKE_REQUIRED_LIBRARIES ${_PARMETIS_LIBRARIES}
          ${_CMAKE_REQUIRED_LIBRARIES})
        unset(HAVE_PARMETIS CACHE)
        check_function_exists(ParMETIS_V3_PartKway HAVE_PARMETIS)
      endif()
    endif()
  endif()
    set(CMAKE_REQUIRED_LIBRARIES "${_CMAKE_REQUIRED_LIBRARIES}") # get backup
endif()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ParMETIS"
  DEFAULT_MSG
  PARMETIS_INCLUDE_DIR
  PARMETIS_LIBRARY
  HAVE_PARMETIS
)

mark_as_advanced(PARMETIS_INCLUDE_DIR METIS_LIBRARY PARMETIS_LIBRARY METIS_LIB_NAME PARMETIS_LIB_NAME)

#restore old values
cmake_pop_check_state()

if(PARMETIS_FOUND)
  set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR})
  set(PARMETIS_LIBRARIES "${_PARMETIS_LIBRARIES}"
      CACHE FILEPATH "ParMETIS libraries needed for linking")
  set(PARMETIS_LINK_FLAGS "${DUNE_MPI_LINK_FLAGS}"
      CACHE STRING "ParMETIS link flags")
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ParMETIS succeeded:\n"
    "Include directory: ${PARMETIS_INCLUDE_DIRS}\n"
    "Library directory: ${PARMETIS_LIBRARIES}\n\n")
  # deprecate versions < 4
  file(READ "${PARMETIS_INCLUDE_DIR}/parmetis.h" parmetisheader)
  string(REGEX MATCH "#define PARMETIS_MAJOR_VERSION[ ]+[0-9]+" versionMacroDef "${parmetisheader}")
  string(REGEX MATCH "[0-9]+" ParMetisMajorVersion "${versionMacroDef}")
  if("${versionMacroDef}" STREQUAL "" OR "${ParMetisMajorVersion}" LESS 4)
    message(AUTHOR_WARNING "Support for METIS older than version 4.x is deprecated in Dune 2.7")
  endif()
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of ParMETIS failed:\n"
    "Include directory: ${PARMETIS_INCLUDE_DIR}\n"
    "ParMETIS library directory: ${PARMETIS_LIBRARY}\n"
    "Header metis.h: ${METIS_FOUND}\n\n")
endif()

# register all ParMETIS related flags
if(PARMETIS_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_PARMETIS=1"
                              LIBRARIES "${PARMETIS_LIBRARIES}"
                              INCLUDE_DIRS "${PARMETIS_INCLUDE_DIRS}")
endif()

