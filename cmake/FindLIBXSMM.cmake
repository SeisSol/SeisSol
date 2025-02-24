# SPDX-FileCopyrightText: 2022-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

if(LIBXSMM_INCLUDE_DIRS AND LIBXSMM_INCLUDE_DIRS AND LIBXSMM_LIBRARIES AND LIBXSMM_VERSION)
    set(LIBXSMM_FOUND TRUE)
else()
    find_path(LIBXSMM_INCLUDE_DIRS libxsmm.h
            PATH_SUFFIXES include
            DOC "Directory where libxsmm header files are located")
    find_library(LIBXSMM_LIBRARIES
            NAMES xsmm
            PATH_SUFFICES lib
            DOC "Path to libxsmm library"
            )

    find_library(LIBDL_LIBRARIES
            NAMES dl
            PATH_SUFFICES lib
            DOC "Path to libdl library"
            )

    list(APPEND LIBXSMM_LIBRARIES ${LIBDL_LIBRARIES})

    if (LIBXSMM_INCLUDE_DIRS STREQUAL "LIBXSMM_INCLUDE_DIRS-NOTFOUND" OR
        LIBXSMM_LIBARIES STREQUAL "LIBXSMM_LIBRARIES-NOTFOUND")
        find_package(PkgConfig QUIET)
        pkg_check_modules(_LIBXSMM QUIET libxsmm)
        if (_LIBXSMM_FOUND)
            set(LIBXSMM_INCLUDE_DIRS "${_LIBXSMM_INCLUDE_DIRS}")
            set(LIBXSMM_LIBRARIES "${_LIBXSMM_LIBRARIES}")
            set(LIBXSMM_VERSION "${_LIBXSMM_VERSION}")
        endif()
    endif()

    # Find version file from header
    if (LIBXSMM_INCLUDE_DIRS)
        file(READ "${LIBXSMM_INCLUDE_DIRS}/libxsmm_version.h" LIBXSMM_VERSION_HEADER)
        string(REGEX MATCH "#define LIBXSMM_CONFIG_VERSION_MAJOR [0-9]+" LIBXSMM_VERSION_MAJOR "${LIBXSMM_VERSION_HEADER}")
        string(REGEX MATCH "[0-9]+" LIBXSMM_VERSION_MAJOR "${LIBXSMM_VERSION_MAJOR}")
        string(REGEX MATCH "#define LIBXSMM_CONFIG_VERSION_MINOR [0-9]+" LIBXSMM_VERSION_MINOR "${LIBXSMM_VERSION_HEADER}")
        string(REGEX MATCH "[0-9]+" LIBXSMM_VERSION_MINOR "${LIBXSMM_VERSION_MINOR}")
        set(LIBXSMM_VERSION "${LIBXSMM_VERSION_MAJOR}.${LIBXSMM_VERSION_MINOR}")
        unset(LIBXSMM_VERSION_HEADER)
    endif()

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(LIBXSMM REQUIRED_VARS LIBXSMM_INCLUDE_DIRS LIBXSMM_LIBRARIES
            VERSION_VAR LIBXSMM_VERSION)

    mark_as_advanced(LIBXSMM_INCLUDE_DIRS LIBXSMM_LIBRARIES LIBXSMM_VERSION_MAJOR LIBXSMM_VERSION_MINOR LIBXSMM_VERSION)
endif()


if (LIBXSMM_FOUND AND NOT TARGET LIBXSMM::LIBXSMM)
    add_library(LIBXSMM::LIBXSMM INTERFACE IMPORTED)
    set_target_properties(LIBXSMM::LIBXSMM PROPERTIES
            INTERFACE_LINK_LIBRARIES "${LIBXSMM_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${LIBXSMM_INCLUDE_DIRS}")
endif()


