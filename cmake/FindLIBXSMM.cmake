include(FindPackageHandleStandardArgs)

find_path(LIBXSMM_INCLUDE_DIRS libxsmm.h
        PATH_SUFFIXES include
        DOC "Directory where libxsmm header files are located")
find_library(LIBXSMM_LIBRARIES
        NAMES xsmm
        PATH_SUFFICES lib
        DOC "Path to libxsmm library"
        )

if (NOT(
        LIBXSMM_INCLUDE_DIRS STREQUAL "LIBXSMM_INCLUDE_DIRS-NOTFOUND" OR
        LIBXSMM_LIBARIES STREQUAL "LIBXSMM_LIBRARIES-NOTFOUND"))
    add_library(LIBXSMM::LIBXSMM INTERFACE IMPORTED)
    set_target_properties(LIBXSMM::LIBXSMM PROPERTIES
            INTERFACE_LINK_LIBRARIES "${LIBXSMM_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${LIBXSMM_INCLUDE_DIRS}")
else()
    find_package(PkgConfig QUIET)
    pkg_check_modules(_LIBXSMM QUIET libxsmm IMPORTED_TARGET)
    if (_LIBXSMM_FOUND)
        # Forward the variables in a consistent way.
        set(LIBXSMM_FOUND "${_LIBXSMM_FOUND}")
        set(LIBXSMM_INCLUDE_DIRS "${_LIBXSMM_INCLUDE_DIRS}")
        set(LIBXSMM_LIBRARIES "${_LIBXSMM_LIBRARIES}")
        set(LIBXSMM_VERSION "${_LIBXSMM_VERSION}")
        if (NOT TARGET LIBXSMM::LIBXSMM)
            add_library(LIBXSMM::LIBXSMM INTERFACE IMPORTED)
            set_target_properties(LIBXSMM::LIBXSMM PROPERTIES
                    INTERFACE_LINK_LIBRARIES "PkgConfig::_LIBXSMM")
        endif()
    endif()
endif()
mark_as_advanced(LIBXSMM_INCLUDE_DIRS LIBXSMM_LIBRARIES)

find_package_handle_standard_args(LIBXSMM REQUIRED_VARS LIBXSMM_INCLUDE_DIRS LIBXSMM_LIBRARIES)