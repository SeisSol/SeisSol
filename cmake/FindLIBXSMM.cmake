include(FindPackageHandleStandardArgs)

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
mark_as_advanced(LIBXSMM_INCLUDE_DIRS LIBXSMM_LIBRARIES)

find_package_handle_standard_args(LIBXSMM REQUIRED_VARS LIBXSMM_INCLUDE_DIRS LIBXSMM_LIBRARIES)