# Try to find ParHIP

# The following definitions are added on success
#
#  ParHIP_FOUND - ParHIP was found
#  ParHIP_INCLUDE_DIR - ParHIP include directory
#  ParHIP_LIBRARY - ParHIP loader library
#
# and the following imported target:
#
#  ParHIP::ParHIP - ParHIP library
#
# The followings hints may be passed in the environment:
#
# ParHIP_ROOT
#

if(ParHIP_INCLUDE_DIR AND ParHIP_LIBRARY)
    set(ParHIP_FOUND TRUE)
else()
    find_path(ParHIP_INCLUDE_DIR NAMES parhip_interface.h
        HINTS
            ENV ParHIP_ROOT
            ENV CPATH
        PATH_SUFFIXES
            include
    )

    find_library(ParHIP_LIBRARY NAMES parhip_interface
        HINTS
        ENV ParHIP_ROOT
        ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
    )

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(ParHIP DEFAULT_MSG ParHIP_INCLUDE_DIR ParHIP_LIBRARY)

    mark_as_advanced(ParHIP_INCLUDE_DIR ParHIP_LIBRARY)
endif()

if(ParHIP_FOUND AND NOT TARGET ParHIP::ParHIP)
    add_library(ParHIP::ParHIP INTERFACE IMPORTED)
    set_target_properties(ParHIP::ParHIP PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${ParHIP_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${ParHIP_LIBRARY}")
endif()
