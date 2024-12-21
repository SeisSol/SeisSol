# SPDX-FileCopyrightText: 2020 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

include(FindPackageHandleStandardArgs)

find_path(MEMKIND_INCLUDE_DIR memkind.h)
find_library(MEMKIND_LIBRARIES memkind)

find_package_handle_standard_args(Memkind DEFAULT_MSG
        MEMKIND_INCLUDE_DIR MEMKIND_LIBRARIES)

mark_as_advanced(MEMKIND_INCLUDE_DIR MEMKIND_LIBRARIES)
