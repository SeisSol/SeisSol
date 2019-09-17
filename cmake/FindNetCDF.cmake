#[==[
Provides the following variables:

  * `NetCDF_FOUND`: Whether NetCDF was found or not.
  * `NetCDF_INCLUDE_DIRS`: Include directories necessary to use NetCDF.
  * `NetCDF_LIBRARIES`: Libraries necessary to use NetCDF.
  * `NetCDF_VERSION`: The version of NetCDF found.
  * `NetCDF::NetCDF`: A target to use with `target_link_libraries`.
#]==]

# Stolen from VTK with some small modifications
# https://github.com/Kitware/VTK/blob/master/CMake/FindNetCDF.cmake
# (BSD 3-Clause License)

#[[=========================================================================

  Program:   Visualization Toolkit
  Module:    Copyright.txt

Copyright (c) 1993-2015 Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================]]
# Note: The version from VTK tries to first find a CMake built version of CMake and then
# a pkgConf built one. This doesn't work on all systems with custom built versions of
# CMake + HDF5.
# Thus, first try to find it manually, then fall back on CMake/PkgConfig

find_path(NetCDF_INCLUDE_DIR
  NAMES netcdf.h
  DOC "netcdf include directories")
mark_as_advanced(NetCDF_INCLUDE_DIR)

find_library(NetCDF_LIBRARY
  NAMES netcdf
  DOC "netcdf library")
mark_as_advanced(NetCDF_LIBRARY)

if (NetCDF_INCLUDE_DIR)
  file(STRINGS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
    REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
  string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
  set(NetCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
  unset(_netcdf_version_major)
  unset(_netcdf_version_minor)
  unset(_netcdf_version_patch)
  unset(_netcdf_version_note)
  unset(_netcdf_version_lines)
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
  REQUIRED_VARS NetCDF_LIBRARY NetCDF_INCLUDE_DIR
  VERSION_VAR NetCDF_VERSION)

if (NetCDF_FOUND)
  set(NetCDF_INCLUDE_DIRS "${NetCDF_INCLUDE_DIR}")
  set(NetCDF_LIBRARIES "${NetCDF_LIBRARY}")

  if (NOT TARGET NetCDF::NetCDF)
    add_library(NetCDF::NetCDF UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF PROPERTIES
      IMPORTED_LOCATION "${NetCDF_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIR}")
  endif ()
else()
  # Try to find a CMake-built NetCDF.
  find_package(netCDF CONFIG QUIET)
  if (netCDF_FOUND)
    # Forward the variables in a consistent way.
    set(NetCDF_FOUND "${netCDF_FOUND}")
    set(NetCDF_INCLUDE_DIRS "${netCDF_INCLUDE_DIR}")
    set(NetCDF_LIBRARIES "${netCDF_LIBRARIES}")
    set(NetCDF_VERSION "${NetCDFVersion}")
    if (NOT TARGET NetCDF::NetCDF)
      add_library(NetCDF::NetCDF INTERFACE IMPORTED)
      set_target_properties(NetCDF::NetCDF PROPERTIES
        INTERFACE_LINK_LIBRARIES "${NetCDF_LIBRARIES}")
    endif ()
    # Skip the rest of the logic in this file.
    return ()
  endif ()

  find_package(PkgConfig QUIET)
  if (PkgConfig_FOUND)
    pkg_check_modules(_NetCDF QUIET netcdf IMPORTED_TARGET)
    if (_NetCDF_FOUND)
      # Forward the variables in a consistent way.
      set(NetCDF_FOUND "${_NetCDF_FOUND}")
      set(NetCDF_INCLUDE_DIRS "${_NetCDF_INCLUDE_DIRS}")
      set(NetCDF_LIBRARIES "${_NetCDF_LIBRARIES}")
      set(NetCDF_VERSION "${_NetCDF_VERSION}")
      if (NOT TARGET NetCDF::NetCDF)
        add_library(NetCDF::NetCDF INTERFACE IMPORTED)
        set_target_properties(NetCDF::NetCDF PROPERTIES
          INTERFACE_LINK_LIBRARIES "PkgConfig::_NetCDF")
      endif ()
      # Skip the rest of the logic in this file.
      return ()
    endif ()
  endif ()
endif ()
