cmake_minimum_required(VERSION 3.7)

project(SeisSol-RConv LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../../cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()

find_package(Eigen3 3.4 REQUIRED)
find_package(NetCDF REQUIRED)
set(CMAKE_INSTALL_RPATH ${NetCDF_LIBRARY_DIR})
find_package(PROJ4 QUIET)
if (NOT PROJ4_FOUND)
    find_package(PkgConfig QUIET)
    if (PkgConfig_FOUND)
        pkg_check_modules(PROJ4 QUIET proj>=4.9)
        if (PROJ4_FOUND)
            set(CMAKE_INSTALL_RPATH "${PROJ4_LIBRARY_DIRS};${CMAKE_INSTALL_RPATH}")
            set(PROJ4_LIBRARIES "${PROJ4_LINK_LIBRARIES}")
        endif()
    endif()
endif()

# rpath has to be defined before add_executable
# https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
add_executable(SeisSol-rconv
  src/main.cpp
  src/XMFWriter.cpp
  src/Map.cpp
  src/NRFWriter.cpp
  src/SRF.cpp)

target_include_directories(SeisSol-rconv PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../src" # SeisSol :(
   "${CMAKE_CURRENT_SOURCE_DIR}/../../../submodules/"
)

target_link_libraries(SeisSol-rconv PUBLIC Eigen3::Eigen)

target_include_directories(SeisSol-rconv PUBLIC ${NetCDF_INCLUDE_DIRS})
target_link_libraries(SeisSol-rconv PUBLIC ${NetCDF_LIBRARY})

if (PROJ4_FOUND)
    message(STATUS "Compiling with proj.4 support.")
    target_include_directories(SeisSol-rconv PUBLIC ${PROJ4_INCLUDE_DIRS})
    target_link_libraries(SeisSol-rconv PUBLIC ${PROJ4_LIBRARIES})
else ()
    target_compile_definitions(SeisSol-rconv PRIVATE -Dnoproj)
    message(WARNING "Compiling WITHOUT proj.4 support.")
endif()

find_package(HDF5 REQUIRED COMPONENTS C HL)
target_link_libraries(SeisSol-rconv PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})

set_target_properties(SeisSol-rconv PROPERTIES OUTPUT_NAME "rconv")

install(TARGETS SeisSol-rconv RUNTIME DESTINATION bin)
