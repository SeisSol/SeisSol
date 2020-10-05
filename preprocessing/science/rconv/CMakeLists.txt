cmake_minimum_required(VERSION 3.7)

project(SeisSol-RConv LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../../cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()


add_executable(SeisSol-rconv
  src/main
  src/XMFWriter.cpp
  src/Map.cpp
  src/NRFWriter.cpp
  src/SRF.cpp)

target_include_directories(SeisSol-rconv PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../src" # SeisSol :(
   "${CMAKE_CURRENT_SOURCE_DIR}/../../../submodules/"
   "${CMAKE_CURRENT_SOURCE_DIR}/../../../submodules/eigen3"

)

find_package(NetCDF REQUIRED)
target_include_directories(SeisSol-rconv PUBLIC ${NetCDF_INCLUDE_DIRS})
target_link_libraries(SeisSol-rconv PUBLIC ${NetCDF_LIBRARY})

find_package(PROJ4)
target_include_directories(SeisSol-rconv PUBLIC ${PROJ4_INCLUDE_DIRS})
target_link_libraries(SeisSol-rconv PUBLIC ${PROJ4_LIBRARIES})

if(NOT PROJ4_FOUND)
    target_compile_definitions(SeisSol-rconv PRIVATE -Dnoproj)
    message(WARNING "Compiling without proj.4 support.")
endif()


set_target_properties(SeisSol-rconv PROPERTIES OUTPUT_NAME "rconv")
