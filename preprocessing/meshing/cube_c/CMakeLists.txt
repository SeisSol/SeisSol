cmake_minimum_required(VERSION 3.10)

# set the project name
project(cube_c)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

cmake_policy(SET CMP0074 NEW)


add_executable(cubeGenerator src/main.cpp)

#logging
set(LOG_LEVEL "warning" CACHE STRING "Log level for the code")
set(LOG_LEVEL_OPTIONS "debug" "info" "warning" "error")
set_property(CACHE LOG_LEVEL PROPERTY STRINGS ${LOG_LEVEL_OPTIONS})
if("${LOG_LEVEL}" STREQUAL "debug")
  target_compile_definitions(cubeGenerator PUBLIC LOG_LEVEL=3)
elseif("${LOG_LEVEL}" STREQUAL "info")
  target_compile_definitions(cubeGenerator PUBLIC LOG_LEVEL=2)
elseif("${LOG_LEVEL}" STREQUAL "warning")
  target_compile_definitions(cubeGenerator PUBLIC LOG_LEVEL=1)
elseif("${LOG_LEVEL}" STREQUAL "error")
  target_compile_definitions(cubeGenerator PUBLIC LOG_LEVEL=0)
endif()

#build and link libraries and executable
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../../cmake" ${CMAKE_MODULE_PATH})
find_package(NetCDF REQUIRED)
target_include_directories(cubeGenerator PUBLIC ${NetCDF_INCLUDE_DIRS})
target_link_libraries(cubeGenerator PUBLIC ${NetCDF_LIBRARY})

find_package(HDF5 REQUIRED
          COMPONENTS C HL)
target_include_directories(cubeGenerator PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(cubeGenerator PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})

find_package(MPI REQUIRED)
target_link_libraries(cubeGenerator PUBLIC MPI::MPI_C)

find_package(OpenMP REQUIRED)
target_link_libraries(cubeGenerator PUBLIC OpenMP::OpenMP_CXX)

#add some compiler specific flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_compile_options(cubeGenerator PUBLIC -fopenmp -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas>)
  target_link_libraries(cubeGenerator PUBLIC "-fopenmp")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options(cubeGenerator PUBLIC -qopenmp -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wall -w3 -diag-disable:remark>)
  target_link_libraries(cubeGenerator PUBLIC "-qopenmp")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp -Wall -Wextra -pedantic")
endif()

target_include_directories(cubeGenerator PUBLIC "${PROJECT_SOURCE_DIR}/../../../submodules")


