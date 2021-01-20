cmake_minimum_required(VERSION 3.10)

# set the project name
project(cube_c)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

cmake_policy(SET CMP0074 NEW)

add_executable(cubeGenerator src/main.cpp)

#build and link libraries and executable
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
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

target_compile_options(cubeGenerator PUBLIC -qopenmp)
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -qopenmp")

target_include_directories(cubeGenerator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules)


