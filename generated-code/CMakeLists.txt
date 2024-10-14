# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/generated-code")

set(GENERATED_FILES_FOR_SEISSOL generated-code/subroutine.h
                                generated-code/tensor.cpp
                                generated-code/subroutine.cpp
                                generated-code/gpulike_subroutine.cpp
                                generated-code/tensor.h
                                generated-code/init.cpp
                                generated-code/init.h
                                generated-code/kernel.h
                                generated-code/kernel.cpp)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/generated-code)

if(PREMULTIPLY_FLUX)
  set(PREMULTIPLY_FLUX_VALUE "--enable_premultiply_flux")
else()
  set(PREMULTIPLY_FLUX_VALUE "--disable_premultiply_flux")
endif()

# TODO: list all Python files automatically?

add_custom_command(
  COMMAND
  "${Python3_EXECUTABLE}" "generate.py"
     "--equations" ${EQUATIONS}
     "--matricesDir" matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/generated-code
     "--host_arch" ${HOST_ARCH_STR}
     "--device_arch" ${DEVICE_ARCH_STR}
     "--device_backend" ${DEVICE_BACKEND}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" ${MEMORY_LAYOUT}
     "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
     "--gemm_tools" ${GEMM_TOOLS_LIST}
     "--drQuadRule" ${DR_QUAD_RULE}
     "--executable_libxsmm=${Libxsmm_executable_PROGRAM}"
     "--executable_pspamm=${PSpaMM_PROGRAM}"
     ${PREMULTIPLY_FLUX_VALUE} # boolean flag
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     DEPENDS
        build-time-make-directory
        generate.py
        config
        matrices
        kernels
    OUTPUT
        ${GENERATED_FILES_FOR_SEISSOL}
    COMMENT "Code generation for tensor operations..."
       )

add_custom_target(seissol-codegen ALL DEPENDS ${GENERATED_FILES_FOR_SEISSOL})

add_library(seissol-genkernel
    ${CMAKE_CURRENT_BINARY_DIR}/generated-code/kernel.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/generated-code/tensor.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/generated-code/subroutine.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/generated-code/init.cpp
)

target_include_directories(seissol-genkernel PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/generated-code
    ${CMAKE_SOURCE_DIR}/submodules/yateto/include
)

target_include_directories(seissol-genkernel INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}
)

if (HAS_REDZONE)
  set_source_files_properties(
      ${CMAKE_CURRENT_BINARY_DIR}/generated-code/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone
  )
endif()

add_dependencies(seissol-genkernel seissol-codegen)

# TODO: move here
# find_package(GemmTools REQUIRED)

if(GemmTools_INCLUDE_DIRS)
  target_include_directories(seissol-genkernel PRIVATE ${GemmTools_INCLUDE_DIRS})
endif()
if(GemmTools_LIBRARIES)
  target_link_libraries(seissol-genkernel PRIVATE ${GemmTools_LIBRARIES})
endif()
if(GemmTools_COMPILER_DEFINITIONS)
  target_compile_definitions(seissol-genkernel PRIVATE ${GemmTools_COMPILER_DEFINITIONS})
endif()

target_link_libraries(seissol-genkernel PRIVATE Eigen3::Eigen)
if (EIGEN3_INCLUDE_DIR)
  target_include_directories(seissol-genkernel PRIVATE ${EIGEN3_INCLUDE_DIR})
endif()

target_compile_options(seissol-genkernel PRIVATE -fPIC)
