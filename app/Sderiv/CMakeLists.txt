# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# sderiv, the command line tool of the programs SeisSol evaluates and of their evaluation records
add_executable(seissol-sderiv-exe Main.cpp)
target_link_libraries(seissol-sderiv-exe PRIVATE seissol-lib)
set_target_properties(seissol-sderiv-exe PROPERTIES OUTPUT_NAME "sderiv")
install(TARGETS seissol-sderiv-exe RUNTIME)

# the Python module of sderiv (import sderiv), on the same library
if (SDERIV_PYTHON)
    find_package(pybind11 2.12 REQUIRED)
    pybind11_add_module(seissol-sderiv-python Python.cpp)
    target_link_libraries(seissol-sderiv-python PRIVATE seissol-lib)
    set_target_properties(seissol-sderiv-python PROPERTIES OUTPUT_NAME "sderiv")
    if (TESTING)
        add_test(NAME "sderiv: Python module"
                 COMMAND ${Python3_EXECUTABLE} -m unittest discover
                         -s ${CMAKE_CURRENT_SOURCE_DIR}/python -p "test_*.py" -v)
        set_tests_properties("sderiv: Python module" PROPERTIES
            ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:seissol-sderiv-python>;PYTHONDONTWRITEBYTECODE=1;SDERIV_EXAMPLE=${CMAKE_CURRENT_SOURCE_DIR}/Example.sderiv")
    endif()
endif()

if (TESTING)
    seissol_add_smoke_test(
        NAME "sderiv: no command"
        TARGET seissol-sderiv-exe
        EXPECT_FAILURE
        EXPECT_OUTPUT "usage: sderiv"
    )
    seissol_add_smoke_test(
        NAME "sderiv: print"
        TARGET seissol-sderiv-exe
        ARGS print ${CMAKE_CURRENT_SOURCE_DIR}/Example.sderiv
        EXPECT_OUTPUT "out rho = "
    )
    seissol_add_smoke_test(
        NAME "sderiv: eval"
        TARGET seissol-sderiv-exe
        ARGS eval ${CMAKE_CURRENT_SOURCE_DIR}/Example.sderiv z=-3000
        EXPECT_OUTPUT "rho = 2700"
    )
endif()
