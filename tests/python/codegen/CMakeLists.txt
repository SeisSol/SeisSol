# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Integrate the Python codegen tests into ctest. The `validation` and
# `ci` subdirectories are intentionally NOT wired here: they have no
# C++ build dependency and are better run via the dedicated
# pytest-python.yml workflow (and, for `ci`, pre-commit).
# Or just by calling `pytest` locally.

find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Probe whether pytest and hypothesis are importable. If not, skip the
# Python test integration with an informative message — don't hard-fail,
# so a minimal C++-only build stays usable.
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c
        "import pytest, hypothesis"
    RESULT_VARIABLE _pytest_probe
    OUTPUT_QUIET ERROR_QUIET
)

if (NOT _pytest_probe EQUAL 0)
    message(STATUS
        "Python test dependencies (pytest, hypothesis) not found; "
        "skipping Python codegen tests. Install them with "
        "'pip install pytest hypothesis' to enable.")
    return()
endif()

# The codegen tests import from codegen/ and its submodules/yateto
# submodule. conftest.py adds the right paths to sys.path based on
# REPO_ROOT = parents[3] (tests/python/codegen/conftest.py -> repo root).
# That works when the test files live in-tree at tests/python/codegen/;
# no extra environment setup needed.

# We invoke pytest once per test file rather than once for the whole
# subtree — this gives ctest per-file granularity, parallelism via
# `ctest -j`, and finer timeout control.
set(_pytest_codegen_files
    test_memlayout.py
    test_aderdg_math.py
    test_equations_structure.py
    test_integration.py
    test_generate_cli.py
    test_generate_goldenfile.py
)

foreach(_test_file IN LISTS _pytest_codegen_files)
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_test(
        NAME pytest-codegen-${_test_name}
        COMMAND ${Python3_EXECUTABLE} -m pytest
                ${CMAKE_CURRENT_SOURCE_DIR}/${_test_file}
                -v --tb=short --color=yes
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(pytest-codegen-${_test_name} PROPERTIES
        LABELS "python;codegen"
        # test_generate_goldenfile invokes generate.py four times and
        # is the slowest by ~30s; the rest are sub-second.
        TIMEOUT 300
    )
endforeach()

# The generate_goldenfile test shells out to generate.py, which imports
# yateto. If the build environment doesn't have yateto installed as a
# package, point PYTHONPATH at the in-tree submodule so the subprocess
# can find it — identical to the mechanism CMake itself uses to call
# generate.py at build time.
set_tests_properties(pytest-codegen-test_generate_goldenfile PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/submodules/yateto"
)
set_tests_properties(pytest-codegen-test_integration PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/submodules/yateto"
)
