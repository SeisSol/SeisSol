# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

function(determine_host_arch arch)
    include(CheckSourceRuns)
    include(CMakePushCheckState)

    set(myarch "noarch")

    cmake_push_check_state(RESET)

    set(CMAKE_REQUIRED_FLAGS "-march=native -mcpu=native -mtune=native")

    if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        check_source_runs(CXX
            "
            int main() {
            #ifdef __AVX__
                return 0;
            #else
                return 1;
            #endif
            }
            "
            HAS_AVX
        )
        check_source_runs(CXX
            "
            int main() {
            #ifdef __AVX2__
                return 0;
            #else
                return 1;
            #endif
            }
            "
            HAS_AVX2
        )
        check_source_runs(CXX
            "
            int main() {
            #ifdef __AVX512F__
                return 0;
            #else
                return 1;
            #endif
            }
            "
            HAS_AVX512
        )

        if (HAS_AVX512)
            set(myarch "avx10-512")
        elseif(HAS_AVX2)
            set(myarch "avx2-256")
        elseif(HAS_AVX)
            set(myarch "snb")
        endif()
    endif()
    if (CMAKE_SYSTEM_PROCESSOR STREQUAL "aarch64")
        set(myarch "neon")

        check_source_runs(CXX
            "
            #ifdef __ARM_FEATURE_SVE
            #include <arm_sve.h>
            int main() {
                return (svcntb() == 16) ? 0 : 1;
            }
            #else
            int main() {
                return 1;
            }
            #endif
            "
            HAS_SVE128
        )
        check_source_runs(CXX
            "
            #ifdef __ARM_FEATURE_SVE
            #include <arm_sve.h>
            int main() {
                return (svcntb() == 32) ? 0 : 1;
            }
            #else
            int main() {
                return 1;
            }
            #endif
            "
            HAS_SVE256
        )
        check_source_runs(CXX
            "
            #ifdef __ARM_FEATURE_SVE
            #include <arm_sve.h>
            int main() {
                return (svcntb() == 64) ? 0 : 1;
            }
            #else
            int main() {
                return 1;
            }
            #endif
            "
            HAS_SVE512
        )

        if (HAS_SVE512)
            set(myarch "sve-512")
        elseif(HAS_SVE256)
            set(myarch "sve-256")
        elseif(HAS_SVE128)
            set(myarch "sve-128")
        endif()
    endif()
    if (CMAKE_SYSTEM_PROCESSOR STREQUAL "riscv64")
        check_source_runs(CXX
            "
            #ifdef __riscv_vector
            #include <riscv_vector.h>
            int main() {
                return (__riscv_vsetvlmax_e8m1() == 16) ? 0 : 1;
            }
            #else
            int main() {
                return 1;
            }
            #endif
            "
            HAS_RVV128
        )
        check_source_runs(CXX
            "
            #ifdef __riscv_vector
            #include <riscv_vector.h>
            int main() {
                return (__riscv_vsetvlmax_e8m1() == 16) ? 0 : 1;
            }
            #else
            int main() {
                return 1;
            }
            #endif
            "
            HAS_RVV256
        )
        check_source_runs(CXX
            "
            #ifdef __riscv_vector
            #include <riscv_vector.h>
            int main() {
                return (__riscv_vsetvlmax_e8m1() == 16) ? 0 : 1;
            }
            #else
            int main() {
                return 1;
            }
            #endif
            "
            HAS_RVV512
        )

        if (HAS_RVV512)
            set(myarch "rvv-512")
        elseif(HAS_RVV256)
            set(myarch "rvv-256")
        elseif(HAS_RVV128)
            set(myarch "rvv-128")
        endif()
    endif()
    if (CMAKE_SYSTEM_PROCESSOR STREQUAL "loongarch64")
        # TODO: lsx lasx
    endif()

    cmake_pop_check_state()

    set(${arch} "${myarch}" PARENT_SCOPE)
endfunction()

function(determine_nvidia_arch arch)
    # TODO: just call CUDA?

    execute_process(COMMAND nvptx-arch OUTPUT_VARIABLE myarch)

    string(REGEX REPLACE "[ \t\r\n]" "" myarch "${myarch}")

    set(${arch} "${myarch}" PARENT_SCOPE)
endfunction()

function(determine_amd_arch arch)
    # TODO: just call HIP?

    execute_process(COMMAND amdgpu-arch OUTPUT_VARIABLE myarch)

    string(REGEX REPLACE "[ \t\r\n]" "" myarch "${myarch}")

    set(${arch} "${myarch}" PARENT_SCOPE)
endfunction()
