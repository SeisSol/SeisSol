# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The runtime options of the sanitizers (cmake/sanitizers.cmake) that the CI runs the tests with.
# In bash, before running the tests of a sanitized build:
#
#   source .ci/sanitizers/environment.sh
#
# Options set before take precedence over the ones here. For the thread sanitizer, it also
# loads Archer, through which the OpenMP runtime of LLVM tells the sanitizer how it synchronizes
# its threads; it looks for it next to the compiler in CXX (default: clang++).

seissol_sanitizer_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export ASAN_OPTIONS="detect_leaks=1:detect_stack_use_after_return=1:check_initialization_order=1:\
strict_init_order=1:${ASAN_OPTIONS:-}"
export LSAN_OPTIONS="suppressions=${seissol_sanitizer_dir}/lsan.supp:print_suppressions=0:\
${LSAN_OPTIONS:-}"
export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1:\
suppressions=${seissol_sanitizer_dir}/ubsan.supp:print_suppressions=0:${UBSAN_OPTIONS:-}"
# The MPI library and the OpenMP runtime are not instrumented: the sanitizer ignores what they
# access, which it could not judge without seeing how they synchronize.
export TSAN_OPTIONS="ignore_noninstrumented_modules=1:history_size=4:second_deadlock_stack=1:\
suppressions=${seissol_sanitizer_dir}/tsan.supp:print_suppressions=0:${TSAN_OPTIONS:-}"
# With the thread sanitizer, PMIx (of Open MPI) cannot attach its shared memory at the address it
# wants, which the sanitizer keeps for itself, and MPI_Init crashes; its hash store needs none.
export PMIX_MCA_gds="${PMIX_MCA_gds:-hash}"

if [ -z "${OMP_TOOL_LIBRARIES:-}" ]; then
  seissol_sanitizer_cxx="$(command -v "${CXX:-clang++}" 2>/dev/null || true)"
  if [ -n "${seissol_sanitizer_cxx}" ]; then
    seissol_sanitizer_archer="$(dirname "$(readlink -f "${seissol_sanitizer_cxx}")")/../lib/libarcher.so"
    if [ -f "${seissol_sanitizer_archer}" ]; then
      OMP_TOOL_LIBRARIES="$(readlink -f "${seissol_sanitizer_archer}")"
      export OMP_TOOL_LIBRARIES
    fi
  fi
  unset seissol_sanitizer_cxx seissol_sanitizer_archer
fi

unset seissol_sanitizer_dir
