# SPDX-License-Identifier: BSD-3-Clause
#
# SPDX-FileCopyrightText: 2025 SeisSol Group

project(unimpi LANGUAGES C)

add_library(unimpi src/mpi.c src/mpitime.c)
target_include_directories(unimpi PUBLIC include)
