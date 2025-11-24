<!--
    SPDX-License-Identifier: BSD-3-Clause

    SPDX-FileCopyrightText: 2025 SeisSol Group
-->

UNI MPI (imported)
------------------

Uni MPI, imported from Petsc; licensed under BSD-2-Clause.
Slightly extended with missing types and functions;
and some defines that Petsc would normally provide.

It provides stubs for MPI library functions; thus we can build SeisSol
with MPI functions, but without an MPI library present.

Sources:

* https://gitlab.com/petsc/petsc/-/blob/main/src/sys/mpiuni/mpi.c
* https://gitlab.com/petsc/petsc/-/blob/main/src/sys/mpiuni/mpitime.c
* https://gitlab.com/petsc/petsc/-/blob/main/include/petsc/mpiuni/mpi.h
