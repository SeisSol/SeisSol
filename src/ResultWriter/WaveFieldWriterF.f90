!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
!!
!! @section LICENSE
!! Copyright (c) 2014, SeisSol Group
!! All rights reserved.
!! 
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!! 
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!! 
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!! 
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products derived from this
!!    software without specific prior written permission.
!! 
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!!
!! @section DESCRIPTION
!! Fortran interface to the wave field writer

module WaveFieldWriter
    use TypesDef

    use iso_c_binding

    implicit none

    ! C functions
    interface
        subroutine wavefield_hdf_init(rank, outputPrefix, dofs, numVars, order, &
                refinement, timestep) bind(C, name="wavefield_hdf_init")
            use, intrinsic :: iso_c_binding

            integer( kind=c_int ), value                       :: rank
            character( kind=c_char ), dimension(*), intent(in) :: outputPrefix
            real( kind=c_double ), dimension(*), intent(in)    :: dofs
            integer( kind=c_int ), value                       :: numVars
            integer( kind=c_int ), value                       :: order
            integer( kind=c_int ), value                       :: refinement
            integer( kind=c_int ), value                       :: timestep
        end subroutine wavefield_hdf_init

        subroutine waveFieldWriterClose() bind(C, name="wavefield_hdf_close")
            use, intrinsic :: iso_c_binding
            ! rename this to wavefield_hdf_close and create an new subroutine
            ! if we need to do additional convertions in Fortran
        end subroutine waveFieldWriterClose

        subroutine wavefield_hdf_write_step(time) bind(C, name="wavefield_hdf_write_step")
            use, intrinsic :: iso_c_binding

            real( kind=c_double ), value                    :: time
        end subroutine wavefield_hdf_write_step
    end interface

contains
    subroutine waveFieldWriterInit(timestep, disc, eqn, io, mesh, mpi)
        implicit none

        integer, intent(in)    :: timestep
        type (tDiscretization) :: disc
        type (tEquations)      :: eqn
        type (tInputOutput)    :: io
        type (tUnstructMesh)   :: mesh
        type (tMPI)            :: mpi

        call wavefield_hdf_init(mpi%myRank, trim(io%OutputFile) // c_null_char, &
            disc%galerkin%dgvar(:, :, :, 1), &
            eqn%nVarTotal, disc%spaceorder, &
            io%Refinement, timestep)
    end subroutine waveFieldWriterInit

    subroutine waveFieldWriterWriteStep(time, disc, mesh, mpi)
        implicit none

        real, intent(in)                   :: time
        type (tDiscretization), intent(in) :: disc
        type (tUnstructMesh), intent(in)   :: mesh
        type (tMPI), intent(in)            :: mpi

        call wavefield_hdf_write_step(time)
    end subroutine waveFieldWriterWriteStep

end module WaveFieldWriter
