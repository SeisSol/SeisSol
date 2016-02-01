!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2015, SeisSol Group
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
!! Velocity field reader Fortran interface

module VelocityFieldReader
    use TypesDef

    use iso_c_binding

    implicit none

    ! C function
    interface
        subroutine read_velocity_field(file, numElements, baryCenters, &
                defaultRho, defaultMu, defaultLambda, materialValues) bind(C, name="read_velocity_field")
            use, intrinsic :: iso_c_binding

            character( kind=c_char ), dimension(*), intent(in) :: file
            integer( kind=c_int ), value                       :: numElements
            real( kind=c_double ), dimension(*), intent(in)    :: baryCenters
            real( kind=c_double ), value                       :: defaultRho
            real( kind=c_double ), value                       :: defaultMu
            real( kind=c_double ), value                       :: defaultLambda
            real( kind=c_double ), dimension(*), intent(out)   :: materialValues

         end subroutine read_velocity_field
    end interface

contains
    subroutine readVelocityField(eqn, mesh, materialValues)
        implicit none

        type (tEquations), intent(in)    :: eqn
        type (tUnstructMesh), intent(in) :: mesh
        real, intent(inout)              :: materialValues(mesh%nElem, 3)

        call read_velocity_field(trim(eqn%MaterialFileName) // c_null_char, mesh%nElem, mesh%elem%xyBary, &
            eqn%rho0, eqn%mu, eqn%lambda, materialValues)
    end subroutine readVelocityField
end module
