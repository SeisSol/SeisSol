!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2016-2017, SeisSol Group
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

module StressReader
	use TypesDef

	use iso_c_binding

	implicit none

	! C function
	interface
		subroutine open_stress_field(file, friction) bind (C, name="open_stress_field")
			use, intrinsic :: iso_c_binding

			character( kind=c_char ), dimension(*), intent(in) :: file
			integer( kind=c_int ), value                       :: friction
		end subroutine open_stress_field

		subroutine readStress(x, y, z, stressValues, frictionValues, stressDefaultValues, frictionDefaultValues) bind(C, name="read_stress")
			use, intrinsic :: iso_c_binding

			real( kind=c_double ), value :: x
			real( kind=c_double ), value :: y
			real( kind=c_double ), value :: z
			real( kind=c_double ), dimension(*), intent(inout) :: stressValues
			real( kind=c_double ), dimension(*), intent(inout) :: frictionValues
			real( kind=c_float ), dimension(*), intent(in)  :: stressDefaultValues
			real( kind=c_float ), dimension(*), intent(in)  :: frictionDefaultValues
		end subroutine readStress

		subroutine closeStressField() bind(C, name="close_stress_field")
			use, intrinsic :: iso_c_binding
		end subroutine closeStressField
	end interface

	contains
	subroutine openStressField(io, friction)
		use, intrinsic :: iso_c_binding
		implicit none

		type(tInputOutput)             :: io
		integer                        :: friction

		call open_stress_field(trim(io%FileName_BackgroundStress) // c_null_char, friction)
	end subroutine openStressField
end module