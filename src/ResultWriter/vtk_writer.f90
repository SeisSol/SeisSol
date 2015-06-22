!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Atanas Atanasov (atanasoa AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Atanas_Atanasov)
!!
!! @section LICENSE
!! Copyright (c) 2013, SeisSol Group
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
!! the vtk plotter class encapsolates the functionality for plotting unstructured grids via vtk file

module vtk
implicit none

interface
    subroutine create_vtk_writer(reference,rank,iteration,basename, binaryoutput)
      use, intrinsic :: iso_c_binding
      integer(8) :: reference
      integer(4) :: rank
      integer(4) :: iteration
      integer(4) :: binaryoutput
      character( kind=c_char ), dimension(*), intent(in) :: basename
    end subroutine create_vtk_writer

	subroutine destroy_vtk_writer(reference)
      integer(8) :: reference
    end subroutine destroy_vtk_writer

    subroutine insert_vertex_vtk_writer(reference,x,y,z)
      integer(8) :: reference
      real(8)::x
      real(8)::y
      real(8)::z
    end subroutine insert_vertex_vtk_writer

    subroutine write_vertices_vtk_writer(reference)
      integer(8) :: reference
    end subroutine write_vertices_vtk_writer

	subroutine start_cell_data_vtk_writer(reference,var_id)
	  integer(8) :: reference
	  integer(4) :: var_id
	end subroutine start_cell_data_vtk_writer

	subroutine end_cell_data_vtk_writer(reference)
	  integer(8) :: reference
	end subroutine end_cell_data_vtk_writer

	subroutine plot_cell_data_vtk_writer(reference,value)
	  integer(8) :: reference
	  real(8) :: value

	end subroutine plot_cell_data_vtk_writer

	subroutine plot_cells_vtk_writer(reference)
	  integer(8) :: reference
	end subroutine plot_cells_vtk_writer

    subroutine open_vtk_writer(reference)
      integer(8) :: reference
    end subroutine open_vtk_writer

    subroutine close_vtk_writer(reference)
      integer(8) :: reference
    end subroutine close_vtk_writer
end interface

end module vtk
