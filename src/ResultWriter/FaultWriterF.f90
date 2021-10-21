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

#include "Initializer/preProcessorMacros.fpp"

module FaultWriter
    use TypesDef

    use iso_c_binding

    implicit none

    ! C functions
    interface
        subroutine fault_hdf_init(cells, vertices, nCells, nVertices, &
                outputMask, dataBuffer, outputPrefix, interval, xdmfWriterBackend) bind(C, name="fault_hdf_init")
            use, intrinsic :: iso_c_binding

            integer( kind=c_int ), dimension(*), intent(in)    :: cells
            real( kind=c_double ), dimension(*), intent(in)    :: vertices
            integer( kind=c_int ), value                       :: nCells
            integer( kind=c_int ), value                       :: nVertices
            integer( kind=c_int ), dimension(*), intent(in)    :: outputMask
            type(c_ptr), dimension(*), intent(in)              :: dataBuffer
            character( kind=c_char ), dimension(*), intent(in) :: outputPrefix
            real( kind=c_double ), value                       :: interval
            character(kind=c_char), dimension(*), intent(in)   :: xdmfWriterBackend
        end subroutine fault_hdf_init

        subroutine fault_hdf_close() bind(C, name="fault_hdf_close")
            use, intrinsic :: iso_c_binding
        end subroutine fault_hdf_close

        subroutine fault_hdf_write(time) bind(C, name="fault_hdf_write")
            use, intrinsic :: iso_c_binding

            real( kind=c_double ), value                    :: time
        end subroutine fault_hdf_write
    end interface

contains
    subroutine initFaultOutput(points, outputMask, dataBuffer, outputPrefix, interval, xdmfWriterBackend)
        implicit none

        type(tUnstructPoint), dimension(:)  :: points
        integer, dimension(:)               :: outputMask
        REAL_TYPE, dimension(:,:,:), target :: dataBuffer
        character(len=60)                   :: outputPrefix
        real                                :: interval
        character(len=64)                   :: xdmfWriterBackend

        integer, dimension(:,:), allocatable :: cells
        real, dimension(:,:), allocatable :: vertices
        integer :: nCells
        integer :: nVertices
        integer :: i, j
        REAL_TYPE, dimension(:), pointer :: dummyBuffer
        type(c_ptr) :: cDataBuffer(SIZE(dataBuffer, 3))

        ! TODO do not dublicate local vertices
        nCells = size(points)
        nVertices = 3*nCells

        allocate(cells(3, nCells), vertices(3, nVertices))

        do i=1,nCells
            cells(1, i) = (i-1)*3
            cells(2, i) = (i-1)*3 + 1
            cells(3, i) = (i-1)*3 + 2

            do j=1,3
                vertices(1, (i-1)*3+j) = points(i)%coordx(j)
                vertices(2, (i-1)*3+j) = points(i)%coordy(j)
                vertices(3, (i-1)*3+j) = points(i)%coordz(j)
            enddo
        enddo

        do i=1,size(dataBuffer, 3)
            dummyBuffer => dataBuffer(:,1,i)
            cDataBuffer(i) = c_loc(dummyBuffer(1))
        enddo

        call fault_hdf_init(cells, vertices, &
            nCells, nVertices, &
            outputMask, cDataBuffer, &
            trim(outputPrefix) // c_null_char, &
            interval, &
            trim(xdmfWriterBackend) // c_null_char)

        deallocate(cells, vertices)
    end subroutine initFaultOutput

    subroutine writeFault(time)
        implicit none

        real :: time

        call fault_hdf_write(time)
    end subroutine writeFault

    subroutine closeFaultOutput()
        implicit none

        call fault_hdf_close()
    end subroutine closeFaultOutput

end module FaultWriter
