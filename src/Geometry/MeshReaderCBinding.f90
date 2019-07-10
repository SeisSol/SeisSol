!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2013-2016, SeisSol Group
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
!! Provides access to the important Fortran arrays in C/C++ and implements read_mesh_fast for SeisSol

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

module MeshReaderCBinding
    use TypesDef
    use QuadPoints_mod

    use iso_c_binding

    implicit none

    type (tEquations), pointer :: m_eqn
    type (tDiscretization), pointer :: m_disc
    type (tUnstructMesh), pointer :: m_mesh
    type (tBoundary), pointer :: m_bnd
    type (tMPI), pointer :: m_mpi

    interface
        subroutine read_mesh_gambitfast_c(rank, meshfile, partitionfile, hasFault, displacement, scalingMatrix) bind(C, name="read_mesh_gambitfast_c")
            use, intrinsic :: iso_c_binding

            integer( kind=c_int ), value                       :: rank
            character( kind=c_char ), dimension(*), intent(in) :: meshfile
            character( kind=c_char ), dimension(*), intent(in) :: partitionfile
            logical( kind=c_bool ), value                      :: hasFault
            real(kind=c_double), dimension(*), intent(in)      :: displacement
            real(kind=c_double), dimension(*), intent(in)      :: scalingMatrix
        end subroutine

        subroutine read_mesh_netcdf_c(rank, nprocs, meshfile, hasFault, displacement, scalingMatrix) bind(C, name="read_mesh_netcdf_c")
            use, intrinsic :: iso_c_binding

            integer( kind=c_int ), value                       :: rank
            integer( kind=c_int ), value                       :: nprocs
            character( kind=c_char ), dimension(*), intent(in) :: meshfile
            logical( kind=c_bool ), value                      :: hasFault
            real(kind=c_double), dimension(*), intent(in)      :: displacement
            real(kind=c_double), dimension(*), intent(in)      :: scalingMatrix
        end subroutine

        subroutine read_mesh_puml_c(meshfile, checkPointFile, hasFault, displacement, scalingMatrix, easiVelocityModel, clusterRate) bind(C, name="read_mesh_puml_c")
            use, intrinsic :: iso_c_binding

            character( kind=c_char ), dimension(*), intent(in) :: meshfile, easiVelocityModel, checkPointFile
            logical( kind=c_bool ), value                      :: hasFault
            real(kind=c_double), dimension(*), intent(in)      :: displacement
            real(kind=c_double), dimension(*), intent(in)      :: scalingMatrix
            integer(kind=c_int), value, intent(in)                :: clusterRate
        end subroutine
    end interface

contains
    !> Calls the C++ function
    subroutine read_mesh_fast(io, eqn, disc, mesh, bnd, mpi)
        implicit none

        type (tInputOutput) :: io
        type (tEquations) :: eqn
        type (tDiscretization) :: disc
        type (tUnstructMesh) :: mesh
        type (tBoundary) :: bnd
        type (tMPI) :: mpi

        character*50 str

        integer i
        integer nVertices
        integer nElements
        integer nBndMPI
        logical( kind=c_bool ) hasFault

        call cBindingStructs(eqn, disc, mesh, bnd, mpi)

        if (eqn%DR .eq. 1) then
            hasFault = .true.
        else
            hasFault = .false.
        endif
        disc%DynRup%DR_output = .false.

        write(str, *) mpi%nCPU
        if (io%meshgenerator .eq. 'Gambit3D-fast') then
            call read_mesh_gambitfast_c(mpi%myRank, trim(io%MeshFile) // c_null_char, \
                trim(io%MetisFile) // '.epart.' // trim(adjustl(str)) // c_null_char, hasFault, MESH%Displacement(:), m_mesh%ScalingMatrix(:,:))
        elseif (io%meshgenerator .eq. 'Netcdf') then
#ifdef PARALLEL
            call read_mesh_netcdf_c(mpi%myRank, mpi%nCPU, trim(io%MeshFile) // c_null_char, hasFault, MESH%Displacement(:), m_mesh%ScalingMatrix(:,:))
#else
            call read_mesh_netcdf_c(0, 1, trim(io%MeshFile) // c_null_char, hasFault, MESH%Displacement(:), m_mesh%ScalingMatrix(:,:))
#endif
        elseif (io%meshgenerator .eq. 'PUML') then
            call read_mesh_puml_c(  trim(io%MeshFile) // c_null_char,           &
                                    trim(io%checkpoint%filename) // c_null_char,&
                                    hasFault,                                   &
                                    MESH%Displacement(:),                       &
                                    m_mesh%ScalingMatrix(:,:),                  &
                                    trim(EQN%MaterialFileName) // c_null_char,  &
                                    disc%galerkin%clusteredLts                  )
        else
            logError(*) 'Unknown mesh reader'
            stop
        endif

        ! Set additional SeisSol variables
        nVertices = size(mesh%VRTX%xyNode, 2)
        nElements = size(mesh%ELEM%Vertex, 2)
        nBndMPI = size(bnd%ObjMPI)

        mesh%nNode = nVertices
        mesh%nElem = nElements
        mesh%nElem_Tet = nElements
        mesh%nElem_Hex = 0

        bnd%NoMPIDomains = nBndMPI

        if (eqn%DR .eq. 1) then
            mesh%Fault%nSide = size(mesh%Fault%Face, 1)
            if (mesh%Fault%nSide .eq. 0) then
                eqn%DR = 0
            endif
        endif

#ifdef USE_DR_CELLAVERAGE
        DISC%Galerkin%nBndGP = 4**ceiling(log( real((DISC%Galerkin%nPoly + 1)*(DISC%Galerkin%nPoly + 2) / 2) )/log(4.))
#else
        DISC%Galerkin%nBndGP = (DISC%Galerkin%nPoly + 2)**2
#endif
        disc%Galerkin%nIntGP = (disc%Galerkin%nPoly+2)**3

        allocate(mesh%LocalVrtxType(nElements))
        mesh%LocalVrtxType(:) = 4
        allocate(mesh%LocalElemType(nElements))
        mesh%LocalElemType(:) = 4

        allocate(mesh%ELEM%BndGP_Tri(2, DISC%Galerkin%nBndGP), &
                 mesh%ELEM%BndGW_Tri(DISC%Galerkin%nBndGP))
#ifdef USE_DR_CELLAVERAGE
        call CellCentresOfSubdivision(DISC%Galerkin%nPoly + 1, mesh%ELEM%BndGP_Tri)
        mesh%ELEM%BndGW_Tri = 1.e99 ! blow up solution if used
#else
        call TriangleQuadraturePoints(                    &
                  nIntGP     = disc%Galerkin%nBndGP,      &
                  IntGaussP  = mesh%ELEM%BndGP_Tri,       &
                  IntGaussW  = mesh%ELEM%BndGW_Tri,       &
                  M          = disc%Galerkin%nPoly+2,     &
                  IO         = io,                        &
                  quiet      = .true.                     )
#endif

        call computeAdditionalMeshInfo()

        allocate(mesh%ELEM%MinDistBarySide(nElements))

        call exchangeBndCellInfo()
    end subroutine read_mesh_fast

    !> Compute additional mesh information (cell local information)
    subroutine computeAdditionalMeshInfo()
        implicit none

        real :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
        integer :: i, nElements, nVertices

        nElements = size(m_mesh%ELEM%Vertex, 2)

        allocate(m_mesh%ELEM%Volume(nElements))
        allocate(m_mesh%ELEM%xyBary(3, nElements))

        do i=1,nElements
            x1 = m_mesh%VRTX%xyNode(1,m_mesh%ELEM%Vertex(1,i))
            y1 = m_mesh%VRTX%xyNode(2,m_mesh%ELEM%Vertex(1,i))
            z1 = m_mesh%VRTX%xyNode(3,m_mesh%ELEM%Vertex(1,i))
            !
            x2 = m_mesh%VRTX%xyNode(1,m_mesh%ELEM%Vertex(2,i))
            y2 = m_mesh%VRTX%xyNode(2,m_mesh%ELEM%Vertex(2,i))
            z2 = m_mesh%VRTX%xyNode(3,m_mesh%ELEM%Vertex(2,i))
            !
            x3 = m_mesh%VRTX%xyNode(1,m_mesh%ELEM%Vertex(3,i))
            y3 = m_mesh%VRTX%xyNode(2,m_mesh%ELEM%Vertex(3,i))
            z3 = m_mesh%VRTX%xyNode(3,m_mesh%ELEM%Vertex(3,i))
            !
            x4 = m_mesh%VRTX%xyNode(1,m_mesh%ELEM%Vertex(4,i))
            y4 = m_mesh%VRTX%xyNode(2,m_mesh%ELEM%Vertex(4,i))
            z4 = m_mesh%VRTX%xyNode(3,m_mesh%ELEM%Vertex(4,i))
            !
            ! Compute with SARRUS's law from the Jacobi determinant of the transformation
            ! from xyz coordinates to the reference tetrahedron
            !
            m_mesh%ELEM%Volume(i) = abs((x2-x1)*(y3-y1)*(z4-z1) +             &
                                        (x3-x1)*(y4-y1)*(z2-z1) +             &
                                        (x4-x1)*(y2-y1)*(z3-z1) -             &
                                        (x4-x1)*(y3-y1)*(z2-z1) -             &
                                        (x2-x1)*(y4-y1)*(z3-z1) -             &
                                        (x3-x1)*(y2-y1)*(z4-z1)         )/6.0D0
            ! Barycenters of the tetrahedrons
            m_mesh%ELEM%xyBary(:,i) = 0.25*(m_mesh%VRTX%xyNode(:,m_mesh%ELEM%Vertex(1,i)) + &
                                            m_mesh%VRTX%xyNode(:,m_mesh%ELEM%Vertex(2,i)) + &
                                            m_mesh%VRTX%xyNode(:,m_mesh%ELEM%Vertex(3,i)) + &
                                            m_mesh%VRTX%xyNode(:,m_mesh%ELEM%Vertex(4,i))   )
        enddo
    end subroutine computeAdditionalMeshInfo

    subroutine exchangeBndCellInfo()
        implicit none
#ifdef PARALLEL
        include 'mpif.h'

        integer :: i, j, k, iElem
        integer :: iErr
        type(tRealMessage), allocatable :: send_messages(:), recv_messages(:)
        integer, allocatable            :: send_requests(:), recv_requests(:)
        integer, allocatable            :: send_status_list(:,:), recv_status_list(:,:)

        allocate(send_messages(m_bnd%NoMPIDomains), &
                 recv_messages(m_bnd%NoMPIDomains))
        allocate(send_requests(m_bnd%NoMPIDomains), &
                 recv_requests(m_bnd%NoMPIDomains))
        allocate(send_status_list(MPI_STATUS_SIZE,m_bnd%NoMPIDomains), &
                 recv_status_list(MPI_STATUS_SIZE,m_bnd%NoMPIDomains))

        ! Recv/Send local vertices
        do i=1,m_bnd%NoMPIDomains
            allocate(recv_messages(i)%content(m_bnd%ObjMPI(i)%nElem*4*3))
            call mpi_irecv(recv_messages(i)%content, m_bnd%ObjMPI(i)%nElem*4*3, m_mpi%MPI_AUTO_REAL, m_bnd%ObjMPI(i)%CPU, 1, m_mpi%commWorld, recv_requests(i), iErr)

            allocate(send_messages(i)%content(m_bnd%ObjMPI(i)%nElem*4*3))

            do j=1,m_bnd%ObjMPI(i)%nElem
                iElem = m_bnd%ObjMPI(i)%DomainElements(j)
                do k=1,4
                    send_messages(i)%content(((j-1)*4+k-1)*3+1:((j-1)*4+k-1)*3+3) = m_mesh%VRTX%xyNode(1:3,m_mesh%ELEM%Vertex(k,iElem))
                enddo
            enddo
            call mpi_isend(send_messages(i)%content, m_bnd%ObjMPI(i)%nElem*4*3, m_mpi%MPI_AUTO_REAL, m_bnd%ObjMPI(i)%CPU, 1, m_mpi%commWorld, send_requests(i), iErr)
        enddo

        ! exchange data
        call mpi_waitall(m_bnd%NoMPIDomains, send_requests, send_status_list, iErr)
        call mpi_waitall(m_bnd%NoMPIDomains, recv_requests, recv_status_list, iErr)

        ! Copy received data/deallocate buffers
        do i=1,m_bnd%NoMPIDomains
            deallocate(send_messages(i)%content)

            allocate(m_bnd%ObjMPI(i)%NeighborCoords(m_eqn%dimension, 4, m_bnd%ObjMPI(i)%nElem))

            do j=1,m_bnd%ObjMPI(i)%nElem
                do k=1,4
                    m_bnd%ObjMPI(i)%NeighborCoords(1:3, k, j) = recv_messages(i)%content(((j-1)*4+k-1)*3+1:((j-1)*4+k-1)*3+3)
                enddo
            enddo

            deallocate(recv_messages(i)%content)
        enddo

        deallocate(send_messages, &
                 recv_messages)
        deallocate(send_requests, &
                 recv_requests)
        deallocate(send_status_list, &
                 recv_status_list)
#endif
    end subroutine exchangeBndCellInfo

    !> Sets the pointer to the global data structure storing all SeisSol
    !! information
    subroutine cBindingStructs(eqn, disc, mesh, bnd, mpi)
        implicit none

        type (tEquations), target :: eqn
        type (tDiscretization), target :: disc
        type (tUnstructMesh), target :: mesh
        type (tBoundary), target :: bnd
        type (tMPI), target :: mpi

        m_eqn => eqn
        m_disc => disc
        m_mesh => mesh
        m_bnd => bnd
        m_mpi => mpi
    end subroutine cBindingStructs

    subroutine allocVertices(n, maxElements) bind(C)
        implicit none

        integer( kind=c_int ), value :: n
        integer( kind=c_int ), value :: maxElements

        allocate(m_mesh%VRTX%xyNode(3, n))
        allocate(m_mesh%VRTX%NrOfElementsConnected(n))
        allocate(m_mesh%VRTX%Element(n, maxElements)) ! TODO swap rows/colomns

        ! Not required in our implementation, but will be freed by SeisSol
        allocate(m_mesh%VRTX%Reference(0))
        allocate(m_mesh%VRTX%BoundaryToObject(0))
    end subroutine allocVertices

    subroutine allocElements(n) bind(C)
        implicit none

        integer( kind=c_int ), value :: n

        allocate(m_mesh%ELEM%Vertex(4, n))
        allocate(m_mesh%ELEM%Reference(0:4, n))
        allocate(m_mesh%ELEM%MPIReference(0:4, n))
        allocate(m_mesh%ELEM%MPINumber(4, n))
        if (m_eqn%DR .eq. 1) then
            allocate(m_mesh%ELEM%MPINumber_DR(4, n))
        endif
        allocate(m_mesh%ELEM%BoundaryToObject(4, n))
        allocate(m_mesh%ELEM%SideNeighbor(4, n))
        allocate(m_mesh%ELEM%LocalNeighborSide(4, n))
        allocate(m_mesh%ELEM%LocalNeighborVrtx(4, n))
    end subroutine allocElements

    subroutine allocBndObjs(n) bind (C)
        implicit none

        integer( kind=c_int ), value :: n

        allocate(m_bnd%ObjMPI(n))
    end subroutine allocBndObjs

    subroutine allocBndObj(i, n) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), value :: n

        allocate(m_bnd%ObjMPI(i)%DomainElements(n))
    end subroutine allocBndObj

    subroutine allocFault(n) bind(C)
        implicit none

        integer( kind=c_int ), value :: n

        allocate(m_mesh%Fault%Face(n, 2, 2))
        allocate(m_mesh%Fault%geoNormals(3, n))
        allocate(m_mesh%Fault%geoTangent1(3, n))
        allocate(m_mesh%Fault%geoTangent2(3, n))
    end subroutine allocFault

    subroutine hasPlusFault() bind(C)
        implicit none

        m_disc%DynRup%DR_output = .true.
    end subroutine hasPlusFault

    subroutine allocBndObjFault(i, n) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), value :: n

        allocate(m_bnd%ObjMPI(i)%Domain_Fault_Elem(n))
    end subroutine allocBndObjFault

    subroutine getVerticesXY(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%VRTX%xyNode, 2)
        ptr = c_loc(m_mesh%VRTX%xyNode(1,1))
    end subroutine getVerticesXY

    subroutine getVerticesNElements(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%VRTX%NrOfElementsConnected)
        ptr = c_loc(m_mesh%VRTX%NrOfElementsConnected(1))
    end subroutine getVerticesNElements

    subroutine getVerticesElements(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%VRTX%Element, 1)
        ptr = c_loc(m_mesh%VRTX%Element(1,1))
    end subroutine getVerticesElements

    subroutine getElementVertices(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%Vertex, 2)
        ptr = c_loc(m_mesh%ELEM%Vertex(1,1))
    end subroutine getElementVertices

    subroutine getReference(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%Reference, 2)
        ptr = c_loc(m_mesh%ELEM%Reference(0,1))
    end subroutine getReference

    subroutine getMPIReference(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%MPIReference, 2)
        ptr = c_loc(m_mesh%ELEM%MPIReference(0,1))
    end subroutine getMPIReference

    subroutine getMPINumber(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%MPINumber, 2)
        ptr = c_loc(m_mesh%ELEM%MPINumber(1,1))
    end subroutine getMPINumber

    subroutine getMPINumberDR(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%MPINumber_DR, 2)
        ptr = c_loc(m_mesh%ELEM%MPINumber_DR(1,1))
    end subroutine getMPINumberDR

    subroutine getBoundaryToObject(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%BoundaryToObject, 2)
        ptr = c_loc(m_mesh%ELEM%BoundaryToObject(1, 1))
    end subroutine getBoundaryToObject

    subroutine getSideNeighbor(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%SideNeighbor, 2)
        ptr = c_loc(m_mesh%ELEM%SideNeighbor(1,1))
    end subroutine getSideNeighbor

    subroutine getLocalNeighborSide(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%LocalNeighborSide, 2)
        ptr = c_loc(m_mesh%ELEM%LocalNeighborSide(1, 1))
    end subroutine getLocalNeighborSide

    subroutine getLocalNeighborVrtx(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%ELEM%LocalNeighborVrtx, 2)
        ptr = c_loc(m_mesh%ELEM%LocalNeighborVrtx(1, 1))
    end subroutine getLocalNeighborVrtx

    subroutine getBndSize(s) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s

        s = size(m_bnd%ObjMPI)
    end subroutine getBndSize

    subroutine getBndRank(i, rank) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), intent(out) :: rank

        rank = m_bnd%ObjMPI(i)%CPU
    end subroutine getBndRank

    subroutine setBndRank(i, rank) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), value :: rank

        m_bnd%ObjMPI(i)%CPU = rank
    end subroutine setBndRank

    subroutine getBndNElem(i, n) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), intent(out) :: n

        n = m_bnd%ObjMPI(i)%nElem
    end subroutine getBndNElem

    subroutine setBndNElem(i, n) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), value :: n

        m_bnd%ObjMPI(i)%nElem = n
    end subroutine setBndNElem

    subroutine getBndDomainElements(i, s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_bnd%ObjMPI(i)%DomainElements)
        ptr = c_loc(m_bnd%ObjMPI(i)%DomainElements(1))
    end subroutine getBndDomainElements

    subroutine getFaultReferencePoint(x, y, z, method) bind(C)
        implicit none

        real( kind=c_double ), intent(out) :: x
        real( kind=c_double ), intent(out) :: y
        real( kind=c_double ), intent(out) :: z
        integer( kind=c_int ), intent(out) :: method

        x = m_eqn%XRef
        y = m_eqn%YRef
        z = m_eqn%ZRef
        method = m_eqn%refPointMethod
    end subroutine getFaultReferencePoint

    subroutine getFaultFace(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%Fault%Face, 1)
        ptr = c_loc(m_mesh%Fault%Face(1,1,1))
    end subroutine getFaultFace

    subroutine getFaultNormals(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%Fault%geoNormals, 2)
        ptr = c_loc(m_mesh%Fault%geoNormals(1,1))
    end subroutine getFaultNormals

    subroutine getFaultTangent1(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%Fault%geoTangent1, 2)
        ptr = c_loc(m_mesh%Fault%geoTangent1(1,1))
    end subroutine getFaultTangent1

    subroutine getFaultTangent2(s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_mesh%Fault%geoTangent2, 2)
        ptr = c_loc(m_mesh%Fault%geoTangent2(1,1))
    end subroutine getFaultTangent2

    subroutine setBndFaultNElem(i, n) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), value :: n

        m_bnd%ObjMPI(i)%nFault_MPI = n
    end subroutine setBndFaultNElem

    subroutine getBndFaultElements(i, s, ptr) bind(C)
        implicit none

        integer( kind=c_int ), value :: i
        integer( kind=c_int ), intent(out) :: s
        type(c_ptr), intent(out) :: ptr

        s = size(m_bnd%ObjMPI(i)%Domain_Fault_Elem)
        ptr = c_loc(m_bnd%ObjMPI(i)%Domain_Fault_Elem(1))
    end subroutine getBndFaultElements

end module MeshReaderCBinding
