!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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

#define mayDeallocateA(value) if (allocated(value)) deallocate(value)
#define mayDeallocateP(value) if (associated(value)) deallocate(value)

module GambitReaderTestTools
    use TypesDef
    use read_mesh_mod
    use MPIExtractMesh_mod
    use MeshReaderCBinding

    use iso_c_binding

    implicit none

    type (tUnstructDomainDescript) :: domain

contains
    subroutine ReadCubeOld(rank) bind(C)
        implicit none

        integer, intent(in), value :: rank

        ! TODO Use SEISSOL_TESTS
        domain%IO%MeshFile = 'src/tests/Geometry/cube4.neu'
        domain%IO%MetisFile = 'src/tests/Geometry/cube4.met'

        domain%BND%NoBndObjects(5) = 1

        call ReadMeshOld(rank)
    end subroutine ReadCubeOld

    subroutine ReadDRBoxOld(rank) bind(C)
        implicit none

        integer, intent(in), value :: rank

        domain%IO%MeshFile = 'src/tests/Geometry/scec.neu'
        domain%IO%MetisFile = 'src/tests/Geometry/scec.met'

        domain%BND%NoBndObjects(1) = 1
        domain%BND%NoBndObjects(3) = 1
        domain%BND%NoBndObjects(5) = 1

        domain%EQN%dr = 1
        domain%EQN%XRef = 1.0e5
        domain%EQN%YRef = -1.0e5
        domain%EQN%ZRef = -10.0e5

        call ReadMeshOld(rank)
    end subroutine ReadDRBoxOld

    !> Reads a test mesh using the old fortran routines.
    subroutine ReadMeshOld(rank)
        implicit none

        integer, intent(in), value :: rank

        ! Free allocated arrays in domain
        mayDeallocateA(domain%MESH%LocalElemType)
        mayDeallocateA(domain%MESH%LocalVrtxType)
        mayDeallocateA(domain%MESH%LocalElemIndex_Tet)
        mayDeallocateP(domain%MESH%VRTX%xyNode)
        mayDeallocateP(domain%MESH%VRTX%NrOfElementsConnected)
        mayDeallocateP(domain%MESH%VRTX%Element)
        mayDeallocateP(domain%MESH%ELEM%Reference)
        mayDeallocateP(domain%MESH%ELEM%SideNeighbor)
        mayDeallocateP(domain%MESH%ELEM%LocalNeighborSide)
        mayDeallocateP(domain%MESH%ELEM%LocalNeighborVrtx)

        ! Set variables required to read the mesh
        domain%IO%meshgenerator = 'Gambit3D-Mixed'
        domain%IO%UNIT%FileIn = 500

        domain%EQN%Dimension = 3

        allocate(domain%MESH%Displacement(3))
        domain%MESH%Displacement(:) = 0
        allocate(domain%MESH%ScalingMatrix(3,3))
        domain%MESH%ScalingMatrix(:,:) = 0
        domain%MESH%ScalingMatrix(1,1) = 1
        domain%MESH%ScalingMatrix(2,2) = 1
        domain%MESH%ScalingMatrix(3,3) = 1
        domain%MESH%Dimension = 3

        myrank = rank
        domain%MPI%myRank = rank
        domain%MPI%nCPU = 2

        call read_mesh(domain%IO, domain%EQN, domain%DISC, domain%MESH, domain%BND, domain%MPI)
        call MPIExtractMesh(domain%EQN, domain%DISC, domain%BND, domain%MESH, domain%IO, domain%MPI)

        ! Free resources
        deallocate(domain%MESH%Displacement)
        deallocate(domain%MESH%ScalingMatrix)

        call cBindingStructs(domain%EQN, domain%DISC, domain%MESH, domain%BND, domain%MPI)
    end subroutine ReadMeshOld

end module GambitReaderTestTools
