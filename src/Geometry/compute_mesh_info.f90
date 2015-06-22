!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) 2010-2014, SeisSol Group
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
!! Compute the information of mesh, e.g. node and element
!! Also compute the neighbour information

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE compute_mesh_info_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  USE common_operators_mod
  USE QuadPoints_mod
  USE DGBasis_mod
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  INTERFACE compute_mesh_Gambit3D_Tetra
     MODULE PROCEDURE compute_mesh_Gambit3D_Tetra
  END INTERFACE

  INTERFACE compute_mesh_ICEMCFD3D_Tetra
     MODULE PROCEDURE compute_mesh_ICEMCFD3D_Tetra
  END INTERFACE

  INTERFACE compute_mesh_Gambit3D_Hexa
     MODULE PROCEDURE compute_mesh_Gambit3D_Hexa
  END INTERFACE

  INTERFACE compute_mesh_Gambit3D_Mixed
     MODULE PROCEDURE compute_mesh_Gambit3D_Mixed
  END INTERFACE

  INTERFACE BND_DG_Periodic3D_Tetra_us
     MODULE PROCEDURE BND_DG_Periodic3D_Tetra_us
  END INTERFACE

  INTERFACE BND_DG_Periodic3D_Hexa_us
     MODULE PROCEDURE BND_DG_Periodic3D_Hexa_us
  END INTERFACE

  !------------------------------------------
  PUBLIC	:: compute_mesh_Gambit3D_Tetra
  PUBLIC	:: compute_mesh_ICEMCFD3D_Tetra
  PUBLIC	:: compute_mesh_Gambit3D_Hexa
  PUBLIC	:: compute_mesh_Gambit3D_Mixed
  PUBLIC	:: BND_DG_Periodic3D_Tetra_us
  PUBLIC	:: BND_DG_Periodic3D_Hexa_us

CONTAINS

!> Compute the mesh and neighbour information of Gambit3D tetrahedrons mesh
!<
  SUBROUTINE compute_mesh_Gambit3D_Tetra(Mesh, IO, EQN, BND)
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH                                          !< Mesh data
    TYPE (tInputOutput)       :: IO                                            !< IO data
    TYPE (tEquations)         :: EQN                                           !< Discretization data
    TYPE (tBoundary)		  :: BND                                           !< Boundary data
    !--------------------------------------------------------------------------
    !local variables
    !--------------------------------------------------------------------------
    INTEGER 				  :: iNode, iVrtx, iElem, iSide, iNeighbor
    INTEGER, POINTER          :: test_array1(:,:)
    INTEGER, POINTER          :: test_array2(:,:)
    REAL                      :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
    REAL                      :: temp(EQN%Dimension)
    INTEGER                   :: test_array_size
    INTEGER,DIMENSION(3)      :: v_tmp,n_tmp
    INTEGER                   :: iRow,iColumns(1),iPosition(1)
    INTEGER					  :: i,j,k,testI, counter
    INTEGER                   :: orientation(4,4)
    !--------------------------------------------------------------------------
    INTENT(INOUT)			  :: Mesh
    INTENT(IN) 				  :: EQN,IO,BND

	    !
    ! Assuming, that faces consist of the three vertices as follow
    !
    !                              Face       Vertices
    !                4              1          1,3,2
    !                *              2          1,2,4
    !               /|\             3          1,4,3
    !             /  |  \           4          2,3,4
    !           /    |    \
    !         /      |      \
    !       1*-------|-------*3
    !         \      |      /       Each face is characterized by the sum of it local vertex numbers,
    !           \    |    /         i.e. side1=6, side2=7, side3=8, side4=9
    !             \  |  /
    !               \|/
    !                *
    !                2

    orientation(:,:)    = 0                ! orientation vector of adjacent sides
    orientation(1,1:4)  = (/1,3,2,0/)      ! Rows identify the local neighbor's face number
    orientation(2,1:4)  = (/1,2,0,3/)      ! Columns identify the local vertex number of the neighbor
    orientation(3,1:4)  = (/1,0,3,2/)      ! Then the entry in the matrix "orientation" gives the
    orientation(4,1:4)  = (/0,1,2,3/)      ! local vertex number of the neighbor's vertex "Column"
    !                                      ! inside the neighbor's face "Row"
    ! Do the node transformation
    DO iNode = 1, MESH%nNode
      MESH%VRTX%xyNode(:,iNode) = MESH%VRTX%xyNode(:,iNode) + MESH%Displacement(:)
      Temp(:) = MATMUL( MESH%ScalingMatrix(:,:),MESH%VRTX%xyNode(:,iNode) )
      MESH%VRTX%xyNode(:,iNode) = Temp(:)
    ENDDO

    ! =======================================================
    ! Computing neigbhorhood information
    ! =======================================================

    logInfo(*) 'Building neighborhood information ... '
    logInfo(*) 'nummmber of elements is ', MESH%nElem
    !
    MESH%VRTX%NrOfElementsConnected(:) = 0
    !
    test_array_size = 1
    ALLOCATE(test_array1(MESH%nNode,test_array_size))
    test_array1(:,:) = 0
    DO iElem = 1, MESH%nElem
        DO iSide = 1, MESH%GlobalElemType
          !
          iVrtx = MESH%ELEM%Vertex(iSide,iElem)
          MESH%VRTX%NrOfElementsConnected(iVrtx) = MESH%VRTX%NrOfElementsConnected(iVrtx) + 1
          !
          ! Enlargen matrix, if necessary
          IF(MESH%VRTX%NrOfElementsConnected(iVrtx).GT.test_array_size) THEN
             test_array_size = MESH%VRTX%NrOfElementsConnected(iVrtx)
             ALLOCATE(test_array2(MESH%nNode,test_array_size))
             ! Do not use array syntax to copy.
             ! For some strange reasons, array syntax does not
             ! run on most supercomputers
             DO i = 1, MESH%nNode
               DO j = 1, test_array_size-1
                 test_array2(i,j) = test_array1(i,j)
               ENDDO
             test_array2(i,test_array_size) = 0
             ENDDO
             DEALLOCATE(test_array1)
             ALLOCATE(test_array1(MESH%nNode,test_array_size))
             DO i = 1, MESH%nNode
                DO j = 1,test_array_size
                   test_array1(i,j) = test_array2(i,j)
                ENDDO
             ENDDO
             DEALLOCATE(test_array2)
          ENDIF
          !
          test_array1(iVrtx,MESH%VRTX%NrOfElementsConnected(iVrtx)) = iElem
          !
        ENDDO
        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem, ' elements done...'
            END IF
        ENDIF
    ENDDO
    !
    MESH%MaxElementsOnVertex = test_array_size
    !
    DEALLOCATE( MESH%VRTX%Element )
    ALLOCATE(MESH%VRTX%Element(MESH%nNode,MESH%MaxElementsOnVertex) )
    !
    DO iNode = 1, MESH%nNode
      DO iVrtx = 1,MESH%MaxElementsOnVertex
        MESH%VRTX%Element(iNode,iVrtx) = test_array1(iNode,iVrtx)
      ENDDO
    ENDDO
    !
    DEALLOCATE(test_array1)
    !
    logInfo(*) 'Maximum number of tetrahedrons sharing a node: ', MESH%MaxElementsOnVertex
    !
    ! =======================================================
    ! START NEIGHBOUR SEARCH
    ! =======================================================

    ALLOCATE(test_array1(2,MESH%MaxElementsOnVertex))
    ALLOCATE(test_array2(2,MESH%MaxElementsOnVertex))

    DO iElem = 1, MESH%nElem
        DO iSide = 1,MESH%GlobalElemType

            SELECT CASE (iSide)

               CASE(1)

                   ! SEARCH NEIGHBOURS FOR SIDE 1 WITH VERTICES 1-3-2
                   v_tmp(1) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(3,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(2,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(2)

                   ! SEARCH NEIGHBOURS FOR SIDE 2 WITH VERTICES 1-2-4
                   v_tmp(1) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(2,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(4,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(3)

                   ! SEARCH NEIGHBOURS FOR SIDE 3 WITH VERTICES 1-4-3
                   v_tmp(1) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(4,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(3,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(4)
                   ! SEARCH NEIGHBOURS FOR SIDE 4 WITH VERTICES 2-3-4
                   v_tmp(1) = MESH%ELEM%Vertex(2,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(3,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(4,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

            END SELECT
            !
            iColumns  = MINVAL(n_tmp)
            iPosition = MINLOC(n_tmp)
            !
            SELECT CASE (iPosition(1))
               CASE(1)
                   iRow              = v_tmp(1)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(2),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
               CASE(2)
                   iRow              = v_tmp(2)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
               CASE(3)
                   iRow              = v_tmp(3)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(2),:)
            END SELECT
            !
            DO i = 1,iColumns(1)

                test_array2(:,:) = 0
                WHERE(test_array1 .EQ. MESH%VRTX%Element(iRow,i) )
                  test_array2 = 1
                ENDWHERE

                testI = SUM(test_array2)

                IF (testI.EQ.2 .AND. MESH%VRTX%Element(iRow,i).NE.iElem) THEN

                    iNeighbor                           = MESH%VRTX%Element(iRow,i)

                    MESH%ELEM%SideNeighbor(iSide,iElem) = iNeighbor

                    counter = 0

                    DO k = 1,MESH%GlobalElemType - 1
                        DO j = 1,MESH%GlobalElemType
                            IF (MESH%ELEM%Vertex(j,iNeighbor).EQ.v_tmp(k)) THEN
                                counter        = counter + 1
                                n_tmp(counter) = j
                                EXIT
                            END IF
                        ENDDO
                    ENDDO

                    iVrtx = SUM(n_tmp) - 5
                    MESH%ELEM%LocalNeighborSide(iSide,iElem) = iVrtx
                    MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = orientation(iVrtx,n_tmp(1))
                    EXIT

                ENDIF

            ENDDO

        ENDDO

        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem,' elements done...'
            END IF
        ENDIF

    ENDDO
    !
    ! Check the consistency of the mesh
    logInfo(*) 'Checking mesh consistency ... '
    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
          counter = counter + 1
        ENDIF
      ENDDO
    ENDDO
    logInfo(*) 'Total element faces without neighbor: ', counter

    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
          IF(MESH%ELEM%Reference(iSide,iElem).EQ.0) THEN
             MESH%ELEM%Reference(iSide,iElem) = 99999
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! Compute the elements' volumes and barycenters
    logInfo(*) 'Computing element volumes and barycenters. '
    !
    DO iElem = 1, MESH%nELEM
      !
      x1 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1,iElem))
      y1 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1,iElem))
      z1 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1,iElem))
      !
      x2 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(2,iElem))
      y2 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(2,iElem))
      z2 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(2,iElem))
      !
      x3 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(3,iElem))
      y3 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(3,iElem))
      z3 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(3,iElem))
      !
      x4 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(4,iElem))
      y4 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(4,iElem))
      z4 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(4,iElem))
      !
      ! Compute with SARRUS's law from the Jacobi determinant of the transformation
      ! from xyz coordinates to the reference tetrahedron
      MESH%ELEM%Volume(iElem) =   abs(  (x2-x1)*(y3-y1)*(z4-z1) +             &
                                        (x3-x1)*(y4-y1)*(z2-z1) +             &
                                        (x4-x1)*(y2-y1)*(z3-z1) -             &
                                        (x4-x1)*(y3-y1)*(z2-z1) -             &
                                        (x2-x1)*(y4-y1)*(z3-z1) -             &
                                        (x3-x1)*(y2-y1)*(z4-z1)               )/6.0
      ! Barycenters of the tetrahedrons
      MESH%ELEM%xyBary(:,iElem) = 0.25*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem))   )
    ENDDO

    ! Build stencil if necessary
    CALL BND_DG_Periodic3D_Tetra_us(MESH,BND,IO)

    ! Temporarily needed until everything's unified
    ALLOCATE(MESH%ELEM%ncIndex(MESH%nElem,MESH%nSideMax))
    MESH%ELEM%ncIndex    = 0

    DEALLOCATE(test_array1)
    DEALLOCATE(test_array2)

  END SUBROUTINE compute_mesh_Gambit3D_Tetra

  SUBROUTINE compute_mesh_ICEMCFD3D_Tetra(Mesh, IO, EQN, BND,BndCond,TmpBndCond, nBndElemTotal,BndFaceNodes, BndFaceNodesCond)
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH                                          !< Mesh data
    TYPE (tInputOutput)       :: IO                                            !< IO data
    TYPE (tEquations)         :: EQN                                           !< Discretization data
    TYPE (tBoundary)		  :: BND                                           !< Boundary data
	!
    INTEGER, POINTER		  :: BndFaceNodes(:,:)                             !< face and node of boundary information
    INTEGER, POINTER		  :: BndFaceNodesCond(:)                           !< face and node of boundary information
    INTEGER, POINTER		  :: BndCond(:)                                    !< face and node of boundary information
    INTEGER, POINTER		  :: TmpBndCond(:)                                 !< face and node of boundary information
    INTEGER					  :: nBndElemTotal                                 !< total number of boundary elements
    !--------------------------------------------------------------------------
    !local variables
    !--------------------------------------------------------------------------
    INTEGER 				  :: iNode, iVrtx, iElem, iSide, iNeighbor,iBndElem
    INTEGER, POINTER          :: test_array1(:,:)
    INTEGER, POINTER          :: test_array2(:,:)
    REAL                      :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
    REAL                      :: temp(EQN%Dimension)
    INTEGER                   :: test_array_size
    INTEGER,DIMENSION(3)      :: v_tmp,n_tmp
    INTEGER                   :: iRow,iColumns(1),iPosition(1)
    INTEGER					  :: i,j,k,testI, counter
    INTEGER                   :: orientation(4,4)
    INTEGER                   :: ICEMToSeisSol(4)
    INTEGER                   :: nBndObj(6)
    !--------------------------------------------------------------------------
    INTENT(INOUT)			  :: Mesh
    INTENT(IN) 				  :: EQN

    !
    MESH%VRTX%NrOfElementsConnected(:) = 0
    !
    test_array_size = 1
    ALLOCATE(test_array1(MESH%nNode,test_array_size))
    test_array1(:,:) = 0
    DO iElem = 1, MESH%nElem
        DO iSide = 1, MESH%GlobalElemType
          !
          iVrtx = MESH%ELEM%Vertex(iSide,iElem)
          MESH%VRTX%NrOfElementsConnected(iVrtx) = MESH%VRTX%NrOfElementsConnected(iVrtx) + 1
          !
          ! Enlargen matrix, if necessary
          !
          IF(MESH%VRTX%NrOfElementsConnected(iVrtx).GT.test_array_size) THEN
             test_array_size = MESH%VRTX%NrOfElementsConnected(iVrtx)
             ALLOCATE(test_array2(MESH%nNode,test_array_size))
             ! Do not use array syntax to copy.
             ! For some strange reasons, array syntax does not
             ! run on most supercomputers
             DO i = 1, MESH%nNode
               DO j = 1, test_array_size-1
                 test_array2(i,j) = test_array1(i,j)
               ENDDO
             test_array2(i,test_array_size) = 0
             ENDDO
             DEALLOCATE(test_array1)
             ALLOCATE(test_array1(MESH%nNode,test_array_size))
             DO i = 1, MESH%nNode
                DO j = 1,test_array_size
                   test_array1(i,j) = test_array2(i,j)
                ENDDO
             ENDDO
             DEALLOCATE(test_array2)
          ENDIF
          !
          test_array1(iVrtx,MESH%VRTX%NrOfElementsConnected(iVrtx)) = iElem
          !
        ENDDO
        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem,' elements done...'
            END IF
        ENDIF
    ENDDO
    !
    MESH%MaxElementsOnVertex = test_array_size
    !
    ALLOCATE(MESH%VRTX%Element(MESH%nNode,MESH%MaxElementsOnVertex) )
    !
    DO iNode = 1, MESH%nNode
      DO iVrtx = 1,MESH%MaxElementsOnVertex
        MESH%VRTX%Element(iNode,iVrtx) = test_array1(iNode,iVrtx)
      ENDDO
    ENDDO
    !
    DEALLOCATE(test_array1)
    !
    logInfo(*) 'Maximum number of tetrahedrons sharing a node: ', MESH%MaxElementsOnVertex
    !
    !
    ICEMToSeisSol(:) = (/1, 2, 4, 3/)   ! Map ICEM ordering convention for boundary conditions to SeisSol convention
    !
    ! Assuming, that faces consist of the three vertices as follow
    !
    !                              Face       Vertices
    !
    !                4              1          1,3,2
    !                *              2          1,2,4
    !               /|\             3          1,4,3
    !             /  |  \           4          2,3,4
    !           /    |    \
    !         /      |      \
    !       1*-------|-------*3
    !         \      |      /       Each face is characterized by the sum of it local vertex numbers,
    !           \    |    /         i.e. side1=6, side2=7, side3=8, side4=9
    !             \  |  /
    !               \|/
    !                *
    !                2
    !
    orientation(:,:)    = 0                ! orientation vector of adjacent sides
    orientation(1,1:4)  = (/1,3,2,0/)      ! Rows identify the local neighbor's face number
    orientation(2,1:4)  = (/1,2,0,3/)      ! Columns identify the local vertex number of the neighbor
    orientation(3,1:4)  = (/1,0,3,2/)      ! Then the entry in the matrix "orientation" gives the
    orientation(4,1:4)  = (/0,1,2,3/)      ! local vertex number of the neighbor's vertex "Column"
    !                                      ! inside the neighbor's face "Row"
    !
    !
    !START NEIGHBOUR SEARCH
    !
    !
    !
    ALLOCATE(test_array1(2,MESH%MaxElementsOnVertex))
    ALLOCATE(test_array2(2,MESH%MaxElementsOnVertex))
    !
    !
    DO iElem = 1, MESH%nElem
        DO iSide = 1,MESH%GlobalElemType

            SELECT CASE (iSide)

               CASE(1)

                   ! SEARCH NEIGHBOURS FOR SIDE 1 WITH VERTICES 1-3-2
                   v_tmp(1) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(3,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(2,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(2)

                   ! SEARCH NEIGHBOURS FOR SIDE 2 WITH VERTICES 1-2-4
                   v_tmp(1) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(2,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(4,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(3)

                   ! SEARCH NEIGHBOURS FOR SIDE 3 WITH VERTICES 1-4-3
                   v_tmp(1) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(4,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(3,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(4)
                   ! SEARCH NEIGHBOURS FOR SIDE 4 WITH VERTICES 2-3-4
                   v_tmp(1) = MESH%ELEM%Vertex(2,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(3,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(4,iElem)
                   DO k = 1,3
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

            END SELECT
            !
            iColumns  = MINVAL(n_tmp)
            iPosition = MINLOC(n_tmp)
            !
            SELECT CASE (iPosition(1))
               CASE(1)
                   iRow              = v_tmp(1)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(2),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
               CASE(2)
                   iRow              = v_tmp(2)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
               CASE(3)
                   iRow              = v_tmp(3)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(2),:)
            END SELECT
            !
            DO i = 1,iColumns(1)

                test_array2(:,:) = 0
                WHERE(test_array1 .EQ. MESH%VRTX%Element(iRow,i) )
                  test_array2 = 1
                ENDWHERE

                testI = SUM(test_array2)

                IF (testI.EQ.2 .AND. MESH%VRTX%Element(iRow,i).NE.iElem) THEN

                    iNeighbor                           = MESH%VRTX%Element(iRow,i)

                    MESH%ELEM%SideNeighbor(iSide,iElem) = iNeighbor

                    counter = 0

                    DO k = 1,MESH%GlobalElemType - 1
                        DO j = 1,MESH%GlobalElemType
                            IF (MESH%ELEM%Vertex(j,iNeighbor).EQ.v_tmp(k)) THEN
                                counter        = counter + 1
                                n_tmp(counter) = j
                                EXIT
                            END IF
                        ENDDO
                    ENDDO

                    iVrtx = SUM(n_tmp) - 5
                    MESH%ELEM%LocalNeighborSide(iSide,iElem) = iVrtx
                    MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = orientation(iVrtx,n_tmp(1))
                    EXIT

                ENDIF

            ENDDO

        ENDDO

        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem,' elements done...'
            END IF
        ENDIF

    ENDDO
    !
    !
    ! COMPUTE BOUNDARY ELEMENTS AND LOCAL SIDES as this is different in ICEM compared to GAMBIT
    !
    DO iElem = 1, nBndElemTotal


           ! SEARCH ELEMENTS WITH SIDE GIVEN BY BOUNDARY TRIANGLE in BndFaceNodes
           v_tmp(1) = BndFaceNodes(1,iElem)
           v_tmp(2) = BndFaceNodes(2,iElem)
           v_tmp(3) = BndFaceNodes(3,iElem)
           DO k = 1,3
                n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
           ENDDO
           !
           iColumns  = MINVAL(n_tmp)
           iPosition = MINLOC(n_tmp)
           !
           SELECT CASE (iPosition(1))
               CASE(1)
                   iRow              = v_tmp(1)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(2),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
               CASE(2)
                   iRow              = v_tmp(2)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
               CASE(3)
                   iRow              = v_tmp(3)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(2),:)
            END SELECT
            !
            DO i = 1,iColumns(1)

                test_array2(:,:) = 0
                WHERE(test_array1 .EQ. MESH%VRTX%Element(iRow,i) )
                  test_array2 = 1
                ENDWHERE

                testI = SUM(test_array2)
                ! find tetrahedral element that contains this boundary triangle as a face
                IF (testI.EQ.2) THEN

                    iBndElem = MESH%VRTX%Element(iRow,i)
                    ! find the local face number of the tetrahedral element that represents this boundary triangle
                    counter = 0
                    DO k = 1,MESH%GlobalElemType - 1
                        DO j = 1,MESH%GlobalElemType
                            IF (MESH%ELEM%Vertex(j,iBndElem).EQ.v_tmp(k)) THEN
                                counter        = counter + 1
                                n_tmp(counter) = j
                                EXIT
                            END IF
                        ENDDO
                    ENDDO

                    iVrtx = SUM(n_tmp) - 5
                    ! assign the Bondary Condition of this Set to the face of the tetrahedral element
                    MESH%ELEM%Reference(iVrtx,iBndElem) = BndFaceNodesCond(iElem)

                    EXIT

                ENDIF

            ENDDO

            IF(MESH%nElem.GT.20) THEN
              IF (MOD(iElem,FLOOR(nBndElemTotal/20.)).EQ.0) THEN
                logInfo(*) iElem,' boundary elements done...'
              END IF
            ENDIF

    ENDDO
    !
    ! Check boundary conditions
    !
    logInfo(*) 'Checking boundary conditions ... '
    nBndObj(:) = 0.
    DO i = 1, 6
      WHERE(MOD(BndCond(:),100).EQ.i)
        TmpBndCond(:) = BndCond(:)/100
      ELSEWHERE
        TmpBndCond(:) = 0
      ENDWHERE
      nBndObj(i) = MAXVAL(TmpBndCond(:))
      IF(nBndObj(i).NE.BND%NoBndObjects(i)) THEN
         logError(*) 'Error in Gambit3D mesh. '
         logError(*) 'You specified ', nBndObj(i), ' boundaries of type ', i, ' in meshfile and '
         logError(*) 'you specified ', BND%NoBndObjects(i), ' boundaries of type ', i, ' in parameterfile. '
         STOP
      ENDIF
    ENDDO

    logInfo(*) 'Boundary conditions checked. '
    !
    ! Deallocate temp. variables
    !
    DEALLOCATE(BndCond, TmpBndCond)
    !
    !
    ! Check the consistency of the mesh
    !
    logInfo(*) 'Checking mesh consistency ... '
    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
          counter = counter + 1
        ENDIF
      ENDDO
    ENDDO
    logInfo(*) 'Total element faces without neighbor: ', counter

    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
          IF(MESH%ELEM%Reference(iSide,iElem).EQ.0) THEN
             logError(*) 'Error on tetrahedron / side ', iElem, iSide
             logError(*) 'Sideneighbor not found and no boundary condition set!'
             logError(*) 'Perhaps your grid is non-conforming or boundary      '
             logError(*) 'conditions were not set on all the boundaries.       '
             logError(*) 'Please check the BC and the conformity of your mesh. '
             STOP
          ENDIF
          IF(MESH%ELEM%LocalNeighborSide(iSide,iElem).NE.0) THEN
             logError(*) 'Error on tetrahedron / side ', iElem, iSide
             logError(*) 'Sideneighbor not found but localneighborside non zero!'
             STOP
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    logInfo(*) 'Mesh seems to be consistent.'
    !
    ! Compute the elements' volumes and barycenters
    !
    logInfo(*) 'Computing element volumes and barycenters.'
    !
    DO iElem = 1, MESH%nELEM
      !
      x1 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1,iElem))
      y1 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1,iElem))
      z1 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1,iElem))
      !
      x2 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(2,iElem))
      y2 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(2,iElem))
      z2 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(2,iElem))
      !
      x3 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(3,iElem))
      y3 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(3,iElem))
      z3 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(3,iElem))
      !
      x4 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(4,iElem))
      y4 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(4,iElem))
      z4 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(4,iElem))
      !
      ! Compute with SARRUS's law from the Jacobi determinant of the transformation
      ! from xyz coordinates to the reference tetrahedron
      !
      MESH%ELEM%Volume(iElem) =   abs(  (x2-x1)*(y3-y1)*(z4-z1) +             &
                                        (x3-x1)*(y4-y1)*(z2-z1) +             &
                                        (x4-x1)*(y2-y1)*(z3-z1) -             &
                                        (x4-x1)*(y3-y1)*(z2-z1) -             &
                                        (x2-x1)*(y4-y1)*(z3-z1) -             &
                                        (x3-x1)*(y2-y1)*(z4-z1)               )/6.0D0
      ! Barycenters of the tetrahedrons
      MESH%ELEM%xyBary(:,iElem) = 0.25D0*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem))   )
      !

    ENDDO
    !
    ! Build stencil if necessary
    CALL BND_DG_Periodic3D_Tetra_us(MESH,BND,IO)

    ! Temporarily needed until everything's unified
    ALLOCATE(MESH%ELEM%ncIndex(MESH%nElem,MESH%nSideMax))
    MESH%ELEM%ncIndex    = 0

  END SUBROUTINE compute_mesh_ICEMCFD3D_Tetra

  SUBROUTINE compute_mesh_Gambit3D_Hexa(Mesh, IO, EQN, BND)
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH                                          !< Mesh data
    TYPE (tInputOutput)       :: IO                                            !< IO data
    TYPE (tEquations)         :: EQN                                           !< Discretization data
    TYPE (tBoundary)		  :: BND                                           !< Boundary data
    !--------------------------------------------------------------------------
    !local variables
    !--------------------------------------------------------------------------
    INTEGER 				  :: iNode, iVrtx, iElem, iSide, iLocalFace, iNeighbor
    INTEGER, POINTER          :: test_array1(:,:)
    INTEGER, POINTER          :: test_array2(:,:)
    REAL                      :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
    REAL                      :: x5,y5,z5, x6,y6,z6, x7,y7,z7, x8,y8,z8
    REAL                      :: xt1,yt1,zt1, xt2,yt2,zt2, xt3,yt3,zt3, xt4,yt4,zt4
    REAL                      :: xt5,yt5,zt5, xt6,yt6,zt6, xt7,yt7,zt7, xt8,yt8,zt8
    REAL                      :: temp(EQN%Dimension)
    INTEGER                   :: test_array_size
    INTEGER,DIMENSION(4)      :: v_tmp,n_tmp
    INTEGER                   :: iRow,iColumns(1),iPosition(1)
    INTEGER					  :: i,j,k,testI, counter
    INTEGER                   :: orientation(6,8)
    INTEGER					  :: m
    INTEGER,DIMENSION(6)      :: side_sum
    !--------------------------------------------------------------------------
    INTENT(INOUT)			  :: Mesh
    INTENT(IN) 				  :: EQN

    !
    ! Assuming, that faces consist of the four vertices as defined by Gambit
    !
    !         7           8         1          1,2,6,5
    !         *-----------*         2          2,4,8,6
    !        /|          /|         3          4,3,7,8
    !       / |         / |         4          3,1,5,7
    !     5*-----------*6 |         5          2,1,3,4
    !      |  |        |  |         6          5,6,8,7
    !      |  |        |  |
    !      | 3*--------|--*4        Each face is characterized by the sum of it local vertex numbers,
    !      | /         | /          i.e. side1=14, side2=20, side3=22, side4=16, side5=10, side6=26
    !      |/          |/
    !     1*-----------*2
    !
    !
    orientation(:,:)    = 0                        ! orientation vector of adjacent sides
    orientation(1,1:8)  = (/1,2,0,0,4,3,0,0/)      ! Rows identify the local neighbor's face number
    orientation(2,1:8)  = (/0,1,0,2,0,4,0,3/)      ! Columns identify the local vertex number of the neighbor
    orientation(3,1:8)  = (/0,0,2,1,0,0,3,4/)      ! Then the entry in the matrix "orientation" gives the
    orientation(4,1:8)  = (/2,0,1,0,3,0,4,0/)      ! local vertex number of the neighbor's vertex "Column"
    orientation(5,1:8)  = (/2,1,3,4,0,0,0,0/)      ! inside the neighbor's face "Row"
    orientation(6,1:8)  = (/0,0,0,0,1,2,4,3/)

    side_sum = (/14,20,22,16,10,26/)

    MESH%VRTX%NrOfElementsConnected(:) = 0
    !
    test_array_size = 1
    ALLOCATE(test_array1(MESH%nNode,test_array_size))
    test_array1(:,:) = 0
    DO iElem = 1, MESH%nElem
        DO iSide = 1, MESH%GlobalElemType + 2  !loop goes over the 8 vertices of a hexagonal element
          !
          iVrtx = MESH%ELEM%Vertex(iSide,iElem)
          MESH%VRTX%NrOfElementsConnected(iVrtx) = MESH%VRTX%NrOfElementsConnected(iVrtx) + 1
          !
          ! Enlargen matrix, if necessary
          !
          IF(MESH%VRTX%NrOfElementsConnected(iVrtx).GT.test_array_size) THEN
             test_array_size = MESH%VRTX%NrOfElementsConnected(iVrtx)
             ALLOCATE(test_array2(MESH%nNode,test_array_size))
             ! Do not use array syntax to copy.
             ! For some strange reasons, array syntax does not
             ! run on most supercomputers
             DO i = 1, MESH%nNode
               DO j = 1, test_array_size-1
                 test_array2(i,j) = test_array1(i,j)
               ENDDO
             test_array2(i,test_array_size) = 0
             ENDDO
             DEALLOCATE(test_array1)
             ALLOCATE(test_array1(MESH%nNode,test_array_size))
             DO i = 1, MESH%nNode
                DO j = 1,test_array_size
                   test_array1(i,j) = test_array2(i,j)
                ENDDO
             ENDDO
             DEALLOCATE(test_array2)
          ENDIF
          !
          test_array1(iVrtx,MESH%VRTX%NrOfElementsConnected(iVrtx)) = iElem
          !
        ENDDO
        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem,' elements done...'
            END IF
        ENDIF
    ENDDO
    !
    MESH%MaxElementsOnVertex = test_array_size
    !
    ALLOCATE(MESH%VRTX%Element(MESH%nNode,MESH%MaxElementsOnVertex) )
    !
    DO iNode = 1, MESH%nNode
      DO iVrtx = 1,MESH%MaxElementsOnVertex
        MESH%VRTX%Element(iNode,iVrtx) = test_array1(iNode,iVrtx)
      ENDDO
    ENDDO

    !
    DEALLOCATE(test_array1)
    !
    logInfo(*) 'Maximum number of hexahedrons sharing a node: ', MESH%MaxElementsOnVertex
    !
    !
    !
    !START NEIGHBOUR SEARCH
    !
    !
    !
    ALLOCATE(test_array1(3,MESH%MaxElementsOnVertex))
    ALLOCATE(test_array2(3,MESH%MaxElementsOnVertex))
    !
    !
    !
    DO iElem = 1, MESH%nElem
        DO iSide = 1,MESH%GlobalElemType

            SELECT CASE (iSide)

               CASE(1)

                   ! SEARCH NEIGHBOURS FOR SIDE 1 WITH VERTICES 1-2-6-5
                   v_tmp(1) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(2,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(6,iElem)
                   v_tmp(4) = MESH%ELEM%Vertex(5,iElem)
                   DO k = 1,4
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(2)

                   ! SEARCH NEIGHBOURS FOR SIDE 2 WITH VERTICES 2-4-8-6
                   v_tmp(1) = MESH%ELEM%Vertex(2,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(4,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(8,iElem)
                   v_tmp(4) = MESH%ELEM%Vertex(6,iElem)
                   DO k = 1,4
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(3)

                   ! SEARCH NEIGHBOURS FOR SIDE 3 WITH VERTICES 4-3-7-8
                   v_tmp(1) = MESH%ELEM%Vertex(4,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(3,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(7,iElem)
                   v_tmp(4) = MESH%ELEM%Vertex(8,iElem)
                   DO k = 1,4
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(4)

                   ! SEARCH NEIGHBOURS FOR SIDE 4 WITH VERTICES 3-1-5-7
                   v_tmp(1) = MESH%ELEM%Vertex(3,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(5,iElem)
                   v_tmp(4) = MESH%ELEM%Vertex(7,iElem)
                   DO k = 1,4
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(5)

                   ! SEARCH NEIGHBOURS FOR SIDE 3 WITH VERTICES 2-1-3-4
                   v_tmp(1) = MESH%ELEM%Vertex(2,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(1,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(3,iElem)
                   v_tmp(4) = MESH%ELEM%Vertex(4,iElem)
                   DO k = 1,4
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

               CASE(6)

                   ! SEARCH NEIGHBOURS FOR SIDE 4 WITH VERTICES 5-6-8-7
                   v_tmp(1) = MESH%ELEM%Vertex(5,iElem)
                   v_tmp(2) = MESH%ELEM%Vertex(6,iElem)
                   v_tmp(3) = MESH%ELEM%Vertex(8,iElem)
                   v_tmp(4) = MESH%ELEM%Vertex(7,iElem)
                   DO k = 1,4
                        n_tmp(k) = MESH%VRTX%NrOfElementsConnected(v_tmp(k))
                   ENDDO

            END SELECT
            !
            iColumns  = MINVAL(n_tmp)
            iPosition = MINLOC(n_tmp)
            !
            SELECT CASE (iPosition(1))
               CASE(1)
                   iRow              = v_tmp(1)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(2),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp(4),:)
               CASE(2)
                   iRow              = v_tmp(2)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(3),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp(4),:)
               CASE(3)
                   iRow              = v_tmp(3)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(2),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp(4),:)
               CASE(4)
                   iRow              = v_tmp(4)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp(2),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp(3),:)
            END SELECT
            !
            DO i = 1,iColumns(1)

                test_array2(:,:) = 0
                WHERE(test_array1 .EQ. MESH%VRTX%Element(iRow,i) )
                  test_array2 = 1
                ENDWHERE

                testI = SUM(test_array2)

                IF (testI.EQ.3 .AND. MESH%VRTX%Element(iRow,i).NE.iElem) THEN

                    iNeighbor                           = MESH%VRTX%Element(iRow,i)

                    MESH%ELEM%SideNeighbor(iSide,iElem) = iNeighbor

                    counter = 0
                    DO k = 1,MESH%GlobalElemType - 2         ! loop over all four  face vertices
                        DO j = 1,MESH%GlobalElemType + 2     ! loop over all eight neighbor vertices
                            IF (MESH%ELEM%Vertex(j,iNeighbor).EQ.v_tmp(k)) THEN
                                counter        = counter + 1
                                n_tmp(counter) = j
                                EXIT
                            END IF
                        ENDDO
                    ENDDO
                    iLocalFace = SUM(n_tmp)
                    DO m = 1,6
                       IF (side_sum(m).EQ.iLocalFace) THEN
                           iLocalFace = m
                           EXIT
                       ENDIF
                    ENDDO
                    MESH%ELEM%LocalNeighborSide(iSide,iElem) = iLocalFace
                    MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = orientation(iLocalFace,n_tmp(1))
                    EXIT

                ENDIF

            ENDDO

        ENDDO

        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem,' elements done...'
            END IF
        ENDIF

    ENDDO
    !
    ! Check the consistency of the mesh
    !
    logInfo(*) 'Checking mesh consistency ... '
    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
          counter = counter + 1
        ENDIF
      ENDDO
    ENDDO
    logInfo(*) 'Total element faces without neighbor: ', counter

    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
          IF(MESH%ELEM%Reference(iSide,iElem).EQ.0) THEN
             logError(*) 'Error on hexahedron / side ', iElem, iSide
             logError(*) 'Sideneighbor not found and no boundary condition set!'
             logError(*) 'Perhaps your mesh is non-conforming or boundary      '
             logError(*) 'conditions were not set on all the boundaries.       '
             logError(*) 'Please check the BC and the conformity of your mesh. '
             STOP
          ENDIF
          IF(MESH%ELEM%LocalNeighborSide(iSide,iElem).NE.0) THEN
             logError(*) 'Error on hexahedron / side ', iElem, iSide
             logError(*) 'Sideneighbor not found but localneighborside non zero!'
             STOP
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    logInfo(*) 'Mesh seems to be consistent. '
    !
    ! Compute the elements' volumes and barycenters
    !
    logInfo(*) 'Computing element volumes and barycenters. '
    !
    DO iElem = 1, MESH%nELEM
      !
      x1 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1,iElem))
      y1 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1,iElem))
      z1 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1,iElem))
      !
      x2 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(2,iElem))
      y2 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(2,iElem))
      z2 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(2,iElem))
      !
      x3 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(3,iElem))
      y3 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(3,iElem))
      z3 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(3,iElem))
      !
      x4 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(4,iElem))
      y4 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(4,iElem))
      z4 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(4,iElem))
      !
      x5 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(5,iElem))
      y5 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(5,iElem))
      z5 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(5,iElem))
      !
      x6 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(6,iElem))
      y6 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(6,iElem))
      z6 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(6,iElem))
      !
      x7 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(7,iElem))
      y7 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(7,iElem))
      z7 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(7,iElem))
      !
      x8 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(8,iElem))
      y8 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(8,iElem))
      z8 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(8,iElem))
      !
      ! Compute the hexahedron's volume by the sum of the 5 tetrahedrons,
      ! that subdivide it. The tetrahedron's volume is computed with SARRUS's law
      ! as before
      !
      MESH%ELEM%Volume(iElem) = 0.
      DO m = 1,5

          SELECT CASE (m)
          CASE(1)
             xt1 = x5; xt2 = x7; xt3 = x6; xt4 = x1;
             yt1 = y5; yt2 = y7; yt3 = y6; yt4 = y1;
             zt1 = z5; zt2 = z7; zt3 = z6; zt4 = z1;
          CASE(2)
             xt1 = x7; xt2 = x3; xt3 = x4; xt4 = x1;
             yt1 = y7; yt2 = y3; yt3 = y4; yt4 = y1;
             zt1 = z7; zt2 = z3; zt3 = z4; zt4 = z1;
          CASE(3)
             xt1 = x6; xt2 = x8; xt3 = x4; xt4 = x7;
             yt1 = y6; yt2 = y8; yt3 = y4; yt4 = y7;
             zt1 = z6; zt2 = z8; zt3 = z4; zt4 = z7;
          CASE(4)
             xt1 = x6; xt2 = x4; xt3 = x2; xt4 = x1;
             yt1 = y6; yt2 = y4; yt3 = y2; yt4 = y1;
             zt1 = z6; zt2 = z4; zt3 = z2; zt4 = z1;
          CASE(5)
             xt1 = x6; xt2 = x1; xt3 = x7; xt4 = x4;
             yt1 = y6; yt2 = y1; yt3 = y7; yt4 = y4;
             zt1 = z6; zt2 = z1; zt3 = z7; zt4 = z4;
          END SELECT

          MESH%ELEM%Volume(iElem) = MESH%ELEM%Volume(iElem)                 +    &
                                        abs ( (xt2-xt1)*(yt3-yt1)*(zt4-zt1) +    &
                                              (xt3-xt1)*(yt4-yt1)*(zt2-zt1) +    &
                                              (xt4-xt1)*(yt2-yt1)*(zt3-zt1) -    &
                                              (xt4-xt1)*(yt3-yt1)*(zt2-zt1) -    &
                                              (xt2-xt1)*(yt4-yt1)*(zt3-zt1) -    &
                                              (xt3-xt1)*(yt2-yt1)*(zt4-zt1)      )/6.0D0

      ENDDO
      ! Barycenters of the hexahedrons
      MESH%ELEM%xyBary(:,iElem) = 1./8.*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(5,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(6,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(7,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(8,iElem))   )

    ENDDO
    !
    ! Build stencil if necessary
    CALL BND_DG_Periodic3D_Hexa_us(MESH,BND,IO)

    ! Temporarily needed until everything's unified
    ALLOCATE(MESH%ELEM%ncIndex(MESH%nElem,MESH%nSideMax))
    MESH%ELEM%ncIndex    = 0

  END SUBROUTINE compute_mesh_Gambit3D_Hexa

  !> compute mesh information and search neighbour for Gambit 3D Mixed mesh
  !> dynamic rupture is considered
  !<
  SUBROUTINE compute_mesh_Gambit3D_Mixed(Mesh, IO, EQN, BND, DISC, TmpFaultElements)
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH                                          !< Mesh data
    TYPE (tInputOutput)       :: IO                                            !< IO data
    TYPE (tEquations)         :: EQN                                           !< Discretization data
    TYPE (tBoundary)		  :: BND                                           !< Boundary data
    TYPE (tDiscretization)    :: DISC                                          !< Discretization data

    !< Dynamic Rupture variables
    INTEGER					  :: TmpFaultElements                              !< temporarily number of fault elements
    !--------------------------------------------------------------------------
    !local variables
    !--------------------------------------------------------------------------
    INTEGER 				  :: iNode, iVrtx, iElem, iSide, iNeighbor
    INTEGER, POINTER          :: test_array1(:,:)
    INTEGER, POINTER          :: test_array2(:,:)
    REAL                      :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
    REAL                      :: x5,y5,z5, x6,y6,z6, x7,y7,z7, x8,y8,z8
    REAL                      :: xt1, xt2, xt3, xt4, yt1, yt2, yt3, yt4
    REAL                      :: zt1, zt2, zt3, zt4
    REAL                      :: temp(EQN%Dimension)
    INTEGER                   :: test_array_size
    INTEGER,DIMENSION(3)      :: v_tmp,n_tmp
    INTEGER                   :: iRow,iColumns(1),iPosition(1)
    INTEGER					  :: i,j,k,testI, counter
    INTEGER                   :: LocElemType
    INTEGER,DIMENSION(3)      :: v_tmp_tet, n_tmp_tet
    INTEGER,DIMENSION(4)      :: v_tmp_hex, n_tmp_hex
    INTEGER					  :: iNCboundary, m
    REAL,DIMENSION(8)         :: x,y,z
    INTEGER                   :: VertexSide_Tet(4,3)
    INTEGER                   :: VertexSide_Hex(6,4)
    INTEGER                   :: orientation_tet(4,4)
    INTEGER                   :: orientation_hex(6,8)
    INTEGER                   :: eType, nVrtx
    LOGICAL                   :: exit_loop=.false.
    INTEGER					  :: nIntGP
    INTEGER                   :: nBNDGP, iBNDGP
    REAL                      :: xi, eta, zeta
    REAL                      :: chiGP, tauGP
    LOGICAL                   :: mirrored_point_in_gap
    REAL                      :: p_new(EQN%Dimension)
    REAL                      :: sidevec(EQN%Dimension,2)
    INTEGER,DIMENSION(6)      :: side_sum
    ! Dynamic Rupture variables
    INTEGER					  :: count
    INTEGER, ALLOCATABLE      :: tmp_el(:,:)
    INTEGER                   :: iElem2, iSide2
    REAL                      :: Dis1, Dis2, XRef, YRef, ZRef
    LOGICAL                   :: cycle
    REAL                      :: Bary1(3),Bary2(3)
    REAL                      :: Length
#ifdef GENERATEDKERNELS
    !< statistics about mesh orientations
    integer                   :: l_meshStatistics(0:51)
    integer                   :: l_combinationIndex
#endif
    !--------------------------------------------------------------------------
    INTENT(INOUT)			  :: Mesh
    INTENT(IN) 				  :: EQN


	VertexSide_Tet(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
	VertexSide_Tet(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
	VertexSide_Tet(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
	VertexSide_Tet(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

	VertexSide_Hex(1,:) =  (/ 1,2,6,5 /)   ! Local hex. vertices of side I        !
	VertexSide_Hex(2,:) =  (/ 2,4,8,6 /)   ! Local hex. vertices of side II       !
	VertexSide_Hex(3,:) =  (/ 4,3,7,8 /)   ! Local hex. vertices of side III      !
	VertexSide_Hex(4,:) =  (/ 3,1,5,7 /)   ! Local hex. vertices of side IV       !
	VertexSide_Hex(5,:) =  (/ 2,1,3,4 /)   ! Local hex. vertices of side V        !
	VertexSide_Hex(6,:) =  (/ 5,6,8,7 /)   ! Local hex. vertices of side VI       !
    !
    ! Assuming, that faces consist of the three vertices as follow
    !
    !                              Face       Vertices
    !
    !                4              1          1,3,2
    !                *              2          1,2,4
    !               /|\             3          1,4,3
    !             /  |  \           4          2,3,4
    !           /    |    \
    !         /      |      \
    !       1*-------|-------*3
    !         \      |      /       Each face is characterized by the sum of its local vertex numbers,
    !           \    |    /         i.e. side1=6, side2=7, side3=8, side4=9
    !             \  |  /
    !               \|/
    !                *
    !                2
    !
    orientation_tet(:,:)    = 0                ! orientation vector of adjacent sides
    orientation_tet(1,1:4)  = (/1,3,2,0/)      ! Rows identify the local neighbor's face number
    orientation_tet(2,1:4)  = (/1,2,0,3/)      ! Columns identify the local vertex number of the neighbor
    orientation_tet(3,1:4)  = (/1,0,3,2/)      ! Then the entry in the matrix "orientation" gives the
    orientation_tet(4,1:4)  = (/0,1,2,3/)      ! local vertex number of the neighbor's vertex "Column"
    !                                          ! inside the neighbor's face "Row"
    !
    !

    !
    ! Assuming, that faces consist of the four vertices as defined by Gambit
    !
    !         7           8         1          1,2,6,5
    !         *-----------*         2          2,4,8,6
    !        /|          /|         3          4,3,7,8
    !       / |         / |         4          3,1,5,7
    !     5*-----------*6 |         5          2,1,3,4
    !      |  |        |  |         6          5,6,8,7
    !      |  |        |  |
    !      | 3*--------|--*4        Each face is characterized by the sum of it local vertex numbers,
    !      | /         | /          i.e. side1=14, side2=20, side3=22, side4=16, side5=10, side6=26
    !      |/          |/
    !     1*-----------*2
    !
    !
    orientation_hex(:,:)    = 0                        ! orientation vector of adjacent sides
    orientation_hex(1,1:8)  = (/1,2,0,0,4,3,0,0/)      ! Rows identify the local neighbor's face number
    orientation_hex(2,1:8)  = (/0,1,0,2,0,4,0,3/)      ! Columns identify the local vertex number of the neighbor
    orientation_hex(3,1:8)  = (/0,0,2,1,0,0,3,4/)      ! Then the entry in the matrix "orientation" gives the
    orientation_hex(4,1:8)  = (/2,0,1,0,3,0,4,0/)      ! local vertex number of the neighbor's vertex "Column"
    orientation_hex(5,1:8)  = (/2,1,3,4,0,0,0,0/)      ! inside the neighbor's face "Row"
    orientation_hex(6,1:8)  = (/0,0,0,0,1,2,4,3/)

    side_sum = (/14,20,22,16,10,26/)

	MESH%VRTX%NrOfElementsConnected(:) = 0
    !
    test_array_size = 1
    ALLOCATE(test_array1(MESH%nNode,test_array_size))
    test_array1(:,:) = 0
    DO iElem = 1, MESH%nElem
        DO iVrtx = 1, MESH%LocalVrtxType(iElem)
          !
          iNode = MESH%ELEM%Vertex(iVrtx,iElem)
          MESH%VRTX%NrOfElementsConnected(iNode) = MESH%VRTX%NrOfElementsConnected(iNode) + 1
          !
          ! Enlargen matrix, if necessary
          !
          IF(MESH%VRTX%NrOfElementsConnected(iNode).GT.test_array_size) THEN
             test_array_size = MESH%VRTX%NrOfElementsConnected(iNode)
             ALLOCATE(test_array2(MESH%nNode,test_array_size))
             ! Do not use array syntax to copy.
             ! For some strange reasons, array syntax does not
             ! run on most supercomputers
             DO i = 1, MESH%nNode
               DO j = 1, test_array_size-1
                 test_array2(i,j) = test_array1(i,j)
               ENDDO
             test_array2(i,test_array_size) = 0
             ENDDO
             DEALLOCATE(test_array1)
             ALLOCATE(test_array1(MESH%nNode,test_array_size))
             DO i = 1, MESH%nNode
                DO j = 1,test_array_size
                   test_array1(i,j) = test_array2(i,j)
                ENDDO
             ENDDO
             DEALLOCATE(test_array2)
          ENDIF
          !
          test_array1(iNode,MESH%VRTX%NrOfElementsConnected(iNode)) = iElem
          !
        ENDDO
        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem,' elements done...'
            END IF
        ENDIF
    ENDDO
    !
    MESH%MaxElementsOnVertex = test_array_size
    !
    ALLOCATE(MESH%VRTX%Element(MESH%nNode,MESH%MaxElementsOnVertex) )
    !
    DO iNode = 1, MESH%nNode
      DO iVrtx = 1,MESH%MaxElementsOnVertex
        MESH%VRTX%Element(iNode,iVrtx) = test_array1(iNode,iVrtx)
      ENDDO
    ENDDO
    !
    DEALLOCATE(test_array1)
    !
    logInfo(*) 'Maximum number of elements sharing a node: ', MESH%MaxElementsOnVertex
    !
    !
    !
    !START NEIGHBOUR SEARCH FOR CONFORMING ELEMENTS
    !
    !
    !
    !
    !
    ! Take the node with lowest amount of connections for the neighbor search
    DO iElem = 1, MESH%nElem
        LocElemType = MESH%LocalElemType(iElem)
        ALLOCATE(test_array1(LocElemType/2,MESH%MaxElementsOnVertex))
        ALLOCATE(test_array2(LocElemType/2,MESH%MaxElementsOnVertex))
        !
        DO iSide = 1,MESH%LocalElemType(iElem)
            SELECT CASE(MESH%LocalElemType(iElem))
            CASE(4) ! v = number of node in mesh file, n = number of connections
                v_tmp_tet(1) = MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem)
                v_tmp_tet(2) = MESH%ELEM%Vertex(VertexSide_Tet(iSide,2),iElem)
                v_tmp_tet(3) = MESH%ELEM%Vertex(VertexSide_Tet(iSide,3),iElem)
                DO k = 1,3
                    n_tmp_tet(k) = MESH%VRTX%NrOfElementsConnected(v_tmp_tet(k))
                ENDDO
                !
                iColumns  = MINVAL(n_tmp_tet)
                iPosition = MINLOC(n_tmp_tet)
                !
                SELECT CASE (iPosition(1))
                CASE(1)
                   iRow              = v_tmp_tet(1)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp_tet(2),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp_tet(3),:)
                CASE(2)
                   iRow              = v_tmp_tet(2)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp_tet(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp_tet(3),:)
                CASE(3)
                   iRow              = v_tmp_tet(3)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp_tet(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp_tet(2),:)
            END SELECT
            CASE(6)
                v_tmp_hex(1) = MESH%ELEM%Vertex(VertexSide_Hex(iSide,1),iElem)
                v_tmp_hex(2) = MESH%ELEM%Vertex(VertexSide_Hex(iSide,2),iElem)
                v_tmp_hex(3) = MESH%ELEM%Vertex(VertexSide_Hex(iSide,3),iElem)
                v_tmp_hex(4) = MESH%ELEM%Vertex(VertexSide_Hex(iSide,4),iElem)
                DO k = 1,4
                    n_tmp_hex(k) = MESH%VRTX%NrOfElementsConnected(v_tmp_hex(k))
                ENDDO
                !
                iColumns  = MINVAL(n_tmp_hex)
                iPosition = MINLOC(n_tmp_hex)
                !
                SELECT CASE (iPosition(1))
                CASE(1)
                   iRow              = v_tmp_hex(1)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp_hex(2),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp_hex(3),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp_hex(4),:)
                CASE(2)
                   iRow              = v_tmp_hex(2)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp_hex(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp_hex(3),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp_hex(4),:)
                CASE(3)
                   iRow              = v_tmp_hex(3)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp_hex(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp_hex(2),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp_hex(4),:)
                CASE(4)
                   iRow              = v_tmp_hex(4)
                   test_array1(1,:)  = MESH%VRTX%Element(v_tmp_hex(1),:)
                   test_array1(2,:)  = MESH%VRTX%Element(v_tmp_hex(2),:)
                   test_array1(3,:)  = MESH%VRTX%Element(v_tmp_hex(3),:)
                END SELECT
            END SELECT

            !
            DO i = 1,iColumns(1)

                test_array2(:,:) = 0
                WHERE(test_array1 .EQ. MESH%VRTX%Element(iRow,i) )
                  test_array2 = 1
                ENDWHERE

                testI = SUM(test_array2)
                ! for tets we need two additional equal nodes, for hexas three
                ! (same as LocalElemType/2)
                IF (testI.EQ.0.5*MESH%LocalElemType(iElem) .AND. MESH%VRTX%Element(iRow,i).NE.iElem) THEN

                    iNeighbor = MESH%VRTX%Element(iRow,i)

                    MESH%ELEM%SideNeighbor(iSide,iElem) = iNeighbor                    

                    counter = 0

                    SELECT CASE(MESH%LocalElemType(iElem))
                    CASE(4)
                        DO k = 1,MESH%LocalElemType(iElem) - 1
                            DO j = 1,MESH%LocalVrtxType(iElem)
                                IF (MESH%ELEM%Vertex(j,iNeighbor).EQ.v_tmp_tet(k)) THEN
                                    counter = counter + 1
                                    n_tmp_tet(counter) = j
                                    EXIT
                                END IF
                            ENDDO
                        ENDDO

                        iVrtx = SUM(n_tmp_tet) - 5
                        MESH%ELEM%LocalNeighborSide(iSide,iElem) = iVrtx
                        MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = orientation_tet(iVrtx,n_tmp_tet(1))
                        EXIT
                    CASE(6)
                        DO k = 1,MESH%LocalElemType(iElem) - 2
                            DO j = 1,MESH%LocalVrtxType(iElem)
                                IF (MESH%ELEM%Vertex(j,iNeighbor).EQ.v_tmp_hex(k)) THEN
                                    counter = counter + 1
                                    n_tmp_hex(counter) = j
                                    EXIT
                                END IF
                            ENDDO
                        ENDDO

                        iVrtx = SUM(n_tmp_hex)
                        DO m = 1,6
                            IF(side_sum(m).EQ.iVrtx)THEN
                                iVrtx = m
                                EXIT
                            ENDIF
                        ENDDO
                        MESH%ELEM%LocalNeighborSide(iSide,iElem) = iVrtx
                        MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = orientation_hex(iVrtx,n_tmp_hex(1))
                        EXIT
                    END SELECT
                ENDIF

            ENDDO ! number of connections

        ENDDO ! iElem

        IF(MESH%nElem.GT.20) THEN
            IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                logInfo(*) iElem,' elements done...'
            END IF
        ENDIF
        !
        DEALLOCATE(test_array1)
        DEALLOCATE(test_array2)
        !
    ENDDO

! Already done above?!?
    !Compute local neighbor side (index of side from neighbor's point of view)
    DO iElem = 1, MESH%nElem
       !
       DO iSide = 1, MESH%LocalElemType(iElem)
          !
          iNeighbor  = MESH%ELEM%SideNeighbor(iSide,iElem)
          IF(iNeighbor.GT.MESH%nElem) THEN                         ! If side is on the
             CYCLE                                                 ! boundary, CYCLE
          ENDIF
          !
          DO k = 1, MESH%LocalElemType(iNeighbor)
             IF(MESH%ELEM%SideNeighbor(k,iNeighbor).EQ.iElem) THEN
                MESH%ELEM%LocalNeighborSide(iSide,iElem) = k
             ENDIF
          ENDDO
          !
       ENDDO
    ENDDO


    !
    ! Check the consistency of the mesh
    !
    logInfo(*) 'Checking mesh consistency ...'
        
    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%LocalElemType(iElem)
        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
          counter = counter + 1
        ENDIF
      ENDDO
    ENDDO
    ! Here only really non-conforming boundaries are counted
    logInfo(*) 'Total element faces without neighbors: ', counter

#ifdef GENERATEDKERNELS
#ifdef PARALLEL
    logWarning0(*) 'Mesh statistics currently only working for the serial code'
#endif
    ! initialize to zero
    l_meshStatistics = 0;

    ! compute and print mesh statistics
    do iElem = 1, mesh%nElem
      do iSide = 1, 4
        ! element local matrices
        l_meshStatistics(iSide-1) = l_meshStatistics(iSide-1) + 1

        ! neigboring matrices
        ! default case: We have a neighboring element
        if( (mesh%elem%reference(iSide,iElem) .ne. 101) .and. (mesh%elem%reference(iSide,iElem) .ne. 105) ) then
          ! compute index
          l_combinationIndex =  4&                                                    ! jump over local matrices
                              + ( iSide-1 )*12&                                       ! jump over index \f$i\f$
                              + ( mesh%elem%localNeighborSide(iSide,iElem) - 1 ) * 3& ! jump over index \f$j\f$
                              + mesh%elem%localNeighborVrtx(iSide,iElem) - 1          ! jump over index \f$h\f$
#ifndef NDEBUG
#ifndef PARALLEL
          ! assert we have a valid index
          if( (l_combinationIndex .le. 3) .or. (l_combinationIndex .ge. 52) ) then
              logError(*) mesh%elem%reference(iSide,iElem), iSide, mesh%elem%localNeighborSide(iSide,iElem), mesh%elem%localNeighborVrtx(iSide,iElem)
              logError(*) 'Invalid index for neighboring flux matrix.', l_combinationIndex
              stop
          endif
#endif
#endif
          l_meshStatistics(l_combinationIndex) = l_meshStatistics(l_combinationIndex) + 1
        
        ! free surface boundary conditions: the corresponding local matrix is called
        else if( mesh%elem%reference(iSide,iElem) .eq. 101 ) then
          l_meshStatistics(iSide-1) = l_meshStatistics(iSide-1) + 1

        ! absobing boundary condtions: nothing to do here
#ifndef NDEBUG
        else
          if( mesh%elem%reference(iSide,iElem) .ne. 105 ) then
            logError(*) 'Invalid boundary conditions for generated kernels.', mesh%elem%reference(iSide,iElem)
            stop
          endif
#endif
        endif
      enddo
    enddo

    logDebug(*) 'Generated Kernels: Printing mesh statistics for local and neighboring flux matrices with face and vertex combinations - id (111, 112, 113, [...], 443), #faces.'
    do l_combinationIndex = 0,51
      logDebug(*) l_combinationIndex, l_meshStatistics(l_combinationIndex)
    enddo
#endif

! THIS PART NEEDS TO COMMENTED IF MORE THAN 1 CPU IS USED FOR NEW SWEEPING ALGORITHM (LOW MEM CONSUMPTION)
!    DO iElem = 1, MESH%nElem
!      DO iSide = 1, MESH%LocalElemType(iElem)
!        IF(MESH%ELEM%SideNeighbor(iSide,iElem).GT.MESH%nElem) THEN
!          IF(MESH%ELEM%Reference(iSide,iElem).EQ.0) THEN
!             WRITE(IO%UNIT%errOut,*) 'Error on element / side ', iElem, iSide
!             WRITE(IO%UNIT%errOut,*) 'Sideneighbor not found and no boundary condition set!'
!             WRITE(IO%UNIT%errOut,*) 'Perhaps your boundary      '
!             WRITE(IO%UNIT%errOut,*) 'conditions were not set on all the boundaries.       '
!             WRITE(IO%UNIT%errOut,*) 'Please check the BC and the conformity of your mesh. '
!             STOP
!          ENDIF
!          IF(MESH%ELEM%LocalNeighborSide(iSide,iElem).NE.0) THEN
!             WRITE(IO%UNIT%errOut,*) 'Error on element / side ', iElem, iSide
!             WRITE(IO%UNIT%errOut,*) 'Sideneighbor not found but localneighborside non zero!'
!             STOP
!          ENDIF
!        ENDIF
!      ENDDO
!    ENDDO

    logInfo(*) 'Mesh seems to be consistent.'
    !
    ! Compute the elements' volumes and barycenters
    !
    logInfo(*) 'Computing element volumes and barycenters.'
    !
    DO iElem = 1, MESH%nELEM
        SELECT CASE(MESH%LocalElemType(iElem))
        CASE(4)
            x1 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1,iElem))
            y1 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1,iElem))
            z1 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1,iElem))
            !
            x2 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(2,iElem))
            y2 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(2,iElem))
            z2 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(2,iElem))
            !
            x3 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(3,iElem))
            y3 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(3,iElem))
            z3 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(3,iElem))
            !
            x4 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(4,iElem))
            y4 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(4,iElem))
            z4 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(4,iElem))
            !
            ! Compute with SARRUS's law from the Jacobi determinant of the transformation
            ! from xyz coordinates to the reference tetrahedron
            !
            MESH%ELEM%Volume(iElem) = abs(  (x2-x1)*(y3-y1)*(z4-z1) +             &
                                            (x3-x1)*(y4-y1)*(z2-z1) +             &
                                            (x4-x1)*(y2-y1)*(z3-z1) -             &
                                            (x4-x1)*(y3-y1)*(z2-z1) -             &
                                            (x2-x1)*(y4-y1)*(z3-z1) -             &
                                            (x3-x1)*(y2-y1)*(z4-z1)         )/6.0D0
            ! Barycenters of the tetrahedrons
            MESH%ELEM%xyBary(:,iElem) = 0.25*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem)) + &
                                          MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem))   )
            !
            ! Barycenters of the four sides
            !SideBary(:,1) = 1./3.*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem))   )
            !
            !SideBary(:,2) = 1./3.*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem))   )
            !
            !SideBary(:,3) = 1./3.*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem))   )
            !
            !SideBary(:,4) = 1./3.*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem)) + &
            !                         MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem))   )
            !
        CASE(6)
            x1 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1,iElem))
            y1 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1,iElem))
            z1 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1,iElem))
            !
            x2 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(2,iElem))
            y2 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(2,iElem))
            z2 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(2,iElem))
            !
            x3 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(3,iElem))
            y3 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(3,iElem))
            z3 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(3,iElem))
            !
            x4 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(4,iElem))
            y4 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(4,iElem))
            z4 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(4,iElem))
            !
            x5 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(5,iElem))
            y5 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(5,iElem))
            z5 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(5,iElem))
            !
            x6 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(6,iElem))
            y6 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(6,iElem))
            z6 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(6,iElem))
            !
            x7 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(7,iElem))
            y7 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(7,iElem))
            z7 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(7,iElem))
            !
            x8 = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(8,iElem))
            y8 = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(8,iElem))
            z8 = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(8,iElem))
            !
            ! Compute the hexahedron's volume by the sum of the 5 tetrahedrons,
            ! that subdivide it. The tetrahedron's volume is computed with SARRUS's law
            ! as before
            !
            MESH%ELEM%Volume(iElem) = 0.
            DO m = 1,5
                SELECT CASE (m)
                CASE(1)
                    xt1 = x5; xt2 = x7; xt3 = x6; xt4 = x1;
                    yt1 = y5; yt2 = y7; yt3 = y6; yt4 = y1;
                    zt1 = z5; zt2 = z7; zt3 = z6; zt4 = z1;
                CASE(2)
                    xt1 = x7; xt2 = x3; xt3 = x4; xt4 = x1;
                    yt1 = y7; yt2 = y3; yt3 = y4; yt4 = y1;
                    zt1 = z7; zt2 = z3; zt3 = z4; zt4 = z1;
                CASE(3)
                    xt1 = x6; xt2 = x8; xt3 = x4; xt4 = x7;
                    yt1 = y6; yt2 = y8; yt3 = y4; yt4 = y7;
                    zt1 = z6; zt2 = z8; zt3 = z4; zt4 = z7;
                CASE(4)
                    xt1 = x6; xt2 = x4; xt3 = x2; xt4 = x1;
                    yt1 = y6; yt2 = y4; yt3 = y2; yt4 = y1;
                    zt1 = z6; zt2 = z4; zt3 = z2; zt4 = z1;
                CASE(5)
                    xt1 = x6; xt2 = x1; xt3 = x7; xt4 = x4;
                    yt1 = y6; yt2 = y1; yt3 = y7; yt4 = y4;
                    zt1 = z6; zt2 = z1; zt3 = z7; zt4 = z4;
                END SELECT

                MESH%ELEM%Volume(iElem) = MESH%ELEM%Volume(iElem)           +    &
                                        abs(  (xt2-xt1)*(yt3-yt1)*(zt4-zt1) +    &
                                              (xt3-xt1)*(yt4-yt1)*(zt2-zt1) +    &
                                              (xt4-xt1)*(yt2-yt1)*(zt3-zt1) -    &
                                              (xt4-xt1)*(yt3-yt1)*(zt2-zt1) -    &
                                              (xt2-xt1)*(yt4-yt1)*(zt3-zt1) -    &
                                              (xt3-xt1)*(yt2-yt1)*(zt4-zt1)      )/6.0D0

            ENDDO
            ! Barycenters of the hexahedrons
            MESH%ELEM%xyBary(:,iElem) = 1./8.*(  MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(1,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(2,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(3,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(4,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(5,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(6,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(7,iElem)) + &
                                           MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(8,iElem))   )


        END SELECT
    ENDDO
    !
    ! Build stencil if necessary
    !
    ! Generate periodic information, if necessary
    !

    DO iElem = 1, MESH%nElem
       eType = MESH%LocalElemType(iElem)
       DO iSide = 1, eType
          IF(MESH%ELEM%Reference(iSide,iElem).EQ.106)THEN
             exit_loop = .true.
             EXIT
          ENDIF
       ENDDO !iSide
       IF (exit_loop) EXIT
    ENDDO ! iElem
    ! Based on eType we call the corresponding periodic boundary search
    SELECT CASE(eType)
    CASE(4)
        CALL BND_DG_Periodic3D_Tetra_us(MESH,BND,IO)
    CASE(6)
        CALL BND_DG_Periodic3D_Hexa_us(MESH,BND,IO)
    END SELECT

    logInfo(*) 'Start neighbour search for the non-conforming elements.'
    

    !----------------------------------------------------------------------
    !Start neighbour search for the non-conforming elements in the mesh
    !----------------------------------------------------------------------
    ! Gaussian Quadrature points on each side of a possibly non-conforming
    ! boundary element are computed and the corresponding element that includes
    ! the Gaussian Quadrature point is identified as neighbor.

     ALLOCATE(MESH%ELEM%GP_Tet(3,(DISC%Galerkin%nPoly+2)**3),                             & !could be local
             MESH%ELEM%GW_Tet((DISC%Galerkin%nPoly+2)**3),                                & !could be local
             MESH%ELEM%GP_Hex(3,(DISC%Galerkin%nPoly+2)**3),                              & !could be local
             MESH%ELEM%GW_Hex((DISC%Galerkin%nPoly+2)**3),                                & !could be local
             MESH%ELEM%BndGP_Tri(2,(DISC%Galerkin%nPoly+2)**2),                           & !could be local?
             MESH%ELEM%BndGW_Tri((DISC%Galerkin%nPoly+2)**2),                             &
             MESH%ELEM%BndGP_Quad(2,(DISC%Galerkin%nPoly+2)**2),                          & !could be local?
             MESH%ELEM%BndGW_Quad((DISC%Galerkin%nPoly+2)**2),                            &
             MESH%ELEM%BndGP_Hex(EQN%Dimension,6,(DISC%Galerkin%nPoly+2)**2),             & !could be local?
             MESH%ELEM%BndGP_Tet(EQN%Dimension,4,(DISC%Galerkin%nPoly+2)**2)              ) !could be local?

    ! Allocations regarding non-conforming meshes:
   IF (MESH%nNonConformingEdges.GT.0) THEN
       ALLOCATE(MESH%ELEM%ncBndNeighbor(MESH%nNonConformingEdges,(DISC%Galerkin%nPoly+2)**2),             &
                MESH%ELEM%ncBndGaussP(MESH%nNonConformingEdges,(DISC%Galerkin%nPoly+2)**2,EQN%Dimension), &
	            MESH%ELEM%NCB_IndexList(MESH%nNonConformingEdges,2))
       ALLOCATE(MESH%ELEM%ncIndex(MESH%nElem,MESH%nSideMax))
       MESH%ELEM%ncIndex = 0
   ENDIF


    ! Quadrature points depend on the basis function used
    CALL TetrahedronQuadraturePoints(                     &
                   nIntGP     = DISC%Galerkin%nIntGP,     &
                   IntGaussP  = MESH%ELEM%GP_Tet,         &
                   IntGaussW  = MESH%ELEM%GW_Tet,         &
                   M          = DISC%Galerkin%nPoly+2,    &
                   IO         = IO,                       &
                   quiet      = .TRUE.                    )


    CALL HypercubeQuadraturePoints(                       &
                   nIntGP     = DISC%Galerkin%nIntGP,     &
                   IntGaussP  = MESH%ELEM%GP_Hex,         &
                   IntGaussW  = MESH%ELEM%GW_Hex,         &
                   M          = DISC%Galerkin%nPoly+2,    &
                   nDim       = 3,                        &
                   X1         = (/ 0., 0., 0. /),         &
                   X2         = (/ 1., 1., 1. /),         &
                   IO         = IO,                       &
                   quiet      = .TRUE.                    )

    nIntGP     = DISC%Galerkin%nIntGP


    CALL TriangleQuadraturePoints(                        &
                  nIntGP     = DISC%Galerkin%nBndGP,      &
                  IntGaussP  = MESH%ELEM%BndGP_Tri,       &
                  IntGaussW  = MESH%ELEM%BndGW_Tri,       &
                  M          = DISC%Galerkin%nPoly+2,     &
                  IO         = IO,                        &
                  quiet      = .TRUE.                     )

    CALL HypercubeQuadraturePoints(                       &
                  nIntGP     = DISC%Galerkin%nBndGP,      &
                  IntGaussP  = MESH%ELEM%BndGP_Quad,      &
                  IntGaussW  = MESH%ELEM%BndGW_Quad,      &
                  M          = DISC%Galerkin%nPoly+2,     &
                  nDim       = 2,                         &
                  X1         = (/ 0., 0. /),              &
                  X2         = (/ 1., 1. /),              &
                  IO         = IO,                        &
                  quiet      = .TRUE.                     )

    nBndGP     = DISC%Galerkin%nBndGP

    ! Assign quadrature points on each side
    DO iSide = 1, 4
      DO iBndGP = 1, DISC%Galerkin%nBndGP
        chiGP  = MESH%ELEM%BndGP_Tri(1,iBndGP)
        tauGP  = MESH%ELEM%BndGP_Tri(2,iBndGP)
        CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chiGP,tauGP,iSide,0)
        !---------------------------------------------
        MESH%ELEM%BndGP_Tet(1,iSide,iBndGP) = xi
        MESH%ELEM%BndGP_Tet(2,iSide,iBndGP) = eta
        MESH%ELEM%BndGP_Tet(3,iSide,iBndGP) = zeta
        !---------------------------------------------
      ENDDO
    ENDDO

    DO iSide = 1, 6
      DO iBndGP = 1, DISC%Galerkin%nBndGP
        chiGP  = MESH%ELEM%BndGP_Quad(1,iBndGP)
        tauGP  = MESH%ELEM%BndGP_Quad(2,iBndGP)
        CALL HexaTrafoChiTau2XiEtaZeta(xi,eta,zeta,chiGP,tauGP,iSide,0)
        !---------------------------------------------
        MESH%ELEM%BndGP_Hex(1,iSide,iBndGP) = xi
        MESH%ELEM%BndGP_Hex(2,iSide,iBndGP) = eta
        MESH%ELEM%BndGP_Hex(3,iSide,iBndGP) = zeta
        !---------------------------------------------
      ENDDO
    ENDDO

    ! Generate additional list for nc boundaries where indices of original list (iElem, iSide) are stored.
    counter = 0
    IF (MESH%nNonConformingEdges.GT.0) THEN
        MESH%ELEM%NCB_IndexList(:,:) = 0
        DO iElem = 1, MESH%nElem
          DO iSide = 1, MESH%LocalElemType(iElem)
            IF(MESH%ELEM%Reference(iSide,iElem).EQ.102)THEN
              counter = counter + 1
              MESH%ELEM%NCB_IndexList(counter,1) = iElem
              MESH%ELEM%NCB_IndexList(counter,2) = iSide
            ENDIF
          ENDDO
        ENDDO
    ENDIF

    IF(MESH%nNonConformingEdges .NE. counter)THEN
         logError(*) 'Non conforming boundaries do not match!'
         STOP
    ENDIF

    ALLOCATE( DISC%Galerkin%geoNormals( 3,MESH%GlobalElemType,MESH%nElem),  &
              DISC%Galerkin%geoSurfaces(  MESH%GlobalElemType,MESH%nElem))
    DISC%Galerkin%geoNormals = 0.
    DISC%Galerkin%geoSurfaces = 0.

    DO iNCboundary = 1, MESH%nNonConformingEdges
         iElem = MESH%ELEM%NCB_IndexList(iNCboundary,1)
         iSide = MESH%ELEM%NCB_IndexList(iNCboundary,2)
         eType = MESH%LocalElemType(iElem)
         nVrtx = MESH%LocalVrtxType(iElem)

         !Search neighbor element that contains quadrature point on the BOUNDARY

         x(1:nVrtx) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:eType,iElem))
         y(1:nVrtx) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:eType,iElem))
         z(1:nVrtx) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:eType,iElem))

         DO iBndGP = 1,DISC%Galerkin%nBndGP
            SELECT CASE(eType)
            CASE(4)
               xi   = MESH%ELEM%BndGP_Tet(1,iSide,iBndGP)
               eta  = MESH%ELEM%BndGP_Tet(2,iSide,iBndGP)
               zeta = MESH%ELEM%BndGP_Tet(3,iSide,iBndGP)
               CALL TetraTrafoXiEtaZeta2XYZ(p_new(1),p_new(2),p_new(3),xi,eta,zeta,x(1:nVrtx),y(1:nVrtx),z(1:nVrtx))
            CASE(6)
               xi   = MESH%ELEM%BndGP_Hex(1,iSide,iBndGP)
               eta  = MESH%ELEM%BndGP_Hex(2,iSide,iBndGP)
               zeta = MESH%ELEM%BndGP_Hex(3,iSide,iBndGP)
               CALL HexaTrafoXiEtaZeta2XYZ(p_new(1),p_new(2),p_new(3),xi,eta,zeta,x(1:nVrtx),y(1:nVrtx),z(1:nVrtx))
            END SELECT

            mirrored_point_in_gap = .TRUE.

            ! We search which other element contains the integration point p_new
            DO i = 1,MESH%nNonConformingEdges
               j = MESH%ELEM%NCB_IndexList(i,1) ! corresponds to iElem
               sidevec(:,:) = 0.

                ! Compute normal vector needed in XYZInElement for hexas
                ! In the future subroutine XYZInElement could be independent of
                ! variable DISC%Galerkin%geoNormals
                IF(MESH%LocalElemType(j).EQ.6) THEN
                    DO k = 1, MESH%LocalElemType(j)

                        ! Boundary side vector pointing in chi-direction
                        sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,2),j)) - &
                                       MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,1),j))
                        ! Boundary side vector pointing in tau-direction
                        sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,3),j)) - &
                                       MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,1),j))
                        ! Normal vector computed by cross product
                        DISC%Galerkin%geoNormals(1,k,j) = sidevec(2,1)*sidevec(3,2) - sidevec(3,1)*sidevec(2,2)
                        DISC%Galerkin%geoNormals(2,k,j) = sidevec(3,1)*sidevec(1,2) - sidevec(1,1)*sidevec(3,2)
                        DISC%Galerkin%geoNormals(3,k,j) = sidevec(1,1)*sidevec(2,2) - sidevec(2,1)*sidevec(1,2)
                        ! Triangle's surface = 0.5 * cross_product
                        DISC%Galerkin%geoSurfaces(k,j)  =  0.5*SQRT(                             &
                                                           DISC%Galerkin%geoNormals(1,k,j)**2 +  &
                                                           DISC%Galerkin%geoNormals(2,k,j)**2 +  &
                                                           DISC%Galerkin%geoNormals(3,k,j)**2    )
                        ! Boundary side vector pointing in chi-direction
                        sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,4),j)) - &
                                       MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,3),j))
                        ! Boundary side vector pointing in tau-direction
                        sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,2),j)) - &
                                       MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(k,3),j))
                        ! Normal vector computed by cross product
                        DISC%Galerkin%geoNormals(1,k,j) = sidevec(2,1)*sidevec(3,2) - sidevec(3,1)*sidevec(2,2)
                        DISC%Galerkin%geoNormals(2,k,j) = sidevec(3,1)*sidevec(1,2) - sidevec(1,1)*sidevec(3,2)
                        DISC%Galerkin%geoNormals(3,k,j) = sidevec(1,1)*sidevec(2,2) - sidevec(2,1)*sidevec(1,2)
                        ! Second triangle's surface = 0.5 * cross_product
                        DISC%Galerkin%geoSurfaces(k,j)  =  DISC%Galerkin%geoSurfaces(k,j) +      &
                                                           0.5*SQRT(                             &
                                                           DISC%Galerkin%geoNormals(1,k,j)**2 +  &
                                                           DISC%Galerkin%geoNormals(2,k,j)**2 +  &
                                                           DISC%Galerkin%geoNormals(3,k,j)**2    )
                        ! Normalize normal vector to length 1
                        DISC%Galerkin%geoNormals(:,k,j) =  DISC%Galerkin%geoNormals(:,k,j) / &
                                                           SQRT(SUM(DISC%Galerkin%geoNormals(:,k,j)**2))
                    ENDDO
                ENDIF

               ! Ask if p_new belongs to element j
               IF(XYZInElement(p_new(1),p_new(2),p_new(3),j,1e-5,MESH,DISC)) THEN
               IF(j.NE.iElem)THEN
                  MESH%ELEM%ncBndNeighbor(iNCboundary,iBndGP) = i
                  MESH%ELEM%ncBndGaussP(iNCboundary,iBndGP,:) = p_new(:)
                  mirrored_point_in_gap = .FALSE.
                  EXIT
               ENDIF
               ENDIF
             ENDDO

             IF(mirrored_point_in_gap)THEN
                logError(*) 'Mirrored Gauss point at non-conforming boundary fell in a gap!'
                logError(*) 'Check element ', iElem,' on Side ',iSide
                STOP
             ENDIF
             MESH%ELEM%ncIndex(iElem,iSide) = iNCboundary
          ENDDO
    ENDDO

    ! Building fault information for dynamic rupture
    IF(EQN%DR.EQ.1) THEN
      !
      logInfo(*) 'Start fault element search for Dynamic Rupture.'
      !
      count = 0
      !
      ALLOCATE(tmp_el(MESH%nElem,2))
      !
      tmp_el(:,:) = 0
      !
      ! Find all 103 side references and count how many
      DO iElem = 1, MESH%nElem
        eType = MESH%LocalElemType(iElem)
        DO iSide = 1, eType
          IF(MESH%ELEM%Reference(iSide,iElem).EQ.103) THEN
            count = count+1
            tmp_el(count,1) = iElem
            tmp_el(count,2) = iSide
          ENDIF
        ENDDO
      ENDDO
      !
      IF(TmpFaultElements.NE.count) THEN
         logError(*) 'Number of fault segments doesn.t agree with number of'
         logError(*) 'elements with 103 condition!'
         STOP
      ENDIF
      !
      MESH%Fault%nElem = count
      MESH%Fault%nSide = count / 2
      !
      !
      ALLOCATE(MESH%Fault%Face(MESH%Fault%nSide,2,2))
      MESH%Fault%Face(:,:,:) = 0

      ! create Face list:
      !
      ! store values of tmp_el in Face and do this also for neighbor (cycle out double)
      count = 1
      k = 1
      iElem = tmp_el(k,1)
      i     = tmp_el(k,2)
      MESH%Fault%Face(count,1,1) = iElem
      MESH%Fault%Face(count,2,1) = i
      MESH%Fault%Face(count,1,2) = MESH%ELEM%SideNeighbor(i,iElem)
      MESH%Fault%Face(count,2,2) = MESH%ELEM%LocalNeighborSide(i,iElem)
      DO k = 2, MESH%Fault%nElem
          cycle = .false.
          iElem = tmp_el(k,1)
          i     = tmp_el(k,2)
          DO j = 1,MESH%Fault%nSide
              IF (iElem == MESH%Fault%Face(j,1,2) .AND. i == MESH%Fault%Face(j,2,2)) THEN
                  cycle = .true.
              ENDIF
          ENDDO
          IF (cycle) CYCLE
          count = count + 1
          MESH%Fault%Face(count,1,1) = iElem
          MESH%Fault%Face(count,2,1) = i
          MESH%Fault%Face(count,1,2) = MESH%ELEM%SideNeighbor(i,iElem)
          MESH%Fault%Face(count,2,2) = MESH%ELEM%LocalNeighborSide(i,iElem)
      ENDDO
      !
      IF (count .NE. MESH%Fault%nSide) THEN
          logError(*) 'Problem in cycle out double fault members'
          STOP
      ENDIF
      !
      ! put elements in order: plus - minus side
      XRef =  EQN%XRef
      YRef =  EQN%YRef
      ZRef =  EQN%ZRef
      DO i=1,MESH%Fault%nSide
        iElem  = MESH%Fault%Face(i,1,1)
        iSide  = MESH%Fault%Face(i,2,1)
        iElem2 = MESH%Fault%Face(i,1,2)
        iSide2 = MESH%Fault%Face(i,2,2)

        Bary1(:)=MESH%ELEM%xyBary(:,MESH%Fault%Face(i,1,1))
        Bary2(:)=MESH%ELEM%xyBary(:,MESH%Fault%Face(i,1,2))
        Dis1=SQRT((Bary1(1)-XRef)**2+(Bary1(2)-YRef)**2+(Bary1(3)-ZRef)**2)
        Dis2=SQRT((Bary2(1)-XRef)**2+(Bary2(2)-YRef)**2+(Bary2(3)-ZRef)**2)
        IF(Dis1.GE.Dis2) THEN
          MESH%Fault%Face(i,1,1) = iElem2
          MESH%Fault%Face(i,2,1) = iSide2
          MESH%Fault%Face(i,1,2) = iElem
          MESH%Fault%Face(i,2,2) = iSide
        ENDIF
      ENDDO
      !
      DEALLOCATE(tmp_el)
      !
      DISC%DynRup%DR_output = .TRUE.              ! output "on" for this domain
      MESH%Fault%nPlusSide  = MESH%Fault%nSide
      !
      ! Store normal vector information for iElem = 0 since it is used in friction.f90
      ! (in the first step for all elements since the information to be stored is small)
      ALLOCATE(MESH%Fault%geoNormals(3,MESH%Fault%nSide)  , &
               MESH%Fault%geoTangent1(3,MESH%Fault%nSide) , &
               MESH%Fault%geoTangent2(3,MESH%Fault%nSide)   )
      DO i = 1,MESH%Fault%nSide
         iElem = MESH%Fault%Face(i,1,1)
         iSide = MESH%Fault%Face(i,2,1)

         ! Boundary side vector pointing in chi-direction
         sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,2),iElem)) - &
                        MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
         ! Boundary side vector pointing in tau-direction
         sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,3),iElem)) - &
                        MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
         ! Normal vector computed by cross product
         MESH%Fault%geoNormals(:,i) =  sidevec(:,1).x.sidevec(:,2)
         ! Normalize normal vector to length 1
         Length  =  SQRT(MESH%Fault%geoNormals(1,i)**2 +  &
                         MESH%Fault%geoNormals(2,i)**2 +  &
                         MESH%Fault%geoNormals(3,i)**2    )
         MESH%Fault%geoNormals(:,i) =  MESH%Fault%geoNormals(:,i) / Length

         ! Compute vector inside the triangle's plane for the rotation matrix
         MESH%Fault%geoTangent1(:,i) = sidevec(:,1)
         ! Normalize to 1
         Length = SQRT(MESH%Fault%geoTangent1(1,i)**2 + &
                       MESH%Fault%geoTangent1(2,i)**2 + &
                       MESH%Fault%geoTangent1(3,i)**2   )
         MESH%Fault%geoTangent1(:,i) = MESH%Fault%geoTangent1(:,i)/Length
         ! Compute second vector in the plane, orthogonal to the normal and tangent 1 vectors
         ! using the crossproduct
         MESH%Fault%geoTangent2(:,i) = MESH%Fault%geoNormals(:,i) .x. &
                                       MESH%Fault%geoTangent1(:,i)
      ENDDO ! i = 1,MESH%Fault%nSide
      !
    ENDIF ! Dynamic Rupture
    
    !
    DEALLOCATE(DISC%Galerkin%geoNormals)
    DEALLOCATE(DISC%Galerkin%geoSurfaces)
    !

  END SUBROUTINE compute_mesh_Gambit3D_Mixed

  SUBROUTINE BND_DG_Periodic3D_Tetra_us(MESH,BND,IO)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tBoundary)          :: BND
    TYPE (tInputOutput)       :: IO
    !--------------------------------------------------------------------------
    INTEGER                   :: i, j, k, iDir
    INTEGER                   :: Dir1, Dir2
    INTEGER                   :: iElem, iNeighbor, iSide, iNeighborSide
    INTEGER                   :: counter
    INTEGER                   :: nPeriodic
    INTEGER                   :: VertexSide(4,3)
    INTEGER                   :: nxmin(3), nxmax(3), xmatch(3)
    INTEGER                   :: Directions(3,2)
    INTEGER, POINTER          :: ElemSide(:,:)
    REAL, POINTER             :: MidSide(:,:)
    REAL                      :: xmin(3),xmax(3)
    REAL                      :: Vertex1(3), Vertex2(3)
    REAL                      :: x1(3), x2(3)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: BND, IO
    INTENT(INOUT)             :: MESH
    !--------------------------------------------------------------------------
    IF(BND%periodic.EQ.0) RETURN

    logInfo(*) 'Generating periodic information for tetrahedral mesh....'

    VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
    VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
    VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
    VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

    Directions(1,:) = (/ 2, 3 /)       ! Coordinate numbers for x-periodicity (yz)
    Directions(2,:) = (/ 1, 3 /)       ! Coordinate numbers for y-periodicity (xz)
    Directions(3,:) = (/ 1, 2 /)       ! Coordinate numbers for z-periodicity (xy)

    nPeriodic = 0

    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%LocalElemType(iElem)
        IF(MOD(MESH%ELEM%Reference(iSide,iElem),100).EQ.6) THEN
          nPeriodic = nPeriodic + 1
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE(MidSide(nPeriodic,3), ElemSide(nPeriodic,2))

    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1,  MESH%LocalElemType(iElem)
        IF(MOD(MESH%ELEM%Reference(iSide,iElem),100).EQ.6) THEN
          counter = counter + 1
          MidSide(counter,:) = 1./3.*(MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem)) + &
                                      MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,2),iElem)) + &
                                      MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,3),iElem))   )
          ElemSide(counter,:) = (/ iElem, iSide /)
        ENDIF
      ENDDO
    ENDDO

    DO iDir = 1, 3
        xmin(iDir) = MINVAL(MidSide(:,iDir))
        xmax(iDir) = MAXVAL(MidSide(:,iDir))
        IF( xmin(iDir).CA.xmax(iDir) ) THEN
           logError(*) 'Error in periodic BND. Singular domain: ', iDir, xmin(iDir), xmax(iDir)
           !STOP
        ENDIF
    ENDDO

    nxmin(:) = 0
    nxmax(:) = 0

    DO counter = 1, nPeriodic
        DO iDir = 1, 3
            IF(MidSide(counter,iDir).CA.xmin(iDir))   nxmin(iDir) = nxmin(iDir) + 1
            IF(MidSide(counter,iDir).CA.xmax(iDir))   nxmax(iDir) = nxmax(iDir) + 1
        ENDDO
    ENDDO

    xmatch(:) = 0

    ! Periodicity in direction iDir (x=1,y=2,z=3)
    DO iDir = 1, 3
        Dir1 = Directions(iDir,1)
        Dir2 = Directions(iDir,2)
        IF(BND%DirPeriodic(iDir)) THEN
           IF(nxmin(iDir).NE.nxmax(iDir)) THEN
              logError(*) 'Number of periodic BND sides does not match. ', iDir,nxmin(iDir),nxmax(iDir)
              STOP
           ENDIF
           DO i = 1, nPeriodic
             iElem          = ElemSide(i,1)
             iSide          = ElemSide(i,2)
             x1(:)          = MidSide(i,:)
             IF( (x1(iDir).CA.xmin(iDir)) .OR. (x1(iDir).CA.xmax(iDir)) ) THEN
                 MESH%ELEM%SideNeighbor(iSide,iElem)      = -1
                 MESH%ELEM%LocalNeighborSide(iSide,iElem) = -1
                 DO j = 1, nPeriodic
                     IF(i.EQ.j) THEN
                        CYCLE
                     ENDIF
                     iNeighbor      = ElemSide(j,1)
                     iNeighborSide  = ElemSide(j,2)
                     x2(:)          = MidSide(j,:)
                     ! Boundary triangles match
                     IF( (x1(Dir1).CA.x2(Dir1)).AND.(x1(Dir2).CA.x2(Dir2)) ) THEN
                       xmatch(iDir) = xmatch(iDir) + 1
                       MESH%ELEM%SideNeighbor(iSide,iElem)      = iNeighbor
                       MESH%ELEM%LocalNeighborSide(iSide,iElem) = iNeighborSide
                       Vertex1(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem))
                       MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = -1
                       DO k = 1,3
                         Vertex2(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iNeighborSide,k),iNeighbor))
                         IF( (Vertex2(Dir1).CA.Vertex1(Dir1)) .AND. (Vertex2(Dir2).CA.Vertex1(Dir2)) ) THEN
                           MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = k
                           EXIT
                         ENDIF
                       ENDDO
                       IF(MESH%ELEM%LocalNeighborVrtx(iSide,iElem).EQ.-1) THEN
                          logError(*) 'No local neighbor vertex found. '
                          logError(*) iElem,iSide,iNeighbor,iNeighborSide
                          logError(*) x1(:), x2(:), Vertex1(:)
                          STOP
                       ENDIF
                     ENDIF
                 ENDDO
                 IF(MESH%ELEM%SideNeighbor(iSide,iElem).EQ.-1) THEN
                      logError(*) 'No periodic sideneighbor found. '
                      logError(*) iElem,iSide
                      logError(*) x1(:)
                      STOP
                 ENDIF
              ENDIF
           ENDDO
           IF( xmatch(iDir).NE.(nxmin(iDir)+nxmax(iDir)) ) THEN
              logError(*) 'Not all boundary surfaces cound be matched for direction ', iDir
              logError(*) nxmin(iDir), nxmax(iDir), xmatch(iDir)
              STOP
           ENDIF
        ENDIF
    ENDDO

    DEALLOCATE(MidSide, ElemSide)

  END SUBROUTINE BND_DG_Periodic3D_Tetra_us


  SUBROUTINE BND_DG_Periodic3D_Hexa_us(MESH,BND,IO)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tBoundary)          :: BND
    TYPE (tInputOutput)       :: IO
    !--------------------------------------------------------------------------
    INTEGER                   :: i, j, k, iDir
    INTEGER                   :: Dir1, Dir2
    INTEGER                   :: iElem, iNeighbor, iSide, iNeighborSide
    INTEGER                   :: counter
    INTEGER                   :: nPeriodic
    INTEGER                   :: VertexSide(6,4)
    INTEGER                   :: nxmin(3), nxmax(3), xmatch(3)
    INTEGER                   :: Directions(3,2)
    INTEGER, POINTER          :: ElemSide(:,:)
    REAL, POINTER             :: MidSide(:,:)
    REAL                      :: xmin(3),xmax(3)
    REAL                      :: Vertex1(3), Vertex2(3)
    REAL                      :: x1(3), x2(3)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: BND, IO
    INTENT(INOUT)             :: MESH
    !--------------------------------------------------------------------------
    IF(BND%periodic.EQ.0) RETURN

    logInfo(*) 'Generating periodic information for hexahedral mesh....'

    VertexSide(1,:) =  (/ 1, 2, 6, 5 /)   ! Local hex. vertices of hex. side I   !
    VertexSide(2,:) =  (/ 2, 4, 8, 6 /)   ! Local hex. vertices of hex. side II  !
    VertexSide(3,:) =  (/ 4, 3, 7, 8 /)   ! Local hex. vertices of hex. side III !
    VertexSide(4,:) =  (/ 3, 1, 5, 7 /)   ! Local hex. vertices of hex. side IV  !
    VertexSide(5,:) =  (/ 2, 1, 3, 4 /)   ! Local hex. vertices of hex. side V   !
    VertexSide(6,:) =  (/ 5, 6, 8, 7 /)   ! Local hex. vertices of hex. side VI  !

    Directions(1,:) = (/ 2, 3 /)          ! Coordinate numbers for x-periodicity (yz)
    Directions(2,:) = (/ 1, 3 /)          ! Coordinate numbers for y-periodicity (xz)
    Directions(3,:) = (/ 1, 2 /)          ! Coordinate numbers for z-periodicity (xy)

    nPeriodic = 0

    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%LocalElemType(iElem)
        IF(MOD(MESH%ELEM%Reference(iSide,iElem),100).EQ.6) THEN
          nPeriodic = nPeriodic + 1
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE(MidSide(nPeriodic,3), ElemSide(nPeriodic,2))

    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%LocalElemType(iElem)
        IF(MOD(MESH%ELEM%Reference(iSide,iElem),100).EQ.6) THEN
          counter = counter + 1
          MidSide(counter,:) = 1./4.*(MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem)) + &
                                      MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,2),iElem)) + &
                                      MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,3),iElem)) + &
                                      MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,4),iElem))   )
          ElemSide(counter,:) = (/ iElem, iSide /)
        ENDIF
      ENDDO
    ENDDO

    DO iDir = 1, 3
        xmin(iDir) = MINVAL(MidSide(:,iDir))
        xmax(iDir) = MAXVAL(MidSide(:,iDir))
        IF( xmin(iDir).CA.xmax(iDir) ) THEN
           logError(*) 'Error in periodic BND. Singular domain: ', iDir, xmin(iDir), xmax(iDir)
        ENDIF
    ENDDO

    nxmin(:) = 0
    nxmax(:) = 0

    DO counter = 1, nPeriodic
        DO iDir = 1, 3
            IF(MidSide(counter,iDir).CA.xmin(iDir))   nxmin(iDir) = nxmin(iDir) + 1
            IF(MidSide(counter,iDir).CA.xmax(iDir))   nxmax(iDir) = nxmax(iDir) + 1
        ENDDO
    ENDDO

    xmatch(:) = 0

    ! Periodicity in direction iDir (x=1,y=2,z=3)
    DO iDir = 1, 3
        Dir1 = Directions(iDir,1)
        Dir2 = Directions(iDir,2)
        IF(BND%DirPeriodic(iDir)) THEN
           IF(nxmin(iDir).NE.nxmax(iDir)) THEN
              logError(*) 'Number of periodic BND sides does not match. ', iDir,nxmin(iDir),nxmax(iDir)
              STOP
           ENDIF
           DO i = 1, nPeriodic
             iElem          = ElemSide(i,1)
             iSide          = ElemSide(i,2)
             x1(:)          = MidSide(i,:)
             IF( (x1(iDir).CA.xmin(iDir)) .OR. (x1(iDir).CA.xmax(iDir)) ) THEN
                 MESH%ELEM%SideNeighbor(iSide,iElem)      = -1
                 MESH%ELEM%LocalNeighborSide(iSide,iElem) = -1
                 DO j = 1, nPeriodic
                     IF(i.EQ.j) THEN
                        CYCLE
                     ENDIF
                     iNeighbor      = ElemSide(j,1)
                     iNeighborSide  = ElemSide(j,2)
                     x2(:)          = MidSide(j,:)
                     ! Boundary quadrilaterals match
                     IF( (x1(Dir1).CA.x2(Dir1)).AND.(x1(Dir2).CA.x2(Dir2)) ) THEN
                       xmatch(iDir) = xmatch(iDir) + 1
                       MESH%ELEM%SideNeighbor(iSide,iElem)      = iNeighbor
                       MESH%ELEM%LocalNeighborSide(iSide,iElem) = iNeighborSide
                       Vertex1(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem))
                       MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = -1
                       DO k = 1,4
                         Vertex2(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iNeighborSide,k),iNeighbor))
                         IF( (Vertex2(Dir1).CA.Vertex1(Dir1)) .AND. (Vertex2(Dir2).CA.Vertex1(Dir2)) ) THEN
                           MESH%ELEM%LocalNeighborVrtx(iSide,iElem) = k
                           EXIT
                         ENDIF
                       ENDDO
                       IF(MESH%ELEM%LocalNeighborVrtx(iSide,iElem).EQ.-1) THEN
                          logError(*) 'No local neighbor vertex found. '
                          logError(*) iElem,iSide,iNeighbor,iNeighborSide
                          logError(*) x1(:), x2(:), Vertex1(:)
                          STOP
                       ENDIF
                     ENDIF
                 ENDDO
                 IF(MESH%ELEM%SideNeighbor(iSide,iElem).EQ.-1) THEN
                      logError(*) 'No periodic sideneighbor found. '
                      logError(*) iElem,iSide
                      logError(*) x1(:)
                      STOP
                 ENDIF
              ENDIF
           ENDDO
           IF( xmatch(iDir).NE.(nxmin(iDir)+nxmax(iDir)) ) THEN
              logError(*) 'Not all boundary surfaces cound be matched for direction ', iDir
              logError(*) nxmin(iDir), nxmax(iDir), xmatch(iDir)
              STOP
           ENDIF
        ENDIF
    ENDDO

    DEALLOCATE(MidSide, ElemSide)

  END SUBROUTINE BND_DG_Periodic3D_Hexa_us


END MODULE compute_mesh_info_mod
