!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2010, SeisSol Group
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

#ifdef BG 
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE read_mesh_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  USE allocate_mesh_mod, ONLY : allocate_mesh_level0_1,allocate_mesh_level0_2
  USE compute_mesh_info_mod

  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE read_mesh
     MODULE PROCEDURE read_mesh
  END INTERFACE

  INTERFACE read_mesh_Gambit3D_Tetra
     MODULE PROCEDURE read_mesh_Gambit3D_Tetra
  END INTERFACE
  INTERFACE read_mesh_ICEMCFD3D_Tetra
     MODULE PROCEDURE read_mesh_ICEMCFD3D_Tetra
  END INTERFACE  
  INTERFACE read_mesh_Gambit3D_Hexa
     MODULE PROCEDURE read_mesh_Gambit3D_Hexa
  END INTERFACE

  INTERFACE read_mesh_Gambit3D_Mixed
    MODULE PROCEDURE read_mesh_Gambit3D_Mixed
  END INTERFACE

  INTERFACE ObjectReferences_Side
     MODULE PROCEDURE ObjectReferences_Side
  END INTERFACE

  TYPE tTetraNodeList
        TYPE(tTetraNodeList), POINTER :: next, prev
        INTEGER                  :: Vertex(4) 
  END TYPE tTetraNodeList 

  !----------------------------------------------------------------------------
  PUBLIC  :: read_mesh
  PUBLIC  :: ObjectReferences_Side
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE read_mesh(IO,EQN,DISC,MESH,BND,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tInputOutput)       :: IO
    TYPE (tEquations)         :: EQN
    TYPE (tDiscretization)    :: DISC
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tBoundary)          :: BND
    TYPE (tMPI)               :: MPI
    INTEGER                   :: iElem,iNode,nDomElem 
    !--------------------------------------------------------------------------
    INTENT(IN)        :: EQN,IO,BND
    INTENT(INOUT)     :: MESH,DISC,MPI
    !--------------------------------------------------------------------------
    !                                                                         !
    logInfo0(*) '<--------------------------------------------------------->'        !
    logInfo0(*) '<                 MESH TREATMENT AND INFO                 >'        !
    logInfo0(*) '<--------------------------------------------------------->'        !
    logInfo0(*) ' '
    logInfo0(*) 'Reading mesh information from file ...'  !
    !
    NULLIFY(MPI%CPUDistribution) 
    ! 
    SELECT CASE (IO%meshgenerator)                                            !
    CASE('Gambit3D-Tetra')
       logInfo(*) 'Read a 3-D Tetrahedral Gambit mesh ...'   !
       !
       CALL read_mesh_Gambit3D_Tetra(                                       & !
            EQN         = EQN                                             , & !
            MESH        = MESH                                            , & ! 
            IO          = IO                                              , & ! 
            BND         = BND                                             , & !
            MPI         = MPI                                               ) ! 
       !                                                                      ! 
    CASE('Gambit3D-Hexa')
       logInfo(*) 'Read a 3-D Hexahedral Gambit mesh ...'    !
       !                                                                      !
       CALL read_mesh_Gambit3D_Hexa(                                        & !
            EQN         = EQN                                             , & !
            MESH        = MESH                                            , & ! 
            IO          = IO                                              , & ! 
            BND         = BND                                               ) ! 
       !    
    CASE('Gambit3D-Mixed')
       logInfo(*) 'Read a 3-D Mixed Gambit mesh ...'    !
       !                                                                      !
       CALL read_mesh_Gambit3D_Mixed(                                       & !
            EQN         = EQN                                             , & !
            MESH        = MESH                                            , & ! 
            IO          = IO                                              , & ! 
            BND         = BND                                             , & !
            DISC        = DISC                                            , & !
            MPI         = MPI                                               ) ! 
       !                                                                      ! 
    CASE('ICEMCFD3D-Tetra')
       logInfo(*) 'Read a 3-D ICEM CFD mesh ...'             !
       !                                                                      !
       CALL read_mesh_ICEMCFD3D_Tetra(                                      & !
            EQN         = EQN                                             , & !
            MESH        = MESH                                            , & ! 
            IO          = IO                                              , & ! 
            BND         = BND                                               ) ! 
       !                                                                      ! 

    CASE DEFAULT                                                              ! DEFAULT
       logError(*) 'No Routine to read a ',TRIM(IO%meshgenerator),' Mesh!'
       STOP                                                                   ! DEFAULT
    END SELECT                                                                !
    !                                                                         !
    CALL ObjectReferences_Side(MESH = MESH)                                   !
    !
    IF(.NOT.ASSOCIATED(MPI%CPUDistribution)) THEN
        CALL read_metis_data(nDomElem,MESH%nElem,EQN,MESH,IO,BND,MPI)
        ALLOCATE( MESH%NodeLocal2Global(MESH%nNode) ) 
        ALLOCATE( MESH%NodeGlobal2Local(MESH%nNode) ) 
        ALLOCATE( MESH%ElemLocal2Global(MESH%nElem) ) 
        ALLOCATE( MESH%ElemGlobal2Local(MESH%nElem) ) 
        DO iNode = 1, MESH%nNode
            MESH%NodeLocal2Global(iNode) = iNode 
            MESH%NodeGlobal2Local(iNode) = iNode 
        ENDDO
        DO iElem = 1, MESH%nElem 
            MESH%ElemLocal2Global(iElem) = iElem 
            MESH%ElemGlobal2Local(iElem) = iElem 
        ENDDO
    ENDIF
    !
    logInfo0(*) 'Mesh geometry extracted.'                                       !
    logInfo0(*) '<--------------------------------------------------------->'        !
    !
  END SUBROUTINE read_mesh

!------------------------------------------------------------------------
!> read_mesh_Gambit3D_Tetra reads a tetrahedral mesh produced by GAMBIT
!> in a memory efficient way. This way, each processor just keeps 
!> a) the elements of its subdomain and the corresponding vertices 
!> b) the elements connected to these vertices
!> c) all boundary elements to facilitate the creation of periodic boundaries.
!> The routine consists of three sweeps that loop through the mesh to extract 
!> only the necessary information.
!<
  SUBROUTINE read_mesh_Gambit3D_Tetra(EQN,MESH,IO,BND,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    TYPE (tEquations)         :: EQN
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tInputOutput)       :: IO
    TYPE (tBoundary)          :: BND
    TYPE (tMPI)               :: MPI 
    ! local (evtl. dummy) variables
    INTEGER                   :: ntotelem, ntotnode 
    INTEGER                   :: ielem, iZone                        !< index of each element
    INTEGER                   :: inode, iSide                        !< index of each node and side
    INTEGER                   :: openStat
    INTEGER                   :: idummy(1:10)
    INTEGER                   :: nSet, iSet, nBndEdges
    INTEGER                   :: GambitToSeisSol(4)
    INTEGER                   :: nLine, iLine, nRest, nDomElem
    INTEGER                   :: i, counter
    INTEGER                   :: itemp(4)
    INTEGER, POINTER          :: BndCond(:), TmpBndCond(:)
    INTEGER, POINTER          :: nZoneElements(:)
    LOGICAL, POINTER          :: ElementFlag(:)  
    INTEGER, POINTER          :: VertexFlag(:) 
    REAL                      :: xmin(3), xmax(3), xminglob(3), xmaxglob(3)
    CHARACTER(LEN=200)        :: cdummy, cdummy2     
    REAL                      :: temp(EQN%Dimension)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,IO,BND
    INTENT(INOUT)             :: MESH
    !--------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! FIRST SWEEP 
    ! reads vertices without storing
    ! flags elements of the subdomain and all boundary elements
    ! -------------------------------------------------------------------------
    xminglob =  1e13
    xmaxglob = -1e13

    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !

    IF (openStat .NE. 0) THEN                                                 ! Error control
       logError(*) 'cannot open ',TRIM(IO%MeshFile)        ! Error control
       STOP                                                                   ! Error control
    END IF                                                                    ! Error control

    logInfo(*) 'Assuming a tetrahedral mesh.'

    MESH%GlobalElemType = 4                                                   ! only tetrahedrons

    DO i = 1,6
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO

    READ(IO%UNIT%FileIn,*)   ntotnode,  ntotelem, MESH%nZones, nSet, idummy(1:2) 

    CALL read_metis_data( nDomElem,ntotelem,EQN,MESH,IO,BND,MPI )             ! read mesh partition file

    SELECT CASE(EQN%linType)                                                  ! check if zones coincide with material definition
    CASE(1)
      IF(ntotelem.EQ.EQN%nLayers)THEN
        logInfo(*)  'Assuming individual material for each element.'
      ELSE
        IF(MESH%nZones.NE.EQN%nLayers) THEN
          logError(*) 'Number of mesh zones and material definitions do not match. '
          STOP
        ENDIF
      ENDIF
    END SELECT

    IF(nSet.EQ.0) THEN                                                        ! Error control
      logError(*) 'Error in reading Gambit3D mesh. No boundary conditions specified!'
      STOP
    ELSE
      ALLOCATE(BndCond(nSet),TmpBndCond(nSet))
    ENDIF

    ALLOCATE(nZoneElements(MESH%nZones))

    DO i = 1,2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO

    logInfo(*)  'Reading node coordinates...'
    ! Read vertex coordinates
    DO iNode = 1, ntotnode
       READ(IO%UNIT%FileIn,*) idummy(1), temp(:)
       xminglob = MIN(xminglob,temp)
       xmaxglob = MAX(xmaxglob,temp)
    ENDDO

    DO i = 1,2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO

    logInfo(*)  'Reading connectivity...'
    ! Read connectivity (=elements) and flag vertices that belong to elements of processor MPI%myrank
    ALLOCATE( VertexFlag(ntotnode)  )
    ALLOCATE( ElementFlag(ntotelem) )
    VertexFlag  = 0
    ElementFlag = .FALSE.
    DO iElem = 1, ntotelem
       READ(IO%UNIT%FileIn,*) idummy(1:3), itemp(:)
       IF(MPI%CPUDistribution(iElem).EQ.MPI%myrank) THEN
         ElementFlag(iElem)   = .TRUE.
         VertexFlag(itemp(:)) = 1
       ENDIF
    ENDDO

    ! Read volume zones
    logInfo(*)  'Reading ', MESH%nZones, ' zones... '
    DO iZone = 1,MESH%nZones
        DO i = 1,2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        READ(IO%UNIT%FileIn,'(a28,i10,a40)') cdummy, nZoneElements(iZone), cdummy2 
        logInfo(*) 'Reading Mesh zone       : ', iZone
        logInfo(*) 'Number of zone elements : ', nZoneElements(iZone)
        DO i = 1, 2                                                               
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        nLine = nZoneElements(iZone)/10                                       ! Number of lines containing 10 entries
        nRest = MOD(nZoneElements(iZone),10)                                  ! Number of entries in the last line (rest)
        counter = 0                                                           ! Initialize counter
        ! Read zone elements for zone iZone
        DO iLine = 1, nLine                                                   ! Read full lines
          READ(IO%UNIT%FileIn,*) idummy(1:10)
        ENDDO
        IF(nRest.GT.0) THEN                                                   ! If there is a rest, read last line
          READ(IO%UNIT%FileIn,*) idummy(1:nRest)
        ENDIF
        logInfo(*) 'Zone read.'
    ENDDO

     ! Look for periodic boundary conditions and flag elements that are on boundary
    logInfo(*) 'Check for periodic boundary conditions'

    DO iSet = 1, nSet
        DO i = 1, 2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        READ(IO%UNIT%FileIn,*) BndCond(iSet), idummy(1), nBndEdges, idummy(2:3)
        IF(EQN%Adjoint.EQ.1.AND.BndCond(iSet).EQ.105) MESH%nOutBnd = nBndEdges
        ! Periodic boundaries 
        IF(MOD(BndCond(iSet),100).EQ.6) THEN
        	logInfo(*) 'Boundary condition',iSet,':',BndCond(iSet),'on',nBndEdges,'faces.'
            DO iSide = 1, nBndEdges
               READ(IO%UNIT%FileIn,*) idummy(1:3)
               ElementFlag(idummy(1)) = .TRUE. 
            ENDDO
            EXIT
        ELSE
            DO iSide = 1, nBndEdges
               READ(IO%UNIT%FileIn,*) idummy(1:3)
            ENDDO

        ENDIF
    ENDDO

    ! Deallocate temp. variables
    DEALLOCATE(nZoneElements)
    DEALLOCATE(BndCond, TmpBndCond)

    CLOSE (IO%UNIT%FileIn)                                                    ! closing input 

    ! -------------------------------------------------------------------------
    ! SECOND SWEEP 
    ! flags vertices without storing
    ! flags elements with previously flagged vertices
    ! -------------------------------------------------------------------------

    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !

    DO i = 1, 6
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO
    READ(IO%UNIT%FileIn,*)   ntotnode,  ntotelem, MESH%nZones, nSet, idummy(1:2) 

    DO i = 1, 2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO
    !                                                                         
    logInfo(*)  'Reading node coordinates...'                                !
    ! Read node coordinates and node references
    xmin =  1e13
    xmax = -1e13
    DO iNode = 1, ntotnode 
       READ(IO%UNIT%FileIn,*) idummy(1), temp(:)
       IF(VertexFlag(iNode).EQ.1) THEN
          xmin = MIN(xmin,temp) 
          xmax = MAX(xmax,temp) 
       ENDIF 
    ENDDO
    DO i = 1, 2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO
    logInfo(*)  'Reading connectivity...'
    ! Read connectivity and flag additional neighbor elements that share a flagged vertex
    DO iElem = 1, ntotelem
       READ(IO%UNIT%FileIn,*) idummy(1:3), itemp(:) 
       IF(.NOT.ElementFlag(iElem)) THEN 
           DO i = 1, 4
               IF(VertexFlag(itemp(i)).EQ.1) THEN 
                 ElementFlag(iElem) = .TRUE. 
                 EXIT 
               ENDIF 
           ENDDO
       ENDIF
       IF(ElementFlag(iElem)) THEN 
           DO i = 1, 4
               IF(VertexFlag(itemp(i)).EQ.0) THEN 
                 VertexFlag(itemp(i)) = 2  ! Flag with 2 to avoid growing element neighborhood
               ENDIF 
           ENDDO
       ENDIF
    ENDDO
    CLOSE(IO%UNIT%FileIn) 

    ! Count elements of subdomain, neighborhood and boudnary
    MESH%nElem = 0 
    DO iElem = 1, ntotelem
      IF(ElementFlag(iElem)) THEN
         MESH%nElem = MESH%nElem + 1 
      ENDIF
    ENDDO
    ! Count vertices of these elements
    MESH%nNode = 0 
    DO iElem = 1, ntotnode 
      IF(VertexFlag(iElem).GT.0) THEN
         MESH%nNode = MESH%nNode + 1 
      ENDIF
    ENDDO

    CALL allocate_mesh_level0_1(                                            & ! allocating array for elem
         IO   = IO                                                        , & ! allocating array for elem
         MESH = MESH                                                        ) ! allocating array for elem
    CALL allocate_mesh_level0_2(                                            & ! allocating array for node
         IO   = IO                                                        , & ! allocating array for node
         EQN  = EQN                                                       , & ! allocating array for node
         MESH = MESH                                                        ) ! allocating array for node
    !
    ALLOCATE(MESH%LocalElemType(MESH%nElem))                                  ! remains constant for single mesh type
    ALLOCATE(MESH%LocalVrtxType(MESH%nElem))
    
    MESH%LocalElemType(:) = 4
    MESH%LocalVrtxType(:) = 4
    MESH%nElem_Tet = MESH%nElem
    MESH%nElem_Hex = 0

    ! Initialize
    MESH%ELEM%SideNeighbor(:,:)      = MESH%nElem+1
    MESH%ELEM%LocalNeighborSide(:,:) = 0
    MESH%ELEM%LocalNeighborVrtx(:,:) = 0
    MESH%ELEM%Reference(:,:)         = 0

    ! -------------------------------------------------------------------------
    ! THIRD SWEEP 
    ! reads only flagged vertices
    ! reads only flagged elements
    ! reads volume zones
    ! reads boundaries
    ! -------------------------------------------------------------------------

    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !

    DO i = 1, 6
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO

    READ(IO%UNIT%FileIn,*)   ntotnode,  ntotelem, MESH%nZones, nSet, idummy(1:2) 

    ALLOCATE(BndCond(nSet),TmpBndCond(nSet))
    ALLOCATE(nZoneElements(MESH%nZones))

    DO i = 1, 2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO

    logInfo(*)  'Reading node coordinates...'
    ALLOCATE( MESH%NodeLocal2Global(MESH%nNode) ) 
    ALLOCATE( MESH%NodeGlobal2Local(ntotnode)   ) 
    ! Read node coordinates and store them
    counter = 0 
    DO iNode = 1, ntotnode 
       READ(IO%UNIT%FileIn,*) idummy(1), temp(:) 
       IF(VertexFlag(iNode).GT.0) THEN
         counter = counter + 1 
         MESH%VRTX%xyNode(:,counter) = temp(:) 
         MESH%NodeLocal2Global(counter) = iNode 
         MESH%NodeGlobal2Local(iNode)   = counter 
       ENDIF
    ENDDO
    DO i = 1, 2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO
    logInfo(*)  'Reading connectivity...'
    ! Read connectivity (=elements) and store them
    ALLOCATE( MESH%ElemLocal2Global(MESH%nElem) )
    ALLOCATE( MESH%ElemGlobal2Local(ntotelem)   )
    counter = 0
    DO iElem = 1, ntotelem
       READ(IO%UNIT%FileIn,*) idummy(1:3), itemp(:)
       IF(ElementFlag(iElem)) THEN
            counter = counter + 1
            MESH%ELEM%Vertex(:,counter) = MESH%NodeGlobal2Local( itemp(:) )
            MESH%ElemLocal2Global(counter)   = iElem
            MESH%ElemGlobal2Local(iElem)     = counter
       ENDIF
    ENDDO

    ! Read volume zones
    logInfo(*)  'Reading ', MESH%nZones, ' zones... '
    DO iZone = 1, MESH%nZones
        DO i = 1, 2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        READ(IO%UNIT%FileIn,'(a28,i10,a40)') cdummy, nZoneElements(iZone), cdummy2
        logInfo(*) 'Reading Mesh zone       : ', iZone
        logInfo(*) 'Number of zone elements : ', nZoneElements(iZone)
        DO i = 1, 2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        nLine = nZoneElements(iZone)/10                                       ! Number of lines containing 10 entries
        nRest = MOD(nZoneElements(iZone),10)                                  ! Number of entries in the last line (rest)
        counter = 0                                                           ! Initialize counter
        ! Read zone elements for zone iZone
        DO iLine = 1, nLine                                                   ! Read full lines
          READ(IO%UNIT%FileIn,*) idummy(1:10)
          DO i = 1, 10
            IF(ElementFlag(idummy(i))) THEN
                MESH%ELEM%Reference( 0,MESH%ElemGlobal2Local(idummy(i)) ) = iZone 
            ENDIF 
          ENDDO
        ENDDO
        IF(nRest.GT.0) THEN                                                   ! If there is a rest, read last line
          READ(IO%UNIT%FileIn,*) idummy(1:nRest)
          DO i = 1, nRest
            IF(ElementFlag(idummy(i))) THEN
                MESH%ELEM%Reference( 0,MESH%ElemGlobal2Local(idummy(i)) ) = iZone 
            ENDIF 
          ENDDO
        ENDIF        
        IF(EQN%nLayers.EQ.ntotelem)THEN                                     ! Individual material for each element
           DO iElem = 1, ntotelem
             IF(ElementFlag(iElem)) THEN
                MESH%ELEM%Reference(0,MESH%ElemGlobal2Local(iElem)) = iElem 
             ENDIF
           ENDDO
        ENDIF
        logInfo(*) 'Zone read. '
    ENDDO

    ! Read boundary conditions
    logInfo(*) 'No. of boundary conditions', nSet
    GambitToSeisSol(:) = (/1, 2, 4, 3/)   ! Map Gambit ordering convention for boundary conditions to SeisSol convention 
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
    !
    DO iSet = 1, nSet
        DO i = 1, 2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        READ(IO%UNIT%FileIn,*) BndCond(iSet), idummy(1), nBndEdges, idummy(2:3)
        logInfo(*) 'Boundary condition',iSet,':',BndCond(iSet),'on',nBndEdges,'faces.'
        IF(EQN%Adjoint.EQ.1.AND.BndCond(iSet).EQ.105) MESH%nOutBnd = nBndEdges
        DO iSide = 1, nBndEdges
           READ(IO%UNIT%FileIn,*) idummy(1:3)
           IF( ElementFlag(idummy(1)) ) THEN
                MESH%ELEM%Reference( GambitToSeisSol(idummy(3)),MESH%ElemGlobal2Local(idummy(1)) ) = BndCond(iSet)
           ENDIF
        ENDDO
    ENDDO
    !
    ! Deallocate temp. variables
    DEALLOCATE(nZoneElements)
    DEALLOCATE(BndCond, TmpBndCond)
    DEALLOCATE( VertexFlag   ) 
    DEALLOCATE( ElementFlag  ) 

    CLOSE (IO%UNIT%FileIn)                                                    ! closing input 

	! ----------------------------------------------------------------------------------------------
	! compute mesh and neightbour info
	CALL compute_mesh_Gambit3D_Tetra(Mesh, IO, EQN, BND)

    RETURN                                                                    !

  END SUBROUTINE read_mesh_Gambit3D_Tetra

  SUBROUTINE read_mesh_ICEMCFD3D_Tetra(EQN,MESH,IO,BND)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)         :: EQN
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tInputOutput)       :: IO
    TYPE (tBoundary)          :: BND
    ! local (evtl. dummy) variables
    INTEGER                   :: ielem, iZone                                  !< index of each element
    INTEGER                   :: inode, iSide							       !< index of each node and side
    INTEGER                   :: NodesTmp, GeomTmp, TypeTmp
    INTEGER                   :: i, j, k, testI, counter
    INTEGER                   :: openStat
    INTEGER                   :: idummy(1:10)
    INTEGER                   :: test_array_size
    INTEGER                   :: nSet, iSet, nBndCond, nBndElemTotal
    INTEGER, POINTER          :: BndCond(:), TmpBndCond(:), nBndElem(:)
    INTEGER, POINTER          :: nZoneElements(:)
    INTEGER, POINTER          :: ZoneElements(:)
    INTEGER, POINTER          :: BndFaceNodes(:,:), BndFaceNodesCond(:)
    INTEGER, POINTER          :: tmp1(:), tmp2(:,:), tmp3(:,:)
    INTEGER                   :: nLine, iLine, nRest
    !
    INTEGER                   :: allocStat

    CHARACTER(LEN=200)        :: cdummy
    CHARACTER(LEN=10)         :: cdummy1
    CHARACTER(LEN=12)         :: cdummy2         
    REAL                      :: temp(EQN%Dimension)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,IO,BND
    INTENT(INOUT)             :: MESH
    !--------------------------------------------------------------------------
    !                                                                         !
    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !
    !                                                                         !
    IF (openStat .NE. 0) THEN                                                 ! Fehlerkontrolle
       logError(*) 'cannot open ',TRIM(IO%MeshFile)
       STOP                                                                   ! Fehlerkontrolle
    END IF                                                                    ! Fehlerkontrolle
    !                                                                         !
    logInfo(*) 'Assuming a tetrahedral mesh.'                 !
    !                                                                   !
    MESH%GlobalElemType = 4                                                   ! Nur Tetras
    !                                                                         !
    DO i = 1, 5                                                               !
       READ(IO%UNIT%FileIn,*)                                                 ! Ignore header (read header away :-) )
    ENDDO
    !
    READ(IO%UNIT%FileIn,*)   MESH%nNode,  MESH%nIcemElem, MESH%nZonesTotal, idummy(1:2) 
    !
    ALLOCATE(nZoneElements(MESH%nZonesTotal))                                 ! Allocate maximum size
    !
    CALL allocate_mesh_level0_2(                                            & ! allocating array for node
         IO   = IO                                                        , & ! allocating array for node
         EQN  = EQN                                                       , & ! allocating array for node
         MESH = MESH                                                        ) ! allocating array for node
    
    ! Initialize    !
    DO i = 1, 7                                                               
       READ(IO%UNIT%FileIn,*)                                                 ! Ignore header (read header away :-) )
    ENDDO
    !                                                                         
    logInfo(*)  'Reading node coordinates...'                                !
    ! Read node coordinates and node references
    DO iNode = 1, MESH%nNode
       READ(IO%UNIT%FileIn,*) idummy(1), MESH%VRTX%xyNode(:,iNode)
    ENDDO
    !
    ! Do the node transformation
    !
    DO iNode = 1, MESH%nNode                                                                  !
      MESH%VRTX%xyNode(:,iNode) = MESH%VRTX%xyNode(:,iNode) + MESH%Displacement(:)            !
      Temp(:) = MATMUL( MESH%ScalingMatrix(:,:),MESH%VRTX%xyNode(:,iNode) )                   !
      MESH%VRTX%xyNode(:,iNode) = Temp(:)                                                     !
    ENDDO                                                                                     !
    !
    DO i = 1, 3                                                               
       READ(IO%UNIT%FileIn,*)                                                 ! Ignore header (read header away :-) )
    ENDDO
    !                                                                         !
    logInfo(*)  'Reading connectivity...'

    MESH%nElem = 0                                                            ! Initial mesh size
    ALLOCATE( MESH%ELEM%Vertex(MESH%GlobalVrtxType,MESH%nElem) )    
    ALLOCATE( MESH%ELEM%Reference(0:MESH%GlobalElemType,MESH%nElem) )
    MESH%ELEM%Reference = 0
    MESH%nZones = 0                                                           ! Initial number of zones         
    nSet = 0   
    !
    ! Read connectivity
    !
    DO iZone = 1, MESH%nZonesTotal
        ! Read zone header
        READ(IO%UNIT%FileIn,*) cdummy1, idummy(1), cdummy1, nZoneElements(iZone), cdummy1, NodesTmp, cdummy1, GeomTmp, cdummy1, &
TypeTmp
        READ(IO%UNIT%FileIn,*) cdummy2, cdummy2, cdummy1
        !
        SELECT CASE(NodesTmp)
        CASE(3) ! Boundary Zone of Triangular Elements

           nSet = nSet + 1
           !
           IF(nSet.EQ.1)THEN
               !
               nBndCond      = MESH%nZonesTotal - MESH%nZones          ! number of boundary conditions
               nBndElemTotal = MESH%nIcemElem - MESH%nElem             ! number of total boundary elements
               logInfo(*) 'No. of boundary conditions ',nBndCond
               ALLOCATE(BndFaceNodes(3,nBndElemTotal),BndFaceNodesCond(nBndElemTotal))
               ALLOCATE(nBndElem(nBndCond), BndCond(nBndCond))
               ! 
               counter = 0
               !
           ENDIF
           !
           nBndElem(nSet) = nZoneElements(iZone)      

           logInfo(*) 'Boundary condition ',nSet,' out of ',nBndCond
           !
           SELECT CASE(TRIM(cdummy1))
           CASE('BC101')
               BndCond(nSet) = 101 ! free surface
           CASE('BC103')
               BndCond(nSet) = 103 ! rupture surface
           CASE('BC104')
               BndCond(nSet) = 104 ! inflow
           CASE('BC105')
               BndCond(nSet) = 105 ! outflow
           CASE('BC106')
               BndCond(nSet) = 106 ! periodic
           CASE DEFAULT
               logError(*) 'ICEM-CFD-3D mesh uses non-defined boundary condition !'
               STOP
           END SELECT
         
           DO iElem = 1, nZoneElements(iZone)
               counter = counter + 1
               READ(IO%UNIT%FileIn,*) idummy(1), BndFaceNodes(:,counter)
               BndFaceNodesCond(counter) = BndCond(nSet)
           ENDDO

        CASE(4) ! Volume Zone of Tetrahedral Elements

           logInfo(*) 'Reading Mesh zone       : ', iZone
           ! enlarge arrays by using a temporary arrays
           ALLOCATE(tmp2(MESH%GlobalVrtxType,MESH%nElem ))
           ALLOCATE(tmp3(0:MESH%GlobalElemType,MESH%nElem ))
           DO i = 1, MESH%nElem
             DO j = 1, MESH%GlobalVrtxType
                tmp2(j,i) = MESH%ELEM%Vertex(j,i)
             ENDDO
             DO j = 0, MESH%GlobalElemType
                tmp3(j,i) = MESH%ELEM%Reference(j,i)
             ENDDO
           ENDDO
           DEALLOCATE(MESH%ELEM%Vertex)
           DEALLOCATE(MESH%ELEM%Reference)
           !
           MESH%nZones = MESH%nZones + 1
           MESH%nElem  = MESH%nElem  + nZoneElements(iZone)
           ! copy old array into larger new one
           ALLOCATE(MESH%ELEM%Vertex(MESH%GlobalVrtxType,MESH%nElem) ) 
           ALLOCATE(MESH%ELEM%Reference(0:MESH%GlobalElemType,MESH%nElem) )
           MESH%ELEM%Reference = 0
           DO i = 1, MESH%nElem - nZoneElements(iZone)
             DO j = 1, MESH%GlobalVrtxType
                MESH%ELEM%Vertex(j,i) = tmp2(j,i)
             ENDDO
             DO j = 0, MESH%GlobalElemType
                MESH%ELEM%Reference(j,i) = tmp3(j,i)
             ENDDO
           ENDDO
           DEALLOCATE(tmp2)
           DEALLOCATE(tmp3)
           !
           ! Read new elements
           DO iElem = MESH%nElem - nZoneElements(iZone) + 1, MESH%nElem
                READ(IO%UNIT%FileIn,*) idummy(1), MESH%ELEM%Vertex(:,iElem)
                MESH%ELEM%Reference(0,iElem) = iZone 
           ENDDO 
           !
           logInfo(*) 'Zone read. '
           !
        CASE DEFAULT  
             logError(*) 'ICEM-CFD-3D mesh has non-tetrahedral or non-triangular zone !'
             STOP
        END SELECT

    ENDDO  ! End nZonesTotal
    !
    CLOSE (IO%UNIT%FileIn)  ! closing input 

    ALLOCATE( MESH%ELEM%MPIreference(     0:MESH%nVertexMax,     MESH%nElem), &
              MESH%ELEM%BoundaryToObject( MESH%GlobalElemType,   MESH%nElem), &
              MESH%ELEM%SideNeighbor(     MESH%GlobalElemType,   MESH%nElem), & 
              MESH%ELEM%LocalNeighborSide(MESH%GlobalElemType,   MESH%nElem), & 
              MESH%ELEM%LocalNeighborVrtx(MESH%GlobalElemType,   MESH%nElem), & 
              MESH%ELEM%xyBary(           MESH%Dimension,        MESH%nElem), & 
              MESH%ELEM%Volume(           MESH%nElem                      ), &
              MESH%ELEM%MinDistBarySide(  MESH%nElem                      ), &
              STAT=allocStat                                              )

    IF (allocStat .NE. 0) THEN
        logError(*) 'could not allocate all mesh variables!'
        STOP
    END IF
    !

    MESH%ELEM%SideNeighbor(:,:)      = MESH%nElem+1
    MESH%ELEM%LocalNeighborSide(:,:) = 0
    MESH%ELEM%LocalNeighborVrtx(:,:) = 0

    ALLOCATE(MESH%LocalElemType(MESH%nElem))                              ! remains constant for single mesh type
    ALLOCATE(MESH%LocalVrtxType(MESH%nElem))
    
    MESH%LocalElemType(:) = 4
    MESH%LocalVrtxType(:) = 4
    MESH%nElem_Tet = MESH%nElem
    MESH%nElem_Hex = 0
    
    !
    ! Check zone numbers
    !
    SELECT CASE(EQN%linType)                                   
    CASE(1)
      IF(MESH%nElem.EQ.EQN%nLayers)THEN
         logInfo(*) 'Assuming individual material for each element.'
         DO iElem = 1,MESH%nElem
            MESH%ELEM%Reference(0,iElem) = iElem
         ENDDO
      ELSE                                                                                                            !
         IF(MESH%nZones.NE.EQN%nLayers) THEN
            logError(*) 'Number of mesh zones and material definitions do not match. '
            STOP
         ENDIF
      ENDIF
    END SELECT
    !
    ! Check bounary numbers
    !
    IF(nSet.EQ.0) THEN
      logError(*) 'Error in reading ICEM-CFD-3D mesh. No boundary conditions specified!'
      STOP
    ELSE
      ALLOCATE(TmpBndCond(nSet))
    ENDIF
    !
    !
    ! =======================================================
    ! Computing mesh information
    ! =======================================================
    !
    !
    logInfo(*) 'Building neighborhood information ... '

    CALL compute_mesh_ICEMCFD3D_Tetra(Mesh, IO, EQN, BND, BndCond, TmpBndCond, nBndElemTotal, BndFaceNodes,BndFaceNodesCond)
    
    DEALLOCATE(nZoneElements)
    RETURN                                                                    !
    !                                                                         !
  END SUBROUTINE read_mesh_ICEMCFD3D_Tetra                                    !


  SUBROUTINE read_mesh_Gambit3D_Hexa(EQN,MESH,IO,BND)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)         :: EQN
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tInputOutput)       :: IO
    TYPE (tBoundary)          :: BND
    ! local (evtl. dummy) variables
    INTEGER                   :: ielem, iZone                                   !< index of each element
    INTEGER                   :: inode, iSide					                !< index of each node and side
    INTEGER                   :: i, j, k, m, testI, counter
    INTEGER                   :: openStat
    INTEGER                   :: idummy(1:10)
    INTEGER                   :: test_array_size
    INTEGER                   :: nSet, iSet, nBndEdges
    INTEGER, POINTER          :: BndCond(:), TmpBndCond(:)
    INTEGER                   :: nBndObj(6)
    INTEGER, POINTER          :: nZoneElements(:)
    INTEGER, POINTER          :: ZoneElements(:)
    INTEGER                   :: nLine, iLine, nRest
    !
    CHARACTER(LEN=200)        :: cdummy, cdummy2     
    REAL                      :: temp(EQN%Dimension)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,IO,BND
    INTENT(INOUT)             :: MESH
    !--------------------------------------------------------------------------
    !                                                                         !
    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !
    !                                                                         !
    IF (openStat .NE. 0) THEN                                                 ! Fehlerkontrolle
       logError(*) 'cannot open ',TRIM(IO%MeshFile)
       STOP                                                                   ! Fehlerkontrolle
    END IF                                                                    ! Fehlerkontrolle
    !                                                                         !
    logInfo(*) 'Assuming a hexahedral mesh.'                !
    !                                                                         !
    MESH%GlobalElemType = 6                                                   ! Only Hexahedrons
    !                                                                         !
    DO i = 1, 6                                                               !
       READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
    ENDDO
    !
    READ(IO%UNIT%FileIn,*)   MESH%nnode,  MESH%nElem, MESH%nZones, nSet, idummy(1:2) 
    !       
    
    ALLOCATE(MESH%LocalElemType(MESH%nElem))                                  ! remains constant for single mesh type
    ALLOCATE(MESH%LocalVrtxType(MESH%nElem))
    
    MESH%LocalElemType(:) = 6
    MESH%LocalVrtxType(:) = 8
    MESH%nElem_Tet = 0
    MESH%nElem_Hex = MESH%nElem
    
    
    SELECT CASE(EQN%linType)                                   
    CASE(1)
      IF(MESH%nElem.EQ.EQN%nLayers)THEN
        logInfo(*) 'Assuming individual material for each element.'
      ELSE                                                                !                                                                        !
        IF(MESH%nZones.NE.EQN%nLayers) THEN
          logError(*) 'Number of mesh zones and material definitions do not match. '
          STOP
        ENDIF
      ENDIF
    END SELECT
    !
    IF(nSet.EQ.0) THEN
      logError(*) 'Error in reading Gambit3D mesh. No boundary conditions specified!'
      STOP
    ELSE
      ALLOCATE(BndCond(nSet),TmpBndCond(nSet))
    ENDIF
    ! 
    ALLOCATE(nZoneElements(MESH%nZones))
    !
    CALL allocate_mesh_level0_1(                                            & ! allocating array for elem
         IO   = IO                                                        , & ! allocating array for elem
         MESH = MESH                                                        ) ! allocating array for elem
    CALL allocate_mesh_level0_2(                                            & ! allocating array for node
         IO   = IO                                                        , & ! allocating array for node
         EQN  = EQN                                                       , & ! allocating array for node
         MESH = MESH                                                        ) ! allocating array for node
                                               !
    !                                                                         ! Initialize
    MESH%ELEM%SideNeighbor(:,:)      = MESH%nElem+1
    MESH%ELEM%LocalNeighborSide(:,:) = 0
    MESH%ELEM%LocalNeighborVrtx(:,:) = 0
    MESH%ELEM%Reference(:,:)         = 0
    !
    DO i = 1, 2                                                               
       READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
    ENDDO
    !                                                                         
    logInfo(*) 'Reading node coordinates...'                                !
    ! Read node coordinates and node references
    DO iNode = 1, MESH%nNode
       READ(IO%UNIT%FileIn,*) idummy(1), MESH%VRTX%xyNode(:,iNode)
    ENDDO
    !
    ! Do the node transformation
    !
    DO iNode = 1, MESH%nNode                                                                  !
      MESH%VRTX%xyNode(:,iNode) = MESH%VRTX%xyNode(:,iNode) + MESH%Displacement(:)            !
      Temp(:) = MATMUL( MESH%ScalingMatrix(:,:),MESH%VRTX%xyNode(:,iNode) )
      MESH%VRTX%xyNode(:,iNode) = Temp(:)
    ENDDO                                                                                     !
    !
    DO i = 1, 2                                                               
       READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
    ENDDO
    !                                                                         
    logInfo(*) 'Reading connectivity...'                                !
    ! Read connectivity
    DO iElem = 1, MESH%nElem
       READ(IO%UNIT%FileIn,*) idummy(1:3), MESH%ELEM%Vertex(1:7,iElem)
       READ(IO%UNIT%FileIn,*) MESH%ELEM%Vertex(8,iElem)
    ENDDO
    !
    ! Read volume zones
    !
    logInfo(*) 'Reading ', MESH%nZones, ' zones... '                                !
    DO iZone = 1, MESH%nZones        
        DO i = 1, 2                                                               
           READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
        ENDDO
        READ(IO%UNIT%FileIn,'(a28,i10,a40)') cdummy, nZoneElements(iZone), cdummy2 
        logInfo(*) 'Reading Mesh zone       : ', iZone
        logInfo(*) 'Number of zone elements : ', nZoneElements(iZone)
        DO i = 1, 2                                                               
           READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
        ENDDO
        ALLOCATE(ZoneElements(nZoneElements(iZone)))
        nLine = nZoneElements(iZone)/10                                           ! Number of lines containing 10 entries
        nRest = MOD(nZoneElements(iZone),10)                                      ! Number of entries in the last line (rest)
        counter = 0                                                               ! Initialize counter
        ! Read zone elements for zone iZone
        DO iLine = 1, nLine                                                       ! Read full lines
          READ(IO%UNIT%FileIn,*) idummy(1:10)
          DO i = 1, 10
           counter = counter + 1
           ZoneElements(counter) = idummy(i)
          ENDDO
        ENDDO
        IF(nRest.GT.0) THEN                                                       ! If there is a rest, read last line
          READ(IO%UNIT%FileIn,*) idummy(1:nRest)
          DO i = 1, nRest
             counter = counter + 1
             ZoneElements(counter) = idummy(i)
          ENDDO
        ENDIF
        DO iElem = 1, nZoneElements(iZone)
           MESH%ELEM%Reference(0,ZoneElements(iElem)) = iZone
           IF(EQN%nLayers.EQ.MESH%nElem)THEN                                      ! Individual material for each element
               MESH%ELEM%Reference(0,ZoneElements(iElem)) = ZoneElements(iElem)
           ENDIF
        ENDDO 
        DEALLOCATE(ZoneElements)
        logInfo(*) 'Zone read. '
    ENDDO
    ! 
    ! Read boundary conditions
    !
    logInfo(*) 'No. of boundary conditions', nSet


    DO iSet = 1, nSet
        DO i = 1, 2                                                               
           READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
        ENDDO
        READ(IO%UNIT%FileIn,*) BndCond(iSet), idummy(1), nBndEdges, idummy(2:3)
        logInfo(*) 'Boundary condition',iSet,':',BndCond(iSet),'on',nBndEdges,'faces.'
        DO iSide = 1, nBndEdges
           READ(IO%UNIT%FileIn,*) idummy(1:3)
           MESH%ELEM%Reference(idummy(3),idummy(1)) = BndCond(iSet)
        ENDDO
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
    DEALLOCATE(nZoneElements)
    DEALLOCATE(BndCond, TmpBndCond)
    !
    CLOSE (IO%UNIT%FileIn)                                                    ! closing input 
    !
    ! =======================================================
    ! Computing mesh information
    ! =======================================================
    !
    !
    logInfo(*) 'Building neighborhood information ... '
    !
    CALL compute_mesh_Gambit3D_Hexa(Mesh, IO, EQN, BND)
    
    RETURN                                                                    !
    !                                                                         !
  END SUBROUTINE read_mesh_Gambit3D_Hexa                                           !

SUBROUTINE read_mesh_Gambit3D_Mixed(EQN,MESH,IO,BND,DISC,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif            
    TYPE (tEquations)         :: EQN
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tInputOutput)       :: IO
    TYPE (tBoundary)          :: BND
    TYPE (tDiscretization)    :: DISC  
    TYPE (tMPI)               :: MPI   
    ! local (evtl. dummy) variables
    INTEGER                   :: ielem, iZone                        !< index of each element
    INTEGER                   :: inode, iSide   			         !< index of each node and side
    INTEGER                   :: i, j, k, testI, counter
    INTEGER                   :: openStat
    INTEGER                   :: idummy(1:10)
    INTEGER                   :: test_array_size
    INTEGER                   :: nSet, iSet, nBndEdges
    INTEGER, POINTER          :: BndCond(:), TmpBndCond(:)
    INTEGER                   :: nBndObj(6)
    INTEGER, POINTER          :: nZoneElements(:)
    INTEGER, POINTER          :: ZoneElements(:)
    INTEGER                   :: GambitToSeisSol(4)
    INTEGER                   :: nLine, iLine, nRest
    INTEGER                   :: nDomElem,tot_tetra,tot_hexa,ntotelem,ntotnode
    INTEGER                   :: itemp(8)
    INTEGER, POINTER          :: VertexFlag(:),ElementType(:),ElementTypeIndex(:)
    LOGICAL, POINTER          :: ElementFlag(:)
    !
    CHARACTER(LEN=200)        :: cdummy, cdummy2     
    REAL                      :: temp(EQN%Dimension)
    !
    ! Dynamic Rupture variables
    INTEGER                   :: TmpFaultElements
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,IO,BND
    INTENT(INOUT)             :: MESH, DISC
    !--------------------------------------------------------------------------
    !     

           
    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !
    !                                                                         !
    IF (openStat .NE. 0) THEN                                                 ! Fehlerkontrolle
       logError(*) 'cannot open ',TRIM(IO%MeshFile)
       STOP                                                                   ! Fehlerkontrolle
    END IF                                                                    ! Fehlerkontrolle
    !                                                                         !
    logInfo(*)  '    Assuming a mixed or non-conforming mesh.'   !
    !                                                                   !
    MESH%GlobalElemType = 7                                                   ! Mixed
    !                                                                         !
    DO i = 1, 6                                                               !
       READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
    ENDDO
    !
    READ(IO%UNIT%FileIn,*)   MESH%nNode,  MESH%nElem, MESH%nZones, nSet, idummy(1:2) 
#ifdef HDF
    MESH%nNode_total = MESH%nNode 
    MESH%nElem_total = MESH%nElem
#endif
    !       

    ! Returns number of elements for each processor
    CALL read_metis_data(nDomElem,MESH%nElem,EQN,MESH,IO,BND,MPI)           
        
    SELECT CASE(EQN%linType)                                   
    CASE(1)
      IF(MESH%nElem.EQ.EQN%nLayers)THEN
        logInfo(*)  'Assuming individual material for each element.'
      ELSE                                                                !                                                                        !
        IF(MESH%nZones.NE.EQN%nLayers) THEN
          logError(*) 'Number of mesh zones and material definitions do not match. '
          STOP
        ENDIF
      ENDIF
    END SELECT
    !
    IF(nSet.EQ.0) THEN
      logError(*) 'Error in reading Gambit3D mesh. No boundary conditions specified!'
      STOP
    ELSE
      ALLOCATE(BndCond(nSet),TmpBndCond(nSet))
    ENDIF
    
    ALLOCATE(nZoneElements(MESH%nZones))    
       
    DO i = 1, 2                                                               
       READ(IO%UNIT%FileIn,*)                                          
    ENDDO
                                                                             
    logInfo(*) 'Skipping node coordinates...'
    ! Read node coordinates and node references
    DO iNode = 1, MESH%nNode
       READ(IO%UNIT%FileIn,*)
    ENDDO
                                                                                              
    DO i = 1, 2                                                               
       READ(IO%UNIT%FileIn,*)                                             
    ENDDO
    

    tot_tetra = 0
    tot_hexa = 0
    
    ! CHECK WHETHER ELEMENT IS TETRA OR HEXA AND COMPUTE TOTAL NUMBER OF TET & HEX
    logInfo(*) 'Reading element type...'
    DO iElem = 1, MESH%nElem
       READ(IO%UNIT%FileIn,*) idummy(1:3)
       IF (idummy(3) .EQ. 4) THEN
!          ElementType(iElem) = 4
          tot_tetra = tot_tetra + 1
       ELSEIF (idummy(3) .EQ. 8) THEN
!          ElementType(iElem) = 6
          tot_hexa = tot_hexa + 1
       ELSE
          logError(*) 'Mixed mesh does not contain only tetrahedrons and hexahedrons!'
          STOP
       ENDIF           
    ENDDO   
          
    CLOSE(IO%UNIT%FileIn)       
    
    logInfo(*)  'Total Number of Tets inside entire mesh: ', tot_tetra
    logInfo(*)  'Total Number of Hexas inside entire mesh: ', tot_hexa
       
    
   !-----------------------------------------------
   ! FIRST SWEEP -- FLAG ELEMENTS AND VERTEX CORRESPONDING TO EACH PROCESSOR
    
    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !
                                                                            
    !                                                                         !
    DO i = 1, MESH%nNode + 11                                                      
       READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
    ENDDO
    
    logInfo(*) 'Reading connectivity...'
    
    ! Read connectivity (=elements) and flag vertices that belong to elements of processor MPI%myrank
    ALLOCATE( VertexFlag(MESH%nNode)  )
    ALLOCATE( ElementFlag(MESH%nElem) )
    
    VertexFlag  = 0
    ElementFlag = .FALSE.
    
    ! FLAGGING THOSE ELEMENTS BELONGING TO THE LOCAL SUB-DOMAIN
    DO iElem = 1, MESH%nElem
       READ(IO%UNIT%FileIn,*) idummy(1:3), itemp(1:idummy(3))
       IF(MPI%CPUDistribution(iElem).EQ.MPI%myrank) THEN
          ElementFlag(iElem)   = .TRUE.
          VertexFlag(itemp(1:idummy(3))) = 1       
       ENDIF 
    ENDDO   
      
    ! COUNT ELEMENTS BELONGING TO THE LOCAL SUB-DOMAIN     
    counter=0    
    DO iElem = 1, MESH%nElem 
       IF(ElementFlag(iElem)) THEN
          counter = counter + 1     
       ENDIF
    ENDDO     

    IF (counter .NE. nDomElem) THEN
       logError(*) 'Number of Elements in the Sub-Region Not Corresponding to Read_metis_data'
       STOP
    ENDIF


    ! SKIP MESH ZONES IN THE FIRST SWEEP (i.e. MATERIALS)
    logInfo(*) 'Skipping ', MESH%nZones, ' zones... '                                !
    DO iZone = 1, MESH%nZones        
        DO i = 1, 2                                                               
           READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
        ENDDO
        READ(IO%UNIT%FileIn,'(a28,i10,a40)') cdummy, nZoneElements(iZone), cdummy2 
        logInfo(*) 'Reading Mesh zone       : ', iZone
        logInfo(*) 'Number of zone elements : ', nZoneElements(iZone)
        DO i = 1, 2                                                               
           READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
        ENDDO
        nLine = nZoneElements(iZone)/10                                           ! Number of lines containing 10 entries
        nRest = MOD(nZoneElements(iZone),10)                                      ! Number of entries in the last line (rest)
        counter = 0                                                               ! Initialize counter
        ! Read zone elements for zone iZone
        DO iLine = 1, nLine                                                       ! Read full lines
          READ(IO%UNIT%FileIn,'(10i8)') idummy(1:10)
        ENDDO
        IF(nRest.GT.0) THEN                                                       ! If there is a rest, read last line
          READ(IO%UNIT%FileIn,'(10i8)') idummy(1:nRest)
        ENDIF
        logInfo(*) 'Zone skipped. '
    ENDDO
    
    
    ! FLAGS TETS AND HEX AT BOUNDARIES IF PERIODIC,DYNAMIC OR NON-CONFORMING BC ARE MET     
    logInfo(*) 'Check for periodic, non-conforming and dynamic rupture boundary conditions'

    DO iSet = 1, nSet
        DO i = 1, 2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        READ(IO%UNIT%FileIn,*) BndCond(iSet), idummy(1), nBndEdges, idummy(2:3)
        IF(EQN%Adjoint.EQ.1.AND.BndCond(iSet).EQ.105) MESH%nOutBnd = nBndEdges
        ! Periodic, dynamic rupture or non-conforming boundaries 
        IF( (MOD(BndCond(iSet),100).EQ.6) .OR. (MOD(BndCond(iSet),100).EQ.3)             &
                                          .OR. (MOD(BndCond(iSet),100).EQ.2) ) THEN
        	logInfo(*) 'Boundary condition',iSet,':',BndCond(iSet),'on',nBndEdges,'faces.'
            DO iSide = 1, nBndEdges
               READ(IO%UNIT%FileIn,*) idummy(1:3)
               ElementFlag(idummy(1)) = .TRUE. 
            ENDDO
            EXIT
        ELSE
            DO iSide = 1, nBndEdges
               READ(IO%UNIT%FileIn,*) idummy(1:3)
            ENDDO

        ENDIF
    ENDDO

    ! Deallocate temp. variables
    DEALLOCATE(nZoneElements)
    DEALLOCATE(BndCond, TmpBndCond)

    CLOSE (IO%UNIT%FileIn)   
    ! END FIRST SWEEP -----------------------------
    
    
    !-----------------------------------------------
    ! SECOND SWEEP -- LOOKING FOR NEIGHBORING ELEMENTS AND RELATIVE VERTEXES
    
    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !
                                                                            
    ! Skip nodes coordinates                                                  !
    DO i = 1, MESH%nNode + 11                                                      
       READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
    ENDDO
    
    logInfo(*) 'Reading connectivity...'
    
    ! READ CONNECTIVITY AND FLAG ADDITIONAL NEIGHBOR ELEMENTS THAT SHARE A FLAGGED VERTEX
    ! WITH THOSE ELEMENTS ALREADY FLAGGED (SUB-DOMAIN + EVENTUALLY BC)
    DO iElem = 1, MESH%nElem

       READ(IO%UNIT%FileIn,*) idummy(1:3), itemp(1:idummy(3)) 
       IF(.NOT.ElementFlag(iElem)) THEN 
           DO i = 1, idummy(3)
               IF(VertexFlag(itemp(i)).EQ.1) THEN 
                 ElementFlag(iElem) = .TRUE. 
                 EXIT 
               ENDIF 
           ENDDO
       ENDIF
       IF(ElementFlag(iElem)) THEN 
           DO i = 1, idummy(3)
               IF(VertexFlag(itemp(i)).EQ.0) THEN 
                 VertexFlag(itemp(i)) = 2  ! Flag with 2 to avoid growing element neighborhood
               ENDIF 
           ENDDO
       ENDIF
    ENDDO
    CLOSE(IO%UNIT%FileIn) 
    
    
    ! COUNT ELEMENTS OF SUBDOMAIN + NEIGHBORHOOD + BOUNDARY
    ntotelem = MESH%nElem
    MESH%nElem = 0 
    
    ! THIS VARIABLE IS NEEDED TO DISTINGUISH BETWEEN TETS & HEXAS 
    ALLOCATE(ElementTypeIndex(ntotelem))
    
    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !
                                                                            
    ! Skip nodes coordinates                                                  !
    DO i = 1, MESH%nNode + 11                                                      
       READ(IO%UNIT%FileIn,*)                                                 ! Lese Header weg
    ENDDO
    
    logInfo(*) 'Computing number of Tets & Hexa per processor...'
    DO iElem = 1, ntotelem
       READ(IO%UNIT%FileIn,*) idummy(1:3), itemp(1:idummy(3))    
       ! COMPUTE THE TOTAL NUMBER OF ELEMENTS PER PROCESSOR
       IF(ElementFlag(iElem)) THEN
         MESH%nElem = MESH%nElem + 1  
         ElementTypeIndex(MESH%nElem) = idummy(3)     
      ENDIF
    ENDDO     
    CLOSE(IO%UNIT%FileIn) 

    
    ! COUNT NUMBER OF ELEMENT VERTEXES PER PROCESSOR
    ntotnode = MESH%nNode
    MESH%nNode = 0 
    DO iElem = 1, ntotnode 
      IF(VertexFlag(iElem).GT.0) THEN
         MESH%nNode = MESH%nNode + 1 
      ENDIF
    ENDDO
    
    CALL allocate_mesh_level0_2(                                            & ! allocating array for node
         IO   = IO                                                        , & ! allocating array for node
         EQN  = EQN                                                       , & ! allocating array for node
         MESH = MESH                                                        ) ! allocating array for node
        
              
    ALLOCATE(MESH%LocalElemType(MESH%nElem))
    ALLOCATE(MESH%LocalVrtxType(MESH%nElem))
    MESH%nElem_Tet = 0
    MESH%nElem_Hex = 0  
                          
    
    ! COMPUTE NUMBER OF TET & HEX ELEMENTS PER PROCESSOR
    DO iElem = 1, MESH%nElem
    
       MESH%LocalElemType(iElem) = ElementTypeIndex(iElem)
       
       IF (MESH%LocalElemType(iElem) .EQ. 4) MESH%LocalVrtxType(iElem) = 4
       IF (MESH%LocalElemType(iElem) .EQ. 6) MESH%LocalVrtxType(iElem) = 8
         
       IF ( (MESH%LocalElemType(iElem) .NE. 4) .AND. (MESH%LocalElemType(iElem) .NE. 6) ) THEN
            logError(*) 'Mixed mesh does not contain only tetrahedrons and hexahedrons!'
            STOP 
       ENDIF

       SELECT CASE(MESH%LocalElemType(iElem))
        CASE(4)
          MESH%nElem_Tet = MESH%nElem_Tet + 1
          IF(MESH%nElem_Tet.EQ.1)THEN
            ALLOCATE(MESH%LocalElemIndex_Tet(MESH%nElem))
            MESH%LocalElemIndex_Tet = 0
          ENDIF
          MESH%LocalElemIndex_Tet(iElem) = MESH%nElem_Tet 
        CASE(6)
          MESH%nElem_Hex = MESH%nElem_Hex + 1
          IF(MESH%nElem_Hex.EQ.1)THEN
            ALLOCATE(MESH%LocalElemIndex_Hex(MESH%nElem))
            MESH%LocalElemIndex_Hex = 0
          ENDIF
          MESH%LocalElemIndex_Hex(iElem) = MESH%nElem_Hex 
          READ(IO%UNIT%FileIn,*)
        CASE DEFAULT
          logError(*) 'Mixed mesh does not contain only tetrahedrons and hexahedrons!'
          STOP 
      END SELECT
    ENDDO      

    DEALLOCATE(ElementTypeIndex)
        
    logInfo(*) 'Number of Tets per processor: ', MESH%nElem_Tet
    logInfo(*) 'Number of Hexas per processor: ', MESH%nElem_Hex
    logInfo(*) 'Number of Elements per processor: ', MESH%nElem


    ! Check mesh type again (really mixed or non-conforming)    
    IF(MESH%nElem_Tet.GT.0.AND.MESH%nElem_Hex.EQ.0)THEN     
        MESH%GlobalElemType = 7
        MESH%nVertexMax     = 4
    ENDIF                                                                  
    IF(MESH%nElem_Tet.EQ.0.AND.MESH%nElem_Hex.GT.0)THEN     
        MESH%GlobalElemType = 7
        MESH%nVertexMax     = 8
    ENDIF                                                                  
    IF(MESH%nElem_Tet.GT.0.AND.MESH%nElem_Hex.GT.0)THEN     
        MESH%GlobalElemType = 7
        MESH%nVertexMax     = 8
    ENDIF    
       

    CALL allocate_mesh_level0_1(                                            & ! allocating array for elem
         IO   = IO                                                        , & ! allocating array for elem
         MESH = MESH                                                        ) ! allocating array for elem  
       

    ! Initialize
    MESH%ELEM%SideNeighbor(:,:)      = MESH%nElem+1
    MESH%ELEM%LocalNeighborSide(:,:) = 0
    MESH%ELEM%LocalNeighborVrtx(:,:) = 0
    MESH%ELEM%Reference(:,:)         = 0
    MESH%ELEM%MPIReference(:,:)      = 0
    MESH%ELEM%BoundaryToObject(:,:)  = 0
    
    MESH%ELEM%Vertex(:,:)            = 0.0d0
    MESH%ELEM%xyBary(:,:)            = 0.0d0
    MESH%ELEM%Volume(:)              = 0.0d0
    MESH%ELEM%MinDistBarySide(:)     = 0.0d0
    
    ! END SECOND SWEEP -----------------------------


    ! THIRD SWEEP -------------------------------
    ! READ ONLY FLAGGED ELEMENTS
    
    OPEN(UNIT   = IO%UNIT%FileIn                                          , & ! opening unit 
         FILE   = TRIM(IO%MeshFile)                                       , & ! input mesh file
         STATUS = 'OLD'                                                   , & !
         ACTION = 'READ'                                                  , & !
         ACCESS = 'SEQUENTIAL'                                            , & !
         IOSTAT = openStat                                                  ) !
    !                                                                         !
    
    DO i = 1, 6
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO

    READ(IO%UNIT%FileIn,*)   ntotnode,  ntotelem, MESH%nZones, nSet, idummy(1:2) 

    ALLOCATE(BndCond(nSet),TmpBndCond(nSet))
    ALLOCATE(nZoneElements(MESH%nZones))

    DO i = 1, 2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO

    logInfo(*) 'Reading node coordinates...'

    ALLOCATE(MESH%NodeLocal2Global(MESH%nNode))
    ALLOCATE(MESH%NodeGlobal2Local(ntotnode))

    counter = 0

    DO iNode = 1,ntotnode 
       READ(IO%UNIT%FileIn,*) idummy(1),temp(:)
       IF (VertexFlag(iNode) .GT. 0) THEN
          counter = counter + 1
          MESH%VRTX%xyNode(:,counter) = temp(:)
          
          MESH%NodeLocal2Global(counter) = iNode
          MESH%NodeGlobal2Local(iNode) = counter             
       ENDIF
    ENDDO
    
    ! Do the node transformation
    DO iNode = 1, MESH%nNode
      MESH%VRTX%xyNode(:,iNode) = MESH%VRTX%xyNode(:,iNode) + MESH%Displacement(:)            
      Temp(:) = MATMUL( MESH%ScalingMatrix(:,:),MESH%VRTX%xyNode(:,iNode) )
      MESH%VRTX%xyNode(:,iNode) = Temp(:)
    ENDDO  

    DO i = 1, 2
       READ(IO%UNIT%FileIn,*)                                                 ! read header
    ENDDO
    
    logInfo(*) 'Reading connectivity...'
    
    ! Map global element list to local (MPIrank) element list
    ALLOCATE(MESH%ElemLocal2Global(MESH%nElem))
    ALLOCATE(MESH%ElemGlobal2Local(ntotelem))
    
    counter = 0
        
    DO iElem = 1,ntotelem 
       READ(IO%UNIT%FileIn,*) idummy(1:3),itemp(1:idummy(3))
       IF (ElementFlag(iElem)) THEN
          counter = counter + 1
          MESH%ELEM%Vertex(:,counter) = MESH%NodeGlobal2Local(itemp(1:MESH%LocalVrtxType(counter)))
          
          MESH%ElemLocal2Global(counter) = iElem
          MESH%ElemGlobal2Local(iElem) = counter
       ENDIF
    ENDDO      
    
    
    ! READ VOLUME ZONE (i.e. MATERIALS) -- ONLY FLAGGED ELEMENTS  
    ! Read volume zones
    logInfo(*) 'Reading ', MESH%nZones, ' zones... '
    DO iZone = 1, MESH%nZones
        DO i = 1, 2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        READ(IO%UNIT%FileIn,'(a28,i10,a40)') cdummy, nZoneElements(iZone), cdummy2
        logInfo(*) 'Reading Mesh zone       : ', iZone
        logInfo(*) 'Number of zone elements : ', nZoneElements(iZone)
        DO i = 1, 2
           READ(IO%UNIT%FileIn,*)                                             ! read header
        ENDDO
        nLine = nZoneElements(iZone)/10                                       ! Number of lines containing 10 entries
        nRest = MOD(nZoneElements(iZone),10)                                  ! Number of entries in the last line (rest)
        counter = 0                                                           ! Initialize counter
        ! Read zone elements for zone iZone
        DO iLine = 1, nLine                                                   ! Read full lines
          READ(IO%UNIT%FileIn,'(10i8)') idummy(1:10)
          DO i = 1, 10
            IF(ElementFlag(idummy(i))) THEN
                MESH%ELEM%Reference( 0,MESH%ElemGlobal2Local(idummy(i)) ) = iZone 
            ENDIF 
          ENDDO
        ENDDO
        IF(nRest.GT.0) THEN                                                   ! If there is a rest, read last line
          READ(IO%UNIT%FileIn,'(10i8)') idummy(1:nRest)
          DO i = 1, nRest
            IF(ElementFlag(idummy(i))) THEN
                MESH%ELEM%Reference( 0,MESH%ElemGlobal2Local(idummy(i)) ) = iZone 
            ENDIF 
          ENDDO
        ENDIF        
        IF(EQN%nLayers.EQ.ntotelem)THEN                                     ! Individual material for each element
           DO iElem = 1, ntotelem
             IF(ElementFlag(iElem)) THEN
                MESH%ELEM%Reference(0,MESH%ElemGlobal2Local(iElem)) = iElem 
             ENDIF
           ENDDO
        ENDIF
        logInfo(*) '    Zone read. '
    ENDDO


    ! READ BOUNDARY CONDITIONS -- ONLY FLAGGED ELEMENTS (SPECIFIC BOUNDARY CONDITIONS)
    !
    logInfo(*) 'No. of boundary conditions', nSet
    GambitToSeisSol(:) = (/1, 2, 4, 3/)   ! Map Gambit ordering convention for boundary conditions to SeisSol convention 
    
    MESH%nNonConformingEdges = 0
    MESH%ELEM%Reference(1:MESH%nVertexMax,:) = 0
    
    DO iSet = 1, nSet
        DO i = 1, 2                                                               
           READ(IO%UNIT%FileIn,*)                                             
        ENDDO
        READ(IO%UNIT%FileIn,*) BndCond(iSet), idummy(1), nBndEdges, idummy(2:3)
        logInfo(*) 'Boundary condition',iSet,':',BndCond(iSet),'on',nBndEdges,'faces.'
        IF(BndCond(iSet).EQ.102) MESH%nNonConformingEdges = nBndEdges
        IF(BndCond(iSet).EQ.103) TmpFaultElements = nBndEdges
        IF(EQN%Adjoint.EQ.1.AND.BndCond(iSet).EQ.105) MESH%nOutBnd = nBndEdges
        
        DO iSide = 1, nBndEdges
           READ(IO%UNIT%FileIn,*) idummy(1:3)
  
           IF(idummy(2) .EQ. 6)THEN ! gambit ID for TETS
              IF( ElementFlag(idummy(1)) ) THEN
                 MESH%ELEM%Reference( GambitToSeisSol(idummy(3)),MESH%ElemGlobal2Local(idummy(1)) ) = BndCond(iSet)
                !MESH%ELEM%Reference(GambitToSeisSol(idummy(3)),idummy(1)) = BndCond(iSet)
              ENDIF 
              
           ! FOR HEXA AN UPDATE IS NEEDED (requires Hexa GambitToSeisSol list) 
           ELSEIF(idummy(2) .EQ. 4)THEN ! gambit ID for hexas
              MESH%ELEM%Reference(idummy(3),idummy(1)) = BndCond(iSet) 
           ELSE
              logError(*) 'Error in Gambit3D mesh. Element cannot be identified as tet or hex.'
              stop
           ENDIF
        ENDDO
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
    DEALLOCATE(nZoneElements)
    DEALLOCATE(BndCond, TmpBndCond)
    !
    CLOSE (IO%UNIT%FileIn)                                                    ! closing input 
    !
    ! =======================================================
    ! Computing mesh information
    ! =======================================================
    !
    !
    logInfo(*) 'Building neighborhood information ... '
    !
    CALL compute_mesh_Gambit3D_Mixed(Mesh, IO, EQN, BND, DISC, TmpFaultElements)

    RETURN                                                                    !
    !                                                                         !
  END SUBROUTINE read_mesh_Gambit3D_Mixed 



  SUBROUTINE ObjectReferences_Side(MESH)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH
    !--------------------------------------------------------------------------
    INTEGER                   :: i,j
    !--------------------------------------------------------------------------
    !                                                                         !
    DO j=1,MESH%nElem                                                         !
       DO i=1,MESH%nVertexMax                                                 !
          MESH%ELEM%BoundaryToObject(i,j) = INT(MESH%ELEM%Reference(i,j)/100) ! Mapping to Object Number
          MESH%ELEM%Reference(i,j)     = MOD(MESH%ELEM%Reference(i,j),100) ! Generate basic reference number from 0 to 99
       ENDDO                                                                  !
    ENDDO                                                                     !
    !                                                                         !
  END SUBROUTINE ObjectReferences_Side                                        !



!------------------------------------------------------------------------
!> read_metis_data reads the content of the Metis mesh partition file
!> and counts the number of elements for each processor MPI%myrank
!<
  SUBROUTINE read_metis_data(nDomElem,ntotelem,EQN,MESH,IO,BND,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)         :: EQN
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tBoundary)          :: BND
    TYPE (tMPI)               :: MPI
    TYPE (tInputOutput)       :: IO
    INTEGER                   :: nDomElem
    !--------------------------------------------------------------------------    
    INTEGER                    :: iElem,ntotelem,stat
    CHARACTER(LEN=600)         :: epartFileName
    CHARACTER(LEN=10)          :: cCPU
    !--------------------------------------------------------------------------
    INTENT(IN)        :: EQN,MESH,IO,BND
    INTENT(INOUT)     :: MPI
    INTENT(OUT)       :: nDomElem 
    !--------------------------------------------------------------------------
    !
    IF(MPI%nCPU.GT.1) THEN
            IF(MPI%nCPU.LT.10) THEN
              WRITE(cCPU,'(i1)') MPI%nCPU
            ELSEIF( (MPI%nCPU.GE.10  ) .AND. (MPI%nCPU.LT.100  ) )       THEN
              WRITE(cCPU,'(i2)') MPI%nCPU
            ELSEIF( (MPI%nCPU.GE.100 ) .AND. (MPI%nCPU.LT.1000 ) )       THEN
              WRITE(cCPU,'(i3)') MPI%nCPU
            ELSEIF( (MPI%nCPU.GE.1000) .AND. (MPI%nCPU.LT.10000) )       THEN
              WRITE(cCPU,'(i4)') MPI%nCPU
            ELSEIF( (MPI%nCPU.GE.10000) .AND. (MPI%nCPU.LT.100000) )     THEN
              WRITE(cCPU,'(i5)') MPI%nCPU
            ELSEIF( (MPI%nCPU.GE.100000) .AND. (MPI%nCPU.LT.1000000) )   THEN
              WRITE(cCPU,'(i6)') MPI%nCPU
            ELSEIF( (MPI%nCPU.GE.1000000) .AND. (MPI%nCPU.LT.10000000) ) THEN
              WRITE(cCPU,'(i7)') MPI%nCPU
            ELSE
              PRINT *, 'Character overflow in METIS Filename!'
              STOP
            ENDIF

            epartFileName=TRIM(IO%MetisFile)//'.epart.'//TRIM(cCPU) 

            logInfo(*) 'Reading METIS information from file ', TRIM(epartFileName)

            OPEN(UNIT   = 99                                                      , & ! opening 
                 FILE   = TRIM(epartFileName)                                     , & ! epart file
                 STATUS = 'OLD'                                                   , & ! generated
                 ACTION = 'READ'                                                  , & ! by METIS.
                 ACCESS = 'SEQUENTIAL'                                            , & !
                 IOSTAT = stat                                                    ) !

            IF(stat.NE.0) THEN
                logError(*) 'METIS File ', TRIM(epartFileName), ' not found. '
                STOP
            ENDIF
            
            ALLOCATE(MPI%CPUDistribution(ntotelem))
            ! Read the distribution of the CPUs given by METIS
            nDomElem = 0 
            DO iElem = 1, ntotelem
              READ(99,*) MPI%CPUDistribution(iElem)
              IF(MPI%CPUDistribution(iElem).EQ.MPI%myrank) THEN
                nDomElem = nDomElem + 1 
              ENDIF
            ENDDO
            CLOSE( 99 ) 
    ELSE
        ALLOCATE(MPI%CPUDistribution(ntotelem))
        MPI%CPUDistribution(:) = MPI%myrank        
        nDomElem = ntotelem 
    ENDIF
    ! 
  END SUBROUTINE read_metis_data

END MODULE read_mesh_mod
