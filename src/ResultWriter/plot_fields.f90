!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2010-2016, SeisSol Group
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

MODULE plot_fields_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE 
  !----------------------------------------------------------------------------
  INTERFACE ini_plot_fields
     MODULE PROCEDURE ini_plot_fields
  END INTERFACE

  INTERFACE close_plot_fields
     MODULE PROCEDURE close_plot_fields
  END INTERFACE

  INTERFACE plot_fields
     MODULE PROCEDURE plot_fields
  END INTERFACE

  INTERFACE write_tecplot_node_data
     MODULE PROCEDURE write_tecplot_node_data
  END INTERFACE

  INTERFACE Verlinken
     MODULE PROCEDURE Verlinken1D, Verlinken2D
  END INTERFACE

  INTERFACE InterpolateBaryToNode
     MODULE PROCEDURE InterpolateBaryToNode
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC  :: ini_plot_fields
  PUBLIC  :: close_plot_fields 
  PUBLIC  :: plot_fields
  PUBLIC  :: InterpolateBaryToNode
  !----------------------------------------------------------------------------
  INTEGER                   :: tbody

CONTAINS


  SUBROUTINE ini_plot_fields(pvar, OptionalFields, EQN, DISC, MESH, IO)
    !--------------------------------------------------------------------------
    IMPLICIT NONE     
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tEquations)              :: EQN
    TYPE (tDiscretization), target :: DISC
    TYPE (tUnstructMesh)           :: MESH
    TYPE (tUnstructOptionalFields) :: OptionalFields
    TYPE (tInputOutput)            :: IO
    REAL                           :: pvar(MESH%nElem,EQN%nVar)
    ! local variable declaration
    INTEGER                        :: allocStat,i,n,EndLink
    !--------------------------------------------------------------------------
    INTENT(IN)                     :: pvar, EQN, MESH ,IO
    INTENT(INOUT)                  :: OptionalFields
    !--------------------------------------------------------------------------
    !                                                     !
    ALLOCATE(OptionalFields%FieldMask(IO%nOutputMask) , & ! Allocate weight
         OptionalFields%weight(MESH%nElem)        , & ! Allocate weight
         STAT = allocStat                           ) ! Allocate weight
    !                                                   !
    IF (allocstat .NE. 0 ) THEN                         !
       logError(*)                        & ! Error Handler
            'ALLOCATE ERROR in plot_fields (weight)!'  ! Error Handler
       STOP                                             ! Error Handler
    END IF                                              ! Error Handler
    IF(EQN%CartesianCoordSystem) THEN                   ! Calc weight
       OptionalFields%weight(:) = 1.0                   ! Calc weight
    ELSE                                                ! Calc weight
       OptionalFields%weight(:) = MESH%ELEM%xyBary(2,:) ! Calc weight
    ENDIF                                               ! Calc weight
    !                                                   !
       

          CALL Verlinken(OptionalFields%FieldMask( 1)%PTR,IO%OutputMask( 1),Mesh%ELEM%xyBary,1)
          CALL Verlinken(OptionalFields%FieldMask( 2)%PTR,IO%OutputMask( 2),Mesh%ELEM%xyBary,2)
          CALL Verlinken(OptionalFields%FieldMask( 3)%PTR,IO%OutputMask( 3),Mesh%ELEM%xyBary,3)
          ! @breuera :)
          optionalFields%fieldMask( 4)%ptr => disc%galerkin%dgvar(1,1,:,1)
          optionalFields%fieldMask( 5)%ptr => disc%galerkin%dgvar(1,2,:,1)
          optionalFields%fieldMask( 6)%ptr => disc%galerkin%dgvar(1,3,:,1)
          optionalFields%fieldMask( 7)%ptr => disc%galerkin%dgvar(1,4,:,1)
          optionalFields%fieldMask( 8)%ptr => disc%galerkin%dgvar(1,5,:,1)
          optionalFields%fieldMask( 9)%ptr => disc%galerkin%dgvar(1,6,:,1)
          optionalFields%fieldMask(10)%ptr => disc%galerkin%dgvar(1,7,:,1)
          optionalFields%fieldMask(11)%ptr => disc%galerkin%dgvar(1,8,:,1)
          optionalFields%fieldMask(12)%ptr => disc%galerkin%dgvar(1,9,:,1)
          !
          EndLink = 12
          IF(EQN%Poroelasticity.NE.0) THEN
            CALL Verlinken(OptionalFields%FieldMask(13)%PTR,IO%OutputMask(13),pvar       ,10)
            CALL Verlinken(OptionalFields%FieldMask(14)%PTR,IO%OutputMask(14),pvar       ,11)
            CALL Verlinken(OptionalFields%FieldMask(15)%PTR,IO%OutputMask(15),pvar       ,12)
            CALL Verlinken(OptionalFields%FieldMask(16)%PTR,IO%OutputMask(16),pvar       ,13)
            EndLink = 16
          ENDIF
  
          IF(EQN%Plasticity.NE.0) THEN
             optionalFields%fieldMask( 13)%ptr => disc%galerkin%pstrain(1,:)
             optionalFields%fieldMask( 14)%ptr => disc%galerkin%pstrain(2,:)
             optionalFields%fieldMask( 15)%ptr => disc%galerkin%pstrain(3,:)
             optionalFields%fieldMask( 16)%ptr => disc%galerkin%pstrain(4,:)
             optionalFields%fieldMask( 17)%ptr => disc%galerkin%pstrain(5,:)
             optionalFields%fieldMask( 18)%ptr => disc%galerkin%pstrain(6,:)
             optionalFields%fieldMask( 19)%ptr => disc%galerkin%pstrain(7,:)

            EndLink = 19
          ENDIF
          !
          DO i = 1, EQN%nBackgroundVar-EQN%nMechanisms*4  ! Link constants 
              CALL Verlinken(OptionalFields%FieldMask(EndLink+i)%PTR,IO%OutputMask(EndLink+i),OptionalFields%BackgroundValue,i)
          ENDDO
          !
      IF(DISC%Galerkin%pAdaptivity.GT.0) THEN
         CALL Verlinken(OptionalFields%FieldMask(59)%PTR,IO%OutputMask(59),DISC%Galerkin%LocPoly )
      ENDIF
      IF(DISC%Galerkin%DGMethod.EQ.3) THEN
         CALL Verlinken(OptionalFields%FieldMask(60)%PTR,IO%OutputMask(60),DISC%LocalTime )
      ENDIF

  END SUBROUTINE ini_plot_fields                        !

  SUBROUTINE close_plot_fields(OptionalFields,IO)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE     
    !--------------------------------------------------------------------------
    TYPE (tInputOutput)            :: IO
    TYPE (tUnstructOptionalFields) :: OptionalFields
    INTEGER :: I
    !--------------------------------------------------------------------------
    INTENT(INOUT)                  :: OptionalFields,IO
    !--------------------------------------------------------------------------
    !                                                   !
    DO I = 1,IO%nOutputMask                             !
       NULLIFY(OptionalFields%FieldMask(I)%PTR)         !
    END DO                                              !
    !                                                   !
    IF(ASSOCIATED(OptionalFields%weight)) THEN          !
       DEALLOCATE(OptionalFields%weight)                !
    ENDIF                                               !
! aheineck, @TODO, not referecned in the code, commented!
!    IF(ASSOCIATED(OptionalFields%Mach)) THEN            !
!       DEALLOCATE(OptionalFields%Mach)                  !
!    ENDIF                                               !
!    IF(ASSOCIATED(OptionalFields%Entropie)) THEN        !
!       DEALLOCATE(OptionalFields%Entropie)              !
!    ENDIF                                               !
!    IF(ASSOCIATED(OptionalFields%Temp)) THEN            !
!       DEALLOCATE(OptionalFields%Temp)                  !
!    ENDIF                                               !
!    IF(ASSOCIATED(OptionalFields%Alf)) THEN             !
!       DEALLOCATE(OptionalFields%Alf)                   !
!    ENDIF                                               !
!    IF(ASSOCIATED(OptionalFields%rot)) THEN             !
!       DEALLOCATE(OptionalFields%rot)                   !
!    ENDIF                                               !
!    IF(ASSOCIATED(OptionalFields%div)) THEN             !
!       DEALLOCATE(OptionalFields%div)                   !
!    ENDIF                                               !
!    IF(ASSOCIATED(OptionalFields%rotB)) THEN            !
!       DEALLOCATE(OptionalFields%rotB)                  !
!    ENDIF                                               !
!    IF(ASSOCIATED(OptionalFields%divB)) THEN            !
!       DEALLOCATE(OptionalFields%divB)                  !
!    ENDIF  
!    IF(ASSOCIATED(OptionalFields%grd)) THEN            !
!       DEALLOCATE(OptionalFields%grd)                  !
!    ENDIF                                              !
    !                                                   !
  END SUBROUTINE close_plot_fields                      !

  SUBROUTINE plot_fields(fileop, time, timestep, EQN, DISC, MESH, IO, BND, OptionalFields, MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE     
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tEquations)         :: EQN
    TYPE (tDiscretization)    :: DISC
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tInputOutput)       :: IO
    TYPE (tBoundary)          :: BND
    TYPE (tUnstructOptionalFields) :: OptionalFields
    TYPE (tMPI)               :: MPI
    CHARACTER (LEN=*)         :: fileop     ! output file
    INTEGER                   :: timestep
    REAL                      :: time
    ! local variables declaration
    INTEGER                   :: stat, stat_tet, stat_hex
    CHARACTER (LEN=60)        :: purfile    ! filename without extension
    CHARACTER(LEN=9)          :: citer      ! characters for timestep
    LOGICAL                   :: binoutput  ! Data-Output in binary Format
!    CHARACTER (LEN=350)       :: Filename
    INTEGER                   :: i,j,n,iVar,nValues,allocStat
    REAL                      :: nodeVal
    REAL,POINTER              :: nodes(:,:)
    INTEGER,POINTER           :: cons(:,:)
    REAL,POINTER              :: scalars(:,:), vectors(:,:)
    REAL, PARAMETER           :: rho_ref=1., p_ref=1.
    !--------------------------------------------------------------------------
    INTENT(IN)                :: fileop, time, timestep
    INTENT(INOUT)             :: OptionalFields, MESH, IO
    !--------------------------------------------------------------------------
    !
    IF(IO%Format.eq.1)THEN
            OPEN(UNIT   = IO%UNIT%FileOut_Tet                   , &  ! open file
              FILE   = TRIM(fileop)//'.tet.'//TRIM(IO%Extension), &  ! open file
              STATUS = 'UNKNOWN'                                , &  ! open file
              RECL   = 5000                                     , &  ! open file
              IOSTAT = stat_tet                                   )  ! open file
    END IF
    !                                                   !
    !                                                   !
    !                                                   !
    ! --------------------------------------------------!
    ! Plot cases                                        !
    ! --------------------------------------------------!
    ! IO%Format=1: node data in TECPLOT-format
    ! IO%Format=5: XDMF (requires HDF5)
    ! --------------------------------------------------!
    !                                                   !
    !                                                   !
    SELECT CASE(IO%Format)                              !
    CASE(1)
        !                                             ! Tecplot format
        CALL write_tecplot_node_data(OptionalFields, EQN, DISC, &
            MESH, IO, timestep, time)
         !                                              !
    case(5)
       ! Data in XDMF format (requires HDF5)
       logError(*) 'Output format 5 (legacy HDF5) is no longer supported'
       stop
         !                                              !
    CASE DEFAULT                                        ! Error Handler
       logError(*)                        & ! Error Handler
            'Error in plot_fields.f90!'                 ! Error Handler
       logError(*)                        & ! Error Handler
            ' IO%Format, CASE DEFAULT'                  ! Error Handler
       STOP                                             ! Error Handler
    END SELECT                                          !
    !                                                   !       
    if (IO%Format.eq.1) then
            CLOSE (IO%UNIT%FileOut_Tet)                 ! close data file
    end if
    !                                                   !
    RETURN                                              !
                                                        !
  END SUBROUTINE plot_fields                            !



  SUBROUTINE write_tecplot_node_data(OptionalFields, EQN, DISC, &
       MESH, IO, timestep, time)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE     
    !--------------------------------------------------------------------------
    TYPE (tUnstructOptionalFields) :: OptionalFields
    TYPE (tEquations)         :: EQN
    TYPE (tDiscretization)    :: DISC
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tInputOutput)       :: IO
    INTEGER                   :: timestep
    REAL                      :: time
    ! local variables declaration
    INTEGER                   :: I
    INTEGER                   :: ielem
    INTEGER                   :: iNode
    INTEGER                   :: iOutVar
    INTEGER                   :: OutVarStart
    INTEGER                   :: eType
    REAL                      :: nodeVal
    !--------------------------------------------------------------------------
    INTENT(IN)                :: OptionalFields,MESH,IO,time
    !--------------------------------------------------------------------------
    ! title
    ! case of unsteady calculation (0), header with current time
    ! case of   steady calculation (1), header with current it.number     
    !                                                       !
    OutVarStart = 4

    WRITE(IO%UNIT%FileOut_Tet,*)'TITLE = "CURRENT TIME ', time,' "' ! unsteady case
    WRITE(IO%UNIT%FileOut_Tet,*) IO%Title               !
    WRITE(IO%UNIT%FileOut_Tet,*)                      & ! Tetrahedrons
        'ZONE N=',MESH%nnode ,                        & ! Tetrahedrons
        '  E='   ,MESH%nElem_Tet ,                    & ! Tetrahedrons
        '  F=FEPOINT  ET=TETRAHEDRON '                  ! Tetrahedrons

    DO inode=1,MESH%nnode
        WRITE(IO%UNIT%FileOut_Tet,104,advance='no')   & ! x-Value
            (MESH%VRTX%xyNode(1,inode))                 !
        !                                               !
        WRITE(IO%UNIT%FileOut_Tet,104,advance='no')   & ! y-Value
            (MESH%VRTX%xyNode(2,inode))                 !
        !                                               !
        WRITE(IO%UNIT%FileOut_Tet,104,advance='no')   & ! z-Value
            (MESH%VRTX%xyNode(3,inode))                 !

        DO iOutVar=OutVarStart,IO%nOutputMask               ! all other Variables
            IF(IO%OutputMask(iOutVar)) THEN                 !
                !                                           !
                CALL InterpolateBaryToNode(              &  !
                    BaryField = OptionalFields%FieldMask(iOutVar)%PTR, & !
                    weight    = OptionalFields%weight , &
                    NodeValue = NodeVal               , &   !
                    iNode     = iNode                 , &   !
                    MESH      = MESH                    )   !
                !                                           !
                WRITE(IO%UNIT%FileOut_Tet,104,advance='no') NodeVal ! Node Value of VarNr
                !                                           ! iOutVar
            ENDIF                                           !
        ENDDO ! iOutVar                                     !
       !                                                    !
       WRITE(IO%UNIT%FileOut_Tet,*) ! Next line
    ENDDO
       
    ! Write element - node connectivity
    DO iElem=1,MESH%nElem                               
        eType = MESH%LocalElemType(iElem)
        WRITE(IO%UNIT%FileOut_Tet,*) MESH%ELEM%Vertex((/1,2,3,4/),iElem)                   !
    ENDDO ! iElem
    !                                                   !
104 FORMAT (e25.14)                                     !
    !                                                   ! 
  END SUBROUTINE write_tecplot_node_data                !

  SUBROUTINE Verlinken1D(Zeiger, Maske, Ziel)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    REAL,  POINTER     :: Zeiger(:)
    LOGICAL            :: Maske
    REAL, TARGET       :: Ziel(:)
    !--------------------------------------------------------------------------
    INTENT(IN)         :: Maske, Ziel
    !--------------------------------------------------------------------------
    !                                                   !
    IF(Maske) THEN                                      !
       Zeiger => Ziel(:)                                !
    ELSE                                                !
       NULLIFY(Zeiger)                                  !
    ENDIF                                               !
    !                                                   !
  END SUBROUTINE Verlinken1D

  SUBROUTINE Verlinken2D(Zeiger, Maske, Ziel,n)
    !--------------------------------------------------------------------------
    IMPLICIT NONE     
    !--------------------------------------------------------------------------
    REAL,  POINTER     :: Zeiger(:)
    LOGICAL            :: Maske
    REAL, TARGET       :: Ziel(:,:)
    INTEGER            :: n
    !--------------------------------------------------------------------------
    INTENT(IN)         :: Maske, Ziel
    !--------------------------------------------------------------------------
    !                                                   !
    IF(Maske) THEN                                      !
       Zeiger => Ziel(:,n)                              !
    ELSE                                                !
       NULLIFY(Zeiger)                                  !
    ENDIF                                               !
    !                                                   !
  END SUBROUTINE Verlinken2D                            !

  SUBROUTINE InterpolateBaryToNode(BaryField,weight,NodeValue,iNode,MESH)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE     
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tUnstructMesh)   :: MESH
    REAL                   :: BaryField(MESH%nElem)
    REAL                   :: weight(MESH%nElem)
    REAL                   :: NodeValue
    INTEGER                :: iNode
    ! local variable declaration
    INTEGER                :: iElem
    INTEGER                :: ElemNr
    REAL                   :: vol
    REAL                   :: work
    !--------------------------------------------------------------------------
    INTENT(IN)  :: BaryField, weight, iNode, MESH
    INTENT(OUT) :: NodeValue
    !--------------------------------------------------------------------------
    ! interpolating cell average field into vertex field!
    !                                                   !
    vol  = 0.0                                          !
    work = 0.0                                          !
    !                                                   !
    DO ielem = 1,MESH%VRTX%NrOfElementsConnected(inode) !
       !
       ElemNr = MESH%VRTX%Element(inode,ielem)          !
       !
       vol  = vol                                     & !
            + MESH%ELEM%Volume(ElemNr)                & !
            * weight(ElemNr)                            !
       !                                                !
       work = work                                    & !
            + BaryField(  ElemNr)                     & !
            * MESH%ELEM%Volume(ElemNr)                & !
            * weight(     ElemNr)                       ! cartesian: 
       !                                                ! weight = 1.0 = const.
       !                                                ! cylindrical:
       !                                                ! weight = radius
    ENDDO                                               !
    !                                                   !
    NodeValue = work/vol                                !
    !
    IF(ABS(NodeValue).LT.1e-20) THEN
       NodeValue = 0.
    ENDIF
    !                                                   !
  END SUBROUTINE InterpolateBaryToNode                  !

END MODULE plot_fields_mod
