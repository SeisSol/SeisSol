!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Atanas Atanasov (atanasoa AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Atanas_Atanasov)
!! @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2012-2017, SeisSol Group
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
!! Routines initializing fault output for dynamic rupture simulations

#include "Initializer/preProcessorMacros.fpp"

MODULE ini_faultoutput_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE COMMON_operators_mod
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  PRIVATE :: fill_variables_list
  PRIVATE :: eval_faultreceiver
  PRIVATE :: create_file_headers
  PRIVATE :: write_header_info_to_files
  PRIVATE :: construct_file_name

  PUBLIC  :: ini_fault_subsampled
  PUBLIC  :: ini_fault_receiver
  PUBLIC  :: ini_fault_xdmfwriter
  !---------------------------------------------------------------------------!
  INTERFACE ini_fault_receiver
     MODULE PROCEDURE ini_fault_receiver
  END INTERFACE
  INTERFACE eval_faultreceiver
     MODULE PROCEDURE eval_faultreceiver
  END INTERFACE
  INTERFACE ini_fault_subsampled
     MODULE PROCEDURE ini_fault_subsampled
  END INTERFACE

  INTERFACE ini_fault_xdmfwriter
    MODULE PROCEDURE ini_fault_xdmfwriter
  END INTERFACE
CONTAINS


  !---------------------------------------------------------------------------!
  ! Subroutine checks whick variables are to output and fills the list with
  ! their names
  SUBROUTINE fill_variables_list(DISC, VariableList)

    IMPLICIT NONE

    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    CHARACTER (len=30000)   :: VariableList, VariableList_temp
    CHARACTER(len=3)        :: VName(8)
    INTEGER                 :: j

    ! Full list of possible variable names
    VName = (/'SRs','SRd','T_s','T_d','P_n','u_n','Mud','StV'/)

    ! Prepare second header line
    VariableList = TRIM('VARIABLES = "Time"')
    !
    DO j=1,8
       VariableList_temp=TRIM(VariableList)
       IF (j.EQ.1.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(1).NE.1) CYCLE
       IF (j.EQ.2.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(1).NE.1) CYCLE
       IF (j.EQ.3.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(2).NE.1) CYCLE
       IF (j.EQ.4.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(2).NE.1) CYCLE
       IF (j.EQ.5.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(2).NE.1) CYCLE
       IF (j.EQ.6.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(3).NE.1) CYCLE
       IF (j.EQ.7.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(4).NE.1) CYCLE
       IF (j.EQ.8.AND.DISC%DynRup%DynRup_out_atPickpoint%OutputMask(4).NE.1) CYCLE

       WRITE(VariableList,'(a,a3,a,a1)')                   &
       TRIM(VariableList_temp),',"',TRIM(VName(j)),'"'
    ENDDO
  END SUBROUTINE



  !---------------------------------------------------------------------------!
  SUBROUTINE construct_file_name(DISC, IO, MPI, receiver_index, ptsoutfile)
    IMPLICIT NONE

    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    TYPE(tInputOutput)      :: IO                     ! IO structure          !
    TYPE(tMPI)              :: MPI                    ! MPI                   !
    INTEGER                 :: receiver_index
    CHARACTER (len=200)     :: ptsoutfile

    CHARACTER (LEN=5)       :: cmyrank

#ifdef PARALLEL
    WRITE(cmyrank,'(I5.5)') MPI%myrank                                   ! myrank -> cmyrank
    WRITE(ptsoutfile, '(a,a15,i5.5,a1,a5,a4)') TRIM(IO%OutputFile),'-faultreceiver-',DISC%DynRup%DynRup_out_atPickpoint%RecPoint(receiver_index)%globalreceiverindex,'-',TRIM(cmyrank),'.dat'!
#else
    WRITE(ptsoutfile, '(a,a15,i5.5,a4)') TRIM(IO%OutputFile),'-faultreceiver-',DISC%DynRup%DynRup_out_atPickpoint%RecPoint(receiver_index)%globalreceiverindex,'.dat'
#endif
  END SUBROUTINE


  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  SUBROUTINE create_file_headers(DISC, IO, MPI, VariableList)

    IMPLICIT NONE

    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    TYPE(tInputOutput)      :: IO                     ! IO structure          !
    TYPE(tMPI)              :: MPI                    ! MPI                   !
    CHARACTER (len=30000)   :: VariableList

    !local variables
    INTEGER                 :: i
    LOGICAL                 :: exist
    INTEGER                 :: stat
    CHARACTER (LEN=5)       :: cmyrank
    CHARACTER (len=200)     :: ptsoutfile

    ALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%VFile(DISC%DynRup%DynRup_out_atPickpoint%nDR_pick))
    !for saving picked temporal values of the variablesiSide
    DO i = 1, DISC%DynRup%DynRup_out_atPickpoint%nDR_pick
      !Giving the VFiles unit file numbers (hopefully unit numbers over 25000 are not occupied)
      DISC%DynRup%DynRup_out_atPickpoint%VFile(i) = 25000 + (i-1)
      !
      CALL construct_file_name(DISC, IO, MPI, i, ptsoutfile)

      logInfo(*) '... open file ', TRIM(ptsoutfile),' to save time signal'
      logInfo(*) '<--------------------------------------------------------->'
      logInfo(*) ' '
      !
      !
      INQUIRE(FILE = ptsoutfile, EXIST = exist)
      IF(exist) THEN
        ! If file exists, then append data
        OPEN(UNIT     = DISC%DynRup%DynRup_out_atPickpoint%VFILE(i)      , & !
             FILE     = ptsoutfile                                       , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'OLD'                                            , & !
             POSITION = 'APPEND'                                         , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        IF(stat.NE.0) THEN                                                   !
           logError(*) 'cannot open ',ptsoutfile          !
           logError(*) 'Error status: ', stat
           STOP                                                              !
        END IF                                                               !
        CLOSE( DISC%DynRup%DynRup_out_atPickpoint%VFILE(i) )
      ELSE
        ! If file does not exist, then write header
        OPEN(UNIT     = DISC%DynRup%DynRup_out_atPickpoint%VFile(i)      , & !
             FILE     = ptsoutfile                                       , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'NEW'                                            , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        !                                                                    !
        IF(stat.NE.0) THEN                                                   !
           logError(*) 'cannot open ',ptsoutfile           !
           logError(*) 'Error status: ', stat
           STOP                                                              !
        END IF                                                               !
        !
        ! Creating the header for the output File
        ! First line of the header
        !
        WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a,I8.8,a)') &
            'TITLE = "Temporal Signal for fault receiver number ', DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%globalreceiverindex, ' "'
        !
        ! Second line of the header (variable names)
        WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a)') TRIM(VariableList)
        !
        ! Third-fifth line of the header (comment, containing x positions)
        !
        WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a4,e25.12)') TRIM('# x1'), DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%X
        WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a4,e25.12)') TRIM('# x2'), DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%Y
        WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a4,e25.12)') TRIM('# x3'), DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%Z
        !
        CLOSE( DISC%DynRup%DynRup_out_atPickpoint%VFile(i) )
        !
      ENDIF
    ENDDO ! i = 1, DISC%DynRup%DynRup_out_atPickpoint%nOutPoints
  END SUBROUTINE


  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  SUBROUTINE write_header_info_to_files(EQN, MESH, DISC, IO, MPI)

    USE create_fault_rotationmatrix_mod

    IMPLICIT NONE

    TYPE(tEquations)        :: EQN                    ! Equation structure    !
    TYPE(tUnstructMesh)     :: MESH                   ! Mesh structure        !
    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    TYPE(tInputOutput)      :: IO                     ! IO structure          !
    TYPE(tMPI)              :: MPI                    ! MPI                   !

    CHARACTER (LEN=5)       :: cmyrank
    CHARACTER (len=200)     :: ptsoutfile
    INTEGER                 :: i
    INTEGER                 :: iBndGP, iFault
    INTEGER                 :: stat
    REAL                    :: tmp_mat(6)       ! temporary stress tensors
    REAL                    :: rotmat(1:6,1:6)  ! rotation matrix for individual fault receiver

    DO i = 1, DISC%DynRup%DynRup_out_atPickpoint%nDR_pick
      ! output distance (in reference units to StdOut)
      logInfo(*) 'Fault pickpoint ',i,' of this MPI domain has its next GP at', &
                 DISC%DynRup%DynRup_out_atPickpoint%OutInt_dist(i),'.'

      ! add background stress values to fault pickpoint header

      CALL construct_file_name(DISC, IO, MPI, i, ptsoutfile)

      OPEN(UNIT     = DISC%DynRup%DynRup_out_atPickpoint%VFile(i)      , & !
           FILE     = ptsoutfile                                       , & !
           FORM     = 'FORMATTED'                                      , & !
           STATUS   = 'OLD'                                            , & !
           POSITION = 'APPEND'                                         , & !
           RECL     = 80000                                            , & !
           IOSTAT = stat                                                 ) !
      IF( stat.NE.0) THEN
        logError(*) 'cannot open ',ptsoutfile
        logError(*) 'Error status: ', stat
        STOP
      END IF
      !
      ! load stress values of nearest boundary GP: iBndGP
      iBndGP = DISC%DynRup%DynRup_out_atPickpoint%OutInt(i,1)
      iFault = DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%index
      !
      ! Create rotation matrix that rotates solution to strike and dip directions on arbitrary shaped faults
      CALL create_fault_rotationmatrix(rotmat,iFault,EQN,MESH)
      DISC%DynRup%DynRup_out_atPickpoint%rotmat(i,1:6,1:6) = rotmat
      !
      ! create tmp_mat = full tensor for rotation
      tmp_mat = 0.0D0
      tmp_mat(1) = EQN%IniBulk_xx(iBndGP,iFault)
      tmp_mat(2) = EQN%IniBulk_yy(iBndGP,iFault)
      tmp_mat(3) = EQN%IniBulk_zz(iBndGP,iFault)
      tmp_mat(4) = EQN%IniShearXY(iBndGP,iFault)
      tmp_mat(5) = EQN%IniShearYZ(iBndGP,iFault)
      tmp_mat(6) = EQN%IniShearXZ(iBndGP,iFault)
      !
      ! rotate into fault system
      tmp_mat(:) = MATMUL(rotmat(1:6,1:6),tmp_mat(:))
      !
      ! Sixth-eighth line of the header (comment, containing background stress)
      WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a5,e25.12)') TRIM('# P_0'), tmp_mat(1)
      WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a5,e25.12)') TRIM('# T_s'), tmp_mat(4)
      WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(i),'(a5,e25.12)') TRIM('# T_d'), tmp_mat(6)
      !
      CLOSE( DISC%DynRup%DynRup_out_atPickpoint%VFile(i) )
      !
    ENDDO ! Determine basis functions for DR pickpoint output and store them in OutEval

  END SUBROUTINE

  !=====================================================================!
  !!>                                                                  !!
  !!> ini_faultreceiver tests if the fault receiver are located        !!
  !!> on the fault within a certain tolerance, determines              !!
  !!> the basis functions, opens the regarding output files            !!
  !!> and allocates the necessary variables                            !!
  !!>                                                                  !!
  !!>                     for DYNAMIC RUPTURE                          !!
  !!>                 by C. Pelties, 13.10.10                          !!
  !!===================================================================!!

  SUBROUTINE ini_fault_receiver(EQN,MESH,BND,DISC,IO,MPI)
    !--------------------------------------------------------------------------!
    USE create_fault_rotationmatrix_mod
    USE DGBasis_mod
    USE common_fault_receiver_mod
    !--------------------------------------------------------------------------!
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)        :: EQN                    ! Equation structure    !
    TYPE(tUnstructMesh)     :: MESH                   ! Mesh structure        !
    TYPE(tBoundary)         :: BND
    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    TYPE(tInputOutput)      :: IO                     ! IO structure          !
    TYPE(tMPI)              :: MPI                    ! MPI                   !
    ! Local variable declaration
    INTEGER                 :: i,j,l,number_of_inner_receivers,iPick                   ! loop index
    INTEGER                 :: iElem, iSide, iNeighbor
    INTEGER                 :: stat

    REAL                    :: xi,eta,zeta,xi_n,eta_n,zeta_n,xi_GP,eta_GP,zeta_GP
    REAL                    :: xV(MESH%nVertexMax)
    REAL                    :: yV(MESH%nVertexMax)
    REAL                    :: zV(MESH%nVertexMax)
    CHARACTER (len=30000)   :: VariableList

    ! needed for local copying
    TYPE(tUnstructPoint), ALLOCATABLE       :: LocalRecPoint(:)

    !-----------------------------------------------------------!
    INTENT(IN)              :: BND, EQN, IO, MESH
    INTENT(INOUT)           :: DISC
    !-----------------------------------------------------------!

    ! initialization common to hdf5 ascii output
    CALL ini_common_fault_receiver(DISC, MESH)

    IF (DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output.EQ. .TRUE.) THEN
      CALL fill_variables_list(DISC, VariableList)

      ! After we got indices of the elements,
      ! we have to create the header for the VFiles and give filenames
      CALL create_file_headers(DISC, IO, MPI, VariableList)

      ! call modular subroutine for polynom evaluation
      ! at each receiver given
      CALL eval_faultreceiver(DISC%DynRup%DynRup_out_atPickpoint,MESH,BND,DISC)

      ! fills file header for each receiver with its parameters
      CALL write_header_info_to_files(EQN, MESH, DISC, IO, MPI)
    ENDIF

  END SUBROUTINE ini_fault_receiver


 !> modular subroutine for polynom evaluation at each receiver given
 !<
 SUBROUTINE eval_faultreceiver(DynRup_output,MESH,BND,DISC)
    !--------------------------------------------------------------------------!
    USE DGBasis_mod
    !--------------------------------------------------------------------------!
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tUnstructMesh)     :: MESH                   !< Mesh structure        !
    TYPE(tBoundary)         :: BND                    !< Boundar< structure    !
    TYPE(tDiscretization)   :: DISC                   !< Discretization struct.!
    !--------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    TYPE(tDynRup_output), target  :: DynRup_output                            !< Output data for Dynamic Rupture processes
    INTEGER                 :: i,l,iFault,iBndGP                              ! loop index
    INTEGER                 :: iElem, iSide, iNeighbor
    INTEGER                 :: MPIIndex, iObject
    REAL                    :: xV(MESH%nVertexMax)
    REAL                    :: yV(MESH%nVertexMax)
    REAL                    :: zV(MESH%nVertexMax)
    REAL                    :: io_x, io_y, io_z                               ! temp store of receiver location
    REAL                    :: xi,eta,zeta,xi_n,eta_n,zeta_n,xi_GP,eta_GP,zeta_GP
    REAL                    :: phi1,phi2
    REAL                    :: distance, shortest
    REAL                    :: chi,tau
    REAL                    :: tmp_mat(6)                       ! temporary stress tensors
    !-------------------------------------------------------------------------!
    INTENT(IN)              :: DISC, MESH, BND
    INTENT(INOUT)           :: DynRup_output
    !-------------------------------------------------------------------------!
    !
    ! output just for DynRup_output%nDR_pick at "+" elements
    ALLOCATE(DynRup_output%OutEval(DynRup_output%nDR_pick,1,DISC%Galerkin%nDegFr,2))
    ALLOCATE(DynRup_output%OutInt(DynRup_output%nDR_pick,1)) ! contains BndGP index
    ALLOCATE(DynRup_output%OutInt_dist(DynRup_output%nDR_pick)) ! saves distance to nears BndGP
    DynRup_output%OutEval       = 0.
    !
    ! Determine basis functions for DR pickpoint output and store them in OutEval
    DO i = 1,DynRup_output%nDR_pick
       !
       iFault              = DynRup_output%RecPoint(i)%index
       iElem               = MESH%Fault%Face(iFault,1,1)
       iSide               = MESH%Fault%Face(iFault,2,1)
       iNeighbor           = MESH%Fault%Face(iFault,1,2)          ! iNeighbor denotes "-" side
       !
       ! get the location in the reference element of Neighbor
       IF (iNeighbor.NE.0) THEN
           xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iNeighbor))
           yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iNeighbor))
           zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iNeighbor))
       ELSE
           ! The neighbor element belongs to a different  domain
           iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
           MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
           !
           xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
           yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
           zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
       ENDIF
       !
       xi                  = DynRup_output%RecPoint(i)%xi
       eta                 = DynRup_output%RecPoint(i)%eta
       zeta                = DynRup_output%RecPoint(i)%zeta
       io_x                = DynRup_output%RecPoint(i)%X
       io_y                = DynRup_output%RecPoint(i)%Y
       io_z                = DynRup_output%RecPoint(i)%Z
       !
       CALL TrafoXYZ2XiEtaZeta(xi_n,eta_n,zeta_n,io_x,io_y,io_z,xV,yV,zV,MESH%LocalVrtxType(iElem))
       !
       ! store basis functions
       DO l=1,DISC%Galerkin%nDegFr
          CALL BaseFunc3D(phi1,l,xi,eta,zeta,DISC%Galerkin%nPoly,       &
          DISC%Galerkin%cPoly3D_Tet,                                    &
          DISC%Galerkin%NonZeroCPoly_Tet,                               &
          DISC%Galerkin%NonZeroCPolyIndex_Tet                           )
          DynRup_output%OutEval(i,1,l,1) = phi1
          CALL BaseFunc3D(phi2,l,xi_n,eta_n,zeta_n,DISC%Galerkin%nPoly, &
          DISC%Galerkin%cPoly3D_Tet,                                    &
          DISC%Galerkin%NonZeroCPoly_Tet,                               &
          DISC%Galerkin%NonZeroCPolyIndex_Tet                           )
          DynRup_output%OutEval(i,1,l,2) = phi2
      ENDDO
      !
      ! Find nearest GP from which we will take MuVal, P_0, S_XY, S_XZ
      shortest = +2.0D0
      DO iBndGP = 1,DISC%Galerkin%nBndGP
         !
         chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
         tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
!~          write(6,*) chi, ' ', tau, ' ', DISC%Galerkin%bndGaussP_Tet(1,iBndGP), ' ', DISC%Galerkin%bndGaussP_Tet(2,iBndGP)
         CALL TrafoChiTau2XiEtaZeta(xi_GP,eta_GP,zeta_GP,chi,tau,iSide,0)
         !
         distance = SQRT((xi-xi_GP)**2 + (eta-eta_GP)**2 + (zeta-zeta_GP)**2)
         IF (distance .LT. shortest) THEN
            shortest = distance
            DynRup_output%OutInt_dist(i) = shortest
            DynRup_output%OutInt(i,1) = iBndGP
         ENDIF
      ENDDO ! iBndGP - to find nearest GP for DR pickpoint
      !
   ENDDO ! Determine basis functions for DR pickpoint output and store them in OutEval
  END SUBROUTINE eval_faultreceiver
 !
 !> modular subroutine for finding the middle point of an element
 !<
  SUBROUTINE findMiddlePoint( &
    point_1_x, &
    point_1_y, &
    point_1_z, &
    point_2_x, &
    point_2_y, &
    point_2_z, &
    point_3 )
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    REAL:: point_1_x,point_1_y,point_1_z
    REAL:: point_2_x,point_2_y,point_2_z
    REAL:: point_3(3)
    !-------------------------------------------------------------------------!
    !
    point_3(1)=0.5*point_1_x+0.5*point_2_x
    point_3(2)=0.5*point_1_y+0.5*point_2_y
    point_3(3)=0.5*point_1_z+0.5*point_2_z
    !
  END SUBROUTINE
 !
 !> recursive function for 2D refinement strategy 2 initialization
 !> element edge middle points are connected
 !<
  RECURSIVE FUNCTION refineFaultOutput_strategy2(MESH,LocalRecPoint,element_index,iFault,element_x,element_y,element_z,xp,yp,zp,xV,yV,zV,refinement,SubElem) &
                     result(local_index)
     !-------------------------------------------------------------------------!
     USE DGBasis_mod
     !-------------------------------------------------------------------------!
     ! Argument list declaration
     TYPE(tUnstructMesh)     :: MESH
     TYPE(tUnstructPoint), dimension(*) :: LocalRecPoint
     INTEGER,VALUE :: iFault
     INTEGER :: element_index
     REAL :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
     REAL :: element_x(MESH%GlobalElemType),element_y(MESH%GlobalElemType),element_z(MESH%GlobalElemType)
     REAL :: xV(MESH%GlobalElemType), yV(MESH%GlobalElemType), zV(MESH%GlobalElemType)
     !-------------------------------------------------------------------------!
     ! Local variable declaration
     REAL :: local_element_x(MESH%GlobalElemType),local_element_y(MESH%GlobalElemType),local_element_z(MESH%GlobalElemType)
     REAL :: local_global_element_x(MESH%GlobalElemType),local_global_element_y(MESH%GlobalElemType),local_global_element_z(MESH%GlobalElemType)
     REAL :: p1_local(3)
     REAL :: p2_local(3)
     REAL :: p3_local(3)
     REAL :: p1_global(3)
     REAL :: p2_global(3)
     REAL :: p3_global(3)
     INTEGER, VALUE :: refinement
     INTEGER :: local_index, j
     INTEGER, VALUE :: SubElem
     IF (refinement.EQ.0) THEN
        LocalRecPoint(iFault)%inside = .TRUE.        ! Point is located in an element
        LocalRecPoint(iFault)%index  = element_index ! number in Fault%Face list!
        LocalRecPoint(iFault)%globalreceiverindex = iFault
        LocalRecPoint(iFault)%X = (xp(1)+xp(2)+xp(3))/3.0
        LocalRecPoint(iFault)%Y = (yp(1)+yp(2)+yp(3))/3.0
        LocalRecPoint(iFault)%Z = (zp(1)+zp(2)+zp(3))/3.0
        DO j=1,3
        	LocalRecPoint(iFault)%coordx(j)=xp(j)
        	LocalRecPoint(iFault)%coordy(j)=yp(j)
        	LocalRecPoint(iFault)%coordz(j)=zp(j)
        END DO
	CALL TrafoXYZ2XiEtaZeta( LocalRecPoint(iFault)%xi, LocalRecPoint(iFault)%eta, LocalRecPoint(iFault)%zeta, &
                                 LocalRecPoint(iFault)%X, LocalRecPoint(iFault)%Y, LocalRecPoint(iFault)%Z, &
                                 xV,yV,zV,MESH%LocalVrtxType(element_index))
        local_index=iFault+1
        RETURN
     ELSE
        CALL findMiddlePoint(&
            element_x(1), element_y(1), element_z(1), &
            element_x(2), element_y(2), element_z(2), &
            p1_local &
        )
        CALL findMiddlePoint(&
            element_x(2), element_y(2), element_z(2), &
            element_x(3), element_y(3), element_z(3), &
            p2_local &
        )
        CALL findMiddlePoint(&
            element_x(3), element_y(3), element_z(3), &
            element_x(1), element_y(1), element_z(1), &
            p3_local &
        )
        CALL findMiddlePoint(&
            xp(1), yp(1), zp(1), &
            xp(2), yp(2), zp(2), &
            p1_global &
        )
        CALL findMiddlePoint(&
            xp(2), yp(2), zp(2), &
            xp(3), yp(3), zp(3), &
            p2_global &
        )
        CALL findMiddlePoint(&
            xp(3), yp(3), zp(3), &
            xp(1), yp(1), zp(1), &
            p3_global &
        )
        ! 1 traingle (1 - (1-2)/2 - (1-3)/2)
        local_element_x(1)=element_x(1)
        local_element_x(2)=p1_local(1)
        local_element_x(3)=p3_local(1)

        local_global_element_x(1)=xp(1)
        local_global_element_x(2)=p1_global(1)
        local_global_element_x(3)=p3_global(1)

        local_element_y(1)=element_y(1)
        local_element_y(2)=p1_local(2)
        local_element_y(3)=p3_local(2)

        local_global_element_y(1)=yp(1)
        local_global_element_y(2)=p1_global(2)
        local_global_element_y(3)=p3_global(2)

        local_element_z(1)=element_z(1)
        local_element_z(2)=p1_local(3)
        local_element_z(3)=p3_local(3)

        local_global_element_z(1)=zp(1)
        local_global_element_z(2)=p1_global(3)
        local_global_element_z(3)=p3_global(3)

        local_index = refineFaultOutput_strategy2(MESH,&
                                localRecPoint,&
                                element_index,&
                                iFault,&
                                local_element_x,&
                                local_element_y,&
                                local_element_z,&
                                local_global_element_x,&
                                local_global_element_y,&
                                local_global_element_z,&
                                xV,&
                                yV,&
                                zV,&
                                refinement-1,SubElem)
        ! (1-2)/2 - 2 - (2-3)/2
        local_element_x(1)=p1_local(1)
        local_element_x(2)=element_x(2)
        local_element_x(3)=p2_local(1)

        local_global_element_x(1)=p1_global(1)
        local_global_element_x(2)=xp(2)
        local_global_element_x(3)=p2_global(1)

        local_element_y(1)=p1_local(2)
        local_element_y(2)=element_y(2)
        local_element_y(3)=p2_local(2)

        local_global_element_y(1)=p1_global(2)
        local_global_element_y(2)=yp(2)
        local_global_element_y(3)=p2_global(2)

        local_element_z(1)=p1_local(3)
        local_element_z(2)=element_z(2)
        local_element_z(3)=p2_local(3)

        local_global_element_z(1)=p1_global(3)
        local_global_element_z(2)=zp(2)
        local_global_element_z(3)=p2_global(3)

        local_index = refineFaultOutput_strategy2(MESH,&
                                localRecPoint,&
                                element_index,&
                                local_index,&
                                local_element_x,&
                                local_element_y,&
                                local_element_z,&
                                local_global_element_x,&
                                local_global_element_y,&
                                local_global_element_z,&
                                xV,&
                                yV,&
                                zV,&
                                refinement-1,SubElem)
        ! (1-2)/2 - (2-3)/2 - (3-1)/2
        local_element_x(1)=p1_local(1)
        local_element_x(2)=p2_local(1)
        local_element_x(3)=p3_local(1)

        local_global_element_x(1)=p1_global(1)
        local_global_element_x(2)=p2_global(1)
        local_global_element_x(3)=p3_global(1)

        local_element_y(1)=p1_local(2)
        local_element_y(2)=p2_local(2)
        local_element_y(3)=p3_local(2)

        local_global_element_y(1)=p1_global(2)
        local_global_element_y(2)=p2_global(2)
        local_global_element_y(3)=p3_global(2)

        local_element_z(1)=p1_local(3)
        local_element_z(2)=p2_local(3)
        local_element_z(3)=p3_local(3)

        local_global_element_z(1)=p1_global(3)
        local_global_element_z(2)=p2_global(3)
        local_global_element_z(3)=p3_global(3)

        local_index = refineFaultOutput_strategy2(MESH,&
                                localRecPoint,&
                                element_index,&
                                local_index,&
                                local_element_x,&
                                local_element_y,&
                                local_element_z,&
                                local_global_element_x,&
                                local_global_element_y,&
                                local_global_element_z,&
                                xV,&
                                yV,&
                                zV,&
                                refinement-1,SubElem)

        ! (3-1)/2 - (2-3)/2 - 3
        local_element_x(1)=p3_local(1)
        local_element_x(2)=p2_local(1)
        local_element_x(3)=element_x(3)

        local_global_element_x(1)=p3_global(1)
        local_global_element_x(2)=p2_global(1)
        local_global_element_x(3)=xp(3)

        local_element_y(1)=p3_local(2)
        local_element_y(2)=p2_local(2)
        local_element_y(3)=element_y(3)

        local_global_element_y(1)=p3_global(2)
        local_global_element_y(2)=p2_global(2)
        local_global_element_y(3)=yp(3)

        local_element_z(1)=p3_local(3)
        local_element_z(2)=p2_local(3)
        local_element_z(3)=element_z(3)

        local_global_element_z(1)=p3_global(3)
        local_global_element_z(2)=p2_global(3)
        local_global_element_z(3)=zp(3)

        local_index = refineFaultOutput_strategy2(MESH,&
                                localRecPoint,&
                                element_index,&
                                local_index,&
                                local_element_x,&
                                local_element_y,&
                                local_element_z,&
                                local_global_element_x,&
                                local_global_element_y,&
                                local_global_element_z,&
                                xV,&
                                yV,&
                                zV,&
                                refinement-1,SubElem)
     !
     END IF
  END FUNCTION refineFaultOutput_strategy2
 !
 !> recursive function for 2D refinement strategy 1
 !> Barycenters of elements are connected
 !<
  RECURSIVE FUNCTION refineFaultOutput_strategy1(MESH,LocalRecPoint,iFault,element_x,element_y,element_z,xp,yp,zp,refinement,SubElem) &
                     result(local_index)
     !-------------------------------------------------------------------------!
     ! Argument list declaration
     TYPE(tUnstructMesh)                :: MESH
     TYPE(tUnstructPoint),dimension(*)  :: LocalRecPoint
     INTEGER, VALUE                     :: iFault
     REAL :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
     REAL :: element_x(MESH%GlobalElemType),element_y(MESH%GlobalElemType),element_z(MESH%GlobalElemType)
     !-------------------------------------------------------------------------!
     ! Local variable declaration
     REAL           :: local_element_x(MESH%GlobalElemType),local_element_y(MESH%GlobalElemType),local_element_z(MESH%GlobalElemType)
     REAL           :: local_global_element_x(MESH%GlobalElemType),local_global_element_y(MESH%GlobalElemType),local_global_element_z(MESH%GlobalElemType)
     INTEGER, VALUE :: refinement
     INTEGER        :: local_index, j
     INTEGER, VALUE :: SubElem
     !-------------------------------------------------------------------------!
     !
     IF (refinement.EQ.0) THEN
        LocalRecPoint(iFault)%xi   = (element_x(1)+element_x(2)+element_x(3))/3.0
        LocalRecPoint(iFault)%eta  = (element_y(1)+element_y(2)+element_y(3))/3.0
        LocalRecPoint(iFault)%zeta = (element_z(1)+element_z(2)+element_z(3))/3.0
        LocalRecPoint(iFault)%inside = .true. ! Point is located in an element
        LocalRecPoint(iFault)%index  = (iFault-1)/SubElem + 1! number in Fault%Face list!
        LocalRecPoint(iFault)%globalreceiverindex = iFault
        LocalRecPoint(iFault)%X = (xp(1)+xp(2)+xp(3))/3.0
        LocalRecPoint(iFault)%Y = (yp(1)+yp(2)+yp(3))/3.0
        LocalRecPoint(iFault)%Z = (zp(1)+zp(2)+zp(3))/3.0
        DO j=1,3
            LocalRecPoint(iFault)%coordx(j)=xp(j)
            LocalRecPoint(iFault)%coordy(j)=yp(j)
            LocalRecPoint(iFault)%coordz(j)=zp(j)
        END DO
        local_index=iFault+1
        RETURN
     ELSE
        !case 1 : 1,2,mid-point
        local_element_x(1)=element_x(1)
        local_element_x(2)=element_x(2)
        local_element_x(3)=(element_x(1)+element_x(2)+element_x(3))/3.0
        local_global_element_x(1)=xp(1)
        local_global_element_x(2)=xp(2)
        local_global_element_x(3)=(xp(1)+xp(2)+xp(3))/3.0

        local_element_y(1)=element_y(1)
        local_element_y(2)=element_y(2)
        local_element_y(3)=(element_y(1)+element_y(2)+element_y(3))/3.0
        local_global_element_y(1)=yp(1)
        local_global_element_y(2)=yp(2)
        local_global_element_y(3)=(yp(1)+yp(2)+yp(3))/3.0

        local_element_z(1)=element_z(1)
        local_element_z(2)=element_z(2)
        local_element_z(3)=(element_z(1)+element_z(2)+element_z(3))/3.0
        local_global_element_z(1)=zp(1)
        local_global_element_z(2)=zp(2)
        local_global_element_z(3)=(zp(1)+zp(2)+zp(3))/3.0

        local_index = refineFaultOutput_strategy1(MESH,&
                                localRecPoint,&
                                iFault,&
                                local_element_x,&
                                local_element_y,&
                                local_element_z,&
                                local_global_element_x,&
                                local_global_element_y,&
                                local_global_element_z,&
                                refinement-1,SubElem)
        !case 2 : mid-point,2,3(iFault-1)/SubElem + 1
        local_element_x(1)=(element_x(1)+element_x(2)+element_x(3))/3.0
        local_element_x(2)=element_x(2)
        local_element_x(3)=element_x(3)
        local_global_element_x(1)=(xp(1)+xp(2)+xp(3))/3.0
        local_global_element_x(2)=xp(2)
        local_global_element_x(3)=xp(3)

        local_element_y(1)=(element_y(1)+element_y(2)+element_y(3))/3.0
        local_element_y(2)=element_y(2)
        local_element_y(3)=element_y(3)
        local_global_element_y(1)=(yp(1)+yp(2)+yp(3))/3.0
        local_global_element_y(2)=yp(2)
        local_global_element_y(3)=yp(3)

        local_element_z(1)=(element_z(1)+element_z(2)+element_z(3))/3.0
        local_element_z(2)=element_z(2)
        local_element_z(3)=element_z(3)
        local_global_element_z(1)=(zp(1)+zp(2)+zp(3))/3.0
        local_global_element_z(2)=zp(2)
        local_global_element_z(3)=zp(3)
        local_index=refineFaultOutput_strategy1(MESH,&
                                localRecPoint,&
                                local_index,&
                                local_element_x,&
                                local_element_y,&
                                local_element_z,&
                                local_global_element_x,&
                                local_global_element_y,&
                                local_global_element_z,&
                                refinement-1,SubElem)
        !case 3 : 1,mid-point,3
        local_element_x(1)=element_x(1)
        local_element_x(2)=(element_x(1)+element_x(2)+element_x(3))/3.0
        local_element_x(3)=element_x(3)
        local_global_element_x(1)=xp(1)
        local_global_element_x(2)=(xp(1)+xp(2)+xp(3))/3.0
        local_global_element_x(3)=xp(3)

        local_element_y(1)=element_y(1)
        local_element_y(2)=(element_y(1)+element_y(2)+element_y(3))/3.0
        local_element_y(3)=element_y(3)
        local_global_element_y(1)=yp(1)
        local_global_element_y(2)=(yp(1)+yp(2)+yp(3))/3.0
        local_global_element_y(3)=yp(3)

        local_element_z(1)=element_z(1)
        local_element_z(2)=(element_z(1)+element_z(2)+element_z(3))/3.0
        local_element_z(3)=element_z(3)
        local_global_element_z(1)=zp(1)
        local_global_element_z(2)=(zp(1)+zp(2)+zp(3))/3.0
        local_global_element_z(3)=zp(3)
        local_index= refineFaultOutput_strategy1(MESH,&
                                localRecPoint,&
                                local_index,&
                                local_element_x,&
                                local_element_y,&
                                local_element_z,&
                                local_global_element_x,&
                                local_global_element_y,&
                                local_global_element_z,&
                                refinement-1,SubElem)

        RETURN
     END IF
  END FUNCTION refineFaultOutput_strategy1
 !
 !> subroutine handling recursive functions for 2D refinement strategy 1 and 2, initialization
 !<
  SUBROUTINE refineFaultOutput(strategy,MESH,element_index,LocalRecPoint,iFault,element_x,element_y,element_z,xp,yp,zp,xV,yV,zV,refinement,SubElem)
     !-------------------------------------------------------------------------!
     ! Argument list declaration
     INTEGER                            :: strategy
     TYPE(tUnstructMesh)                :: MESH
     TYPE(tUnstructPoint), DIMENSION(*) :: LocalRecPoint
     INTEGER, VALUE :: iFault
     INTEGER        :: element_index
     REAL           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
     REAL           :: xV(MESH%GlobalElemType), yV(MESH%GlobalElemType), zV(MESH%GlobalElemType)
     REAL           :: element_x(MESH%GlobalElemType),element_y(MESH%GlobalElemType),element_z(MESH%GlobalElemType)
     INTEGER, VALUE :: refinement
     INTEGER, VALUE :: SubElem
     !-------------------------------------------------------------------------!
     ! Local variable declaration
     INTEGER    :: r
     !-------------------------------------------------------------------------!
     !
     IF (strategy.EQ.1) THEN
        r=refineFaultOutput_strategy1(MESH,LocalRecPoint,iFault,element_x,element_y,element_z,xp,yp,zp,refinement,SubElem)
     ELSE
        r=refineFaultOutput_strategy2(MESH,LocalRecPoint,element_index,iFault,element_x,element_y,element_z,xp,yp,zp,xV,yV,zV,refinement,SubElem)
     ENDIF

  END SUBROUTINE refineFaultOutput
  !
  !> subroutine setting local coordinate system in each fault elemenet
  !<
  SUBROUTINE set_xi_eta_zeta(MESH,side_index,element_xi,element_eta,element_zeta)
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tUnstructMesh)     :: MESH
    integer                 :: side_index
    REAL                    :: element_xi(MESH%nVertexMax)
    REAL                    :: element_eta(MESH%nVertexMax)
    REAL                    :: element_zeta(MESH%nVertexMax)
    !-------------------------------------------------------------------------!
    !
    ! 1 - 3 - 2 face (zeta=0)
    IF (side_index.EQ.1) THEN
        element_xi(1) = 0.0
        element_xi(2) = 0.0
        element_xi(3) = 1.0
        element_eta(1) = 0.0
        element_eta(2) = 1.0
        element_eta(3) = 0.0
        element_zeta(1) = 0.0
        element_zeta(2) = 0.0
        element_zeta(3) = 0.0
    ! 1 - 2 - 4 face (eta =0)
    ELSE IF (side_index.EQ.2) THEN
        element_xi(1) = 0.0
        element_xi(2) = 1.0
        element_xi(3) = 0.0
        element_eta(1) = 0.0
        element_eta(2) = 0.0
        element_eta(3) = 0.0
        element_zeta(1) = 0.0
        element_zeta(2) = 0.0
        element_zeta(3) = 1.0
    ! 1 - 4 - 3  face (xi=0)
    ELSE IF (side_index.EQ.3) THEN
        element_xi(1) = 0.0
        element_xi(2) = 0.0
        element_xi(3) = 0.0
        element_eta(1) = 0.0
        element_eta(2) = 0.0
        element_eta(3) = 1.0
        element_zeta(1) = 0.0
        element_zeta(2) = 1.0
        element_zeta(3) = 0.0
    ! 2 - 3 - 4
    ELSE
        element_xi(1) = 1.0
        element_xi(2) = 0.0
        element_xi(3) = 0.0
        element_eta(1) = 0.0
        element_eta(2) = 1.0
        element_eta(3) = 0.0
        element_zeta(1) = 0.0
        element_zeta(2) = 0.0
        element_zeta(3) = 1.0
    ENDIF
  END SUBROUTINE set_xi_eta_zeta

  SUBROUTINE ini_fault_xdmfwriter(DISC,IO)
   use FaultWriter
  !-------------------------------------------------------------------------!
  ! Argument list declaration
  TYPE(tDiscretization)   :: DISC
  TYPE(tInputOutput)      :: IO
  !-------------------------------------------------------------------------!

   call initFaultOutput(DISC%DynRup%DynRup_out_elementwise%RecPoint, &
        DISC%DynRup%DynRup_out_elementwise%OutputMask, &
        DISC%DynRup%DynRup_out_elementwise%TmpState, &
        IO%OutputFile, &
        DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec, &
        IO%xdmfWriterBackend)

  END SUBROUTINE
!
!> Subroutine initializing the fault output
!<
  SUBROUTINE ini_fault_subsampled(EQN,MESH,BND,DISC,IO,MPI)
   USE DGBasis_mod
   USE create_fault_rotationmatrix_mod
   use FaultWriter
    !--------------------------------------------------------------------------!
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)        :: EQN                    ! Equation structure    !
    TYPE(tUnstructMesh)     :: MESH                   ! Mesh structure        !
    TYPE(tBoundary)         :: BND
    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    TYPE(tInputOutput)      :: IO                     ! IO structure          !
    TYPE(tMPI)              :: MPI                    ! MPI                   !
    !--------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                 :: i,j,k,l,vIt,in,iPick,iBndGP,iLocalNeighborSide,iOutPoints,iSubTet  ! loop index
    INTEGER                 :: iElem, iSide, iNeighbor, iFault
    INTEGER                 :: MPIIndex, iObject
    INTEGER                 :: OutVars
    INTEGER                 :: stat
    INTEGER                         :: VertexSide(4,3)
    REAL                    :: tolerance                     ! tolerance if the receiver belongs to this element
    REAL                    :: io_x, io_y, io_z              ! temp store of receiver location
    REAL                    :: xmin,xmax,ymin,ymax,zmin,zmax ! Maximum extend of computational domain
    REAL                    :: xi,eta,zeta,xi_n,eta_n,zeta_n,xi_GP,eta_GP,zeta_GP
    REAL                    :: chi,tau
    REAL                    :: xV(MESH%nVertexMax)
    REAL                    :: yV(MESH%nVertexMax)
    REAL                    :: zV(MESH%nVertexMax)
    REAL                    :: element_xi(MESH%nVertexMax)
    REAL                    :: element_eta(MESH%nVertexMax)
    REAL                    :: element_zeta(MESH%nVertexMax)
    REAL                    :: phi1,phi2
    REAL                    :: distance, shortest
    REAL                    :: rotmat(1:6,1:6)                  ! rotation matrix for individual fault receiver
    REAL                    :: tmp_mat(6)                       ! temporary stress tensors
    REAL                    :: x,y,z
    REAL                    :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
    LOGICAL                 :: exist
    INTEGER                 :: number_of_subtriangles, SubElem
    CHARACTER (len=200)     :: ptsoutfile
    CHARACTER (len=30000)   :: VariableList, VariableList_temp
    INTEGER                 :: allocStat
    !------------------------------iFault-------------------------------------------!
    INTENT(IN)              :: EQN, IO
    INTENT(INOUT)           :: DISC, MESH
    INTEGER                 :: r
    TYPE(tUnstructPoint), ALLOCATABLE       :: LocalRecPoint(:)    ! needed for local copying
    !------------------------------
    !
    !------------------------------
    !------------------------TrafoXi-------------------------------------------------!
    ! Prepare second header line
    !VariableList=TRIM('VARIABLES = "Time"')
    VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
    VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
    VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
    VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !
    IF(DISC%DynRup%DynRup_out_elementwise%refinement_strategy.eq.1) then
        number_of_subtriangles=3
    ELSE
        number_of_subtriangles=4
    ENDIF
    SubElem = number_of_subtriangles**DISC%DynRup%DynRup_out_elementwise%refinement
    logInfo(*) "Initialising Fault output. Refinement strategy: ",DISC%DynRup%DynRup_out_elementwise%refinement_strategy, &
               " Number of subtriangles: ",number_of_subtriangles
    ! ALLOCATE local list
    ALLOCATE(LocalRecPoint(MESH%Fault%nSide*(number_of_subtriangles**DISC%DynRup%DynRup_out_elementwise%refinement)))
    LocalRecPoint(:)%inside = .false.
    LocalRecPoint(:)%index = -1
    LocalRecPoint(:)%globalreceiverindex = -1
    LocalRecPoint(:)%X = 0.0D0
    LocalRecPoint(:)%Y = 0.0D0
    LocalRecPoint(:)%Z = 0.0D0
    LocalRecPoint(:)%coordx(1) = 0.0D0
    LocalRecPoint(:)%coordy(1) = 0.0D0
    LocalRecPoint(:)%coordz(1) = 0.0D0
    LocalRecPoint(:)%coordx(2) = 0.0D0
    LocalRecPoint(:)%coordy(2) = 0.0D0
    LocalRecPoint(:)%coordz(2) = 0.0D0
    LocalRecPoint(:)%coordx(3) = 0.0D0
    LocalRecPoint(:)%coordy(3) = 0.0D0
    LocalRecPoint(:)%coordz(3) = 0.0D0
    LocalRecPoint(:)%Y = 0.0D0
    LocalRecPoint(:)%Z = 0.0D0
    !
    in=0
    !
    DO iFault = 1, MESH%Fault%nSide
       iElem  = MESH%Fault%Face(iFault,1,1)
       !
       IF (iElem > 0) THEN
           xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
           yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
           zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
           iSide               = MESH%Fault%Face(iFault,2,1)
        CALL set_xi_eta_zeta(MESH,iSide,element_xi,element_eta,element_zeta)
        iLocalNeighborSide = MESH%Fault%Face(iFault,2,2)
        DO j=1,3
           xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
           yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
           zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
        END DO
        !
        !CALL refineFaultOutput(DISC%DynRup%DynRup_out_elementwise%refinement_strategy,MESH,iFault,LocalRecPoint, &
        !		       (iFault-1)*SubElem, &
        !                       element_xi,element_eta,element_zeta,xp,yp,zp,xV,yV,zV, &
        !                       DISC%DynRup%DynRup_out_elementwise%refinement,SubElem)
        CALL refineFaultOutput(DISC%DynRup%DynRup_out_elementwise%refinement_strategy,MESH,iFault,LocalRecPoint,(iFault-1)*(number_of_subtriangles**DISC%DynRup%DynRup_out_elementwise%refinement)+1,element_xi,element_eta,element_zeta,xp,yp,zp,xV,yV,zV,DISC%DynRup%DynRup_out_elementwise%refinement,SubElem)
        in=in+SubElem
        END IF
    ENDDO !iFault = 1, MESH%Fault%nSide
    !
    IF ( in.EQ.0 ) THEN
       ! If no pickpoint lies inside the domain, switch fault pickpoint output off.
       logInfo(*) 'No fault output receivers in this MPI domain '
       DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .FALSE.
       DISC%DynRup%DynRup_out_elementwise%nDR_pick = 0
       DISC%DynRup%DynRup_out_elementwise%nOutPoints = 0

       ALLOCATE( DISC%DynRup%DynRup_out_elementwise%RecPoint(0) )
       RETURN
    ELSE
       logInfo(*) 'Pick fault output at ',in,' points in this MPI domain.'
       DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .TRUE.
       DISC%DynRup%DynRup_out_elementwise%nDR_pick = in
       DISC%DynRup%DynRup_out_elementwise%nOutPoints = in
       !
       ALLOCATE( DISC%DynRup%DynRup_out_elementwise%RecPoint(DISC%DynRup%DynRup_out_elementwise%nDR_pick) )
       !
       in = 0
       DO i = 1,MESH%Fault%nSide*SubElem
         IF (LocalRecPoint(i)%inside) THEN
            in = in +1
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%X = LocalRecPoint(i)%X
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%Y = LocalRecPoint(i)%Y
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%Z = LocalRecPoint(i)%Z
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%index = LocalRecPoint(i)%index
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%globalreceiverindex = LocalRecPoint(i)%globalreceiverindex
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%xi = LocalRecPoint(i)%xi
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%eta = LocalRecPoint(i)%eta
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%zeta = LocalRecPoint(i)%zeta
            DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%inside = LocalRecPoint(i)%inside
            DO j=1,3
              DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%coordx(j)=LocalRecPoint(i)%coordx(j)
              DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%coordy(j)=LocalRecPoint(i)%coordy(j)
              DISC%DynRup%DynRup_out_elementwise%RecPoint(in)%coordz(j)=LocalRecPoint(i)%coordz(j)
            END DO
         ENDIF
       ENDDO !1,MESH%Fault%nSide*(number_of_subtriangles**DISC%DynRup%DynRup_out_elementwise%refinement)
       !
       ALLOCATE( DISC%DynRup%DynRup_Constants(DISC%DynRup%DynRup_out_elementwise%nDR_pick))
       ALLOCATE( DISC%DynRup%DynRup_Constants_GlobInd(MAXVAL(DISC%DynRup%DynRup_out_elementwise%RecPoint(:)%globalreceiverindex)))
       !
       DISC%DynRup%DynRup_out_elementwise%MaxPickStore = 50
       ! Deallocate global array
       ! options of outputted components
       OutVars = 0

       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(1).EQ.1) OutVars = OutVars + 2
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(2).EQ.1) OutVars = OutVars + 3
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(3).EQ.1) OutVars = OutVars + 1
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(4).EQ.1) OutVars = OutVars + 2
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(5).EQ.1) OutVars = OutVars + 3
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(6).EQ.1) OutVars = OutVars + 2
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(7).EQ.1) OutVars = OutVars + 1
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(8).EQ.1) OutVars = OutVars + 1
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(9).EQ.1) OutVars = OutVars + 1
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(10).EQ.1) OutVars = OutVars + 1
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(11).EQ.1) OutVars = OutVars + 1
       !
       ALLOCATE( DISC%DynRup%DynRup_out_elementwise%CurrentPick(DISC%DynRup%DynRup_out_elementwise%nDR_pick))
       ALLOCATE( DISC%DynRup%DynRup_out_elementwise%TmpTime(DISC%DynRup%DynRup_out_elementwise%MaxPickStore))
       ALLOCATE( DISC%DynRup%DynRup_out_elementwise%TmpState(DISC%DynRup%DynRup_out_elementwise%nDR_pick, 1, OutVars)) ! One TmpState is enough for the element wise output
       ALLOCATE( DISC%DynRup%DynRup_out_elementwise%OutVal(MAXVAL(DISC%DynRup%DynRup_out_elementwise%RecPoint(:)%globalreceiverindex),1,OutVars))
       ALLOCATE( DISC%DynRup%DynRup_out_elementwise%rotmat(DISC%DynRup%DynRup_out_elementwise%nDR_pick/SubElem,1:6,1:6)) ! store rotation matrix only for mother tets
       !
       ALLOCATE (DISC%DynRup%DynRup_out_elementwise%OutputLabel(OutVars))
       !
       k=1
       !
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(1).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(1) = 1
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(2) = 2
        k=k+2
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(2).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 3
        k=k+1
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 4
        k=k+1
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 5
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(3).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 6
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(4).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 7
        k=k+1
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 8
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(5).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 9
        k=k+1
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 10
        k=k+1
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 11
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(6).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 12
        k=k+1
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 13
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(7).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 14
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(8).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 15
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(9).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 16
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(10).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 17
        k=k+1
       ENDIF
       IF (DISC%DynRup%DynRup_out_elementwise%OutputMask(11).EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%OutputLabel(k) = 18
        k=k+1
       ENDIF
       !
       DISC%DynRup%DynRup_out_elementwise%CurrentPick(:)= 0.
       DISC%DynRup%DynRup_out_elementwise%OutVal        = 0.
       DISC%DynRup%DynRup_out_elementwise%rotmat        = 0.
       !
       CALL eval_faultreceiver(DISC%DynRup%DynRup_out_elementwise,MESH,BND,DISC)
       ! loop to evaluate the intial stress field in every subtet
       DO iOutPoints = 1,DISC%DynRup%DynRup_out_elementwise%nDR_pick ! sub tet indices
         iBndGP = DISC%DynRup%DynRup_out_elementwise%OutInt(iOutPoints,1)
         iFault = DISC%DynRup%DynRup_out_elementwise%RecPoint(iOutPoints)%index
         ! create fault rotationmatrix only once for every mother tet
         !
         IF (MOD((iOutPoints-1),SubElem).EQ.0) CALL create_fault_rotationmatrix(rotmat,iFault,EQN,MESH) ! mother index
         !
         ! store fault rotationmatrix only for every mother tet
         DISC%DynRup%DynRup_out_elementwise%rotmat((iOutPoints-1)/(SubElem)+1,1:6,1:6) = rotmat
         !
         tmp_mat = 0.0D0
         tmp_mat(1) = EQN%IniBulk_xx(iBndGP,iFault)
         tmp_mat(2) = EQN%IniBulk_yy(iBndGP,iFault)
         tmp_mat(3) = EQN%IniBulk_zz(iBndGP,iFault)
         tmp_mat(4) = EQN%IniShearXY(iBndGP,iFault)
         tmp_mat(5) = EQN%IniShearYZ(iBndGP,iFault)
         tmp_mat(6) = EQN%IniShearXZ(iBndGP,iFault)
         ! rotate into fault system
         tmp_mat(:) = MATMUL(rotmat(1:6,1:6),tmp_mat(:))
         !
         DISC%DynRup%DynRup_Constants(iOutPoints)%p0 =  tmp_mat(1)
         DISC%DynRup%DynRup_Constants(iOutPoints)%ts0 = tmp_mat(4)
         DISC%DynRup%DynRup_Constants(iOutPoints)%td0 = tmp_mat(6)
       END DO
   ENDIF
  !
  END SUBROUTINE ini_fault_subsampled
!
END MODULE ini_faultoutput_mod
