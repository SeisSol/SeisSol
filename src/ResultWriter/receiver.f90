!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

MODULE receiver_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE ini_receiver
     MODULE PROCEDURE ini_receiver
  END INTERFACE

  INTERFACE PGM_output
     MODULE PROCEDURE PGM_output
  END INTERFACE

  INTERFACE receiver
     MODULE PROCEDURE receiver
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC  :: receiver
  PUBLIC  :: ini_receiver
  PUBLIC  :: PGM_output
  !----------------------------------------------------------------------------

CONTAINS

  !< This routine finds the elements containing record points and creates the header for output files but doesn't enter values.
  !< This is the job of the sub receiver, which passes the values at the points to the output files for the different points.                     !
  SUBROUTINE ini_receiver(EQN,MESH,DISC,SOURCE,IO,MPI)
    !--------------------------------------------------------------------------
    USE DGBasis_mod
    USE common_receiver_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------
    TYPE (tEquations)        :: EQN
    TYPE (tUnstructMesh)     :: MESH
    TYPE (tDiscretization)   :: DISC
    TYPE (tSource)           :: SOURCE
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    !--------------------------------------------------------------------------
    ! local Variables
    !--------------------------------------------------------------------------
    CHARACTER(len=2)         :: VName(EQN%nVar) 
    CHARACTER(len=2)         :: VNameRot(EQN%Dimension)
    CHARACTER(len=3)         :: VNameMTC(9)
    CHARACTER(len=2)         :: VNameCD(4)
    !--------------------------------------------------------------------------
    INTEGER                  :: i, j, h, l, n, ielem, iVar, index_unitnumber, in, iPick, iCPU, iRot
    INTEGER                  :: Vel_first, Vel_last
    INTEGER                  :: AllocStat, stat, nFiles
    INTEGER                  :: nContinuousOutputRecordPoints
    CHARACTER (len=200)      :: ptsoutfile
    CHARACTER (len=30000)    :: VariableList, VariableList_temp
    LOGICAL                  :: exist
    CHARACTER (LEN=5)        :: cmyrank
    !--------------------------------------------------------------------------
    PARAMETER (nContinuousOutputRecordPoints = 50)
    !--------------------------------------------------------------------------
    INTENT(INOUT)            :: IO
    !--------------------------------------------------------------------------
    !
    ! Define required output pickmask
    EQN%nVar_fld = 0
    DO iVar = 1, EQN%nVar
      IF(IO%Outputmask(EQN%Dimension+iVar)) THEN
        EQN%nVar_fld = EQN%nVar_fld + 1
      ENDIF
    ENDDO
    ALLOCATE(IO%pickmask(EQN%nVar_fld))
    EQN%nVar_fld = 0
    DO iVar = 1, EQN%nVar
      IF(IO%Outputmask(EQN%Dimension+iVar)) THEN
        EQN%nVar_fld = EQN%nVar_fld + 1
        IO%pickmask(EQN%nVar_fld) = iVar
      ENDIF
    ENDDO
    !
    ! rotation
    iRot = size(IO%RotationMask)
    IF(IO%Rotation.NE.0) THEN
        EQN%nVar_rot = 0
        DO iVar = 1,iRot
          IF(IO%RotationMask(iVar)) THEN
            EQN%nVar_rot = EQN%nVar_rot + 1
          ENDIF
        ENDDO
        ALLOCATE(IO%pickmask_rot(EQN%nVar_rot))
        EQN%nVar_rot = 0
        DO iVar = 1, EQN%Dimension
          IF(IO%RotationMask(iVar)) THEN
            EQN%nVar_rot = EQN%nVar_rot + 1
            IO%pickmask_rot(EQN%nVar_rot) = iVar
          ENDIF
        ENDDO
    ENDIF

    ! Define header in receiver files
    IF(EQN%Poroelasticity.EQ.0) THEN
           VName= (/'xx','yy','zz','xy','yz','xz','u ','v ','w '/)
           VNameRot= (/'rx','ry','rz'/)
           VNameMTC= (/'Gxx','Gxy','Gxz','Gyx','Gyy','Gyz','Gzx','Gzy','Gzz'/)
           VNameCD= (/'cx','cy','cz','dv'/)
    ELSE
           VName= (/'xx','yy','zz','xy','xz','yz','u ','v ','w ','p ','uf','vf','wf'/)
           VNameRot= (/'rx','ry','rz'/)

    ENDIF
    !
    CALL common_ini_receiver(EQN,MESH,DISC,SOURCE,IO,MPI)
    !                                                                          
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<       Coordinate of picking location                    >'
    logInfo(*) '<--------------------------------------------------------->'
    !                                                                          
    IF ( IO%nlocalRecordPoint.EQ.0 ) THEN
        ! If no receiver lies inside the domain, tell the code the modified number.
        logInfo(*) 'No temporal signal or PGM is picked'
        IO%nRecordPoint      = 0 
        IO%ntotalRecordPoint = 0                     
    ELSE                                                                       
        logInfo(*) 'Pick temporal signal at ',IO%nRecordPoint,' points.'
        !Allocate space for the VFiles                                                 
        nFiles=IO%nRecordPoint
        IO%w_or_wout_ana=9

        !  
        ALLOCATE( IO%UNIT%VFile(nFiles) ) 
        !Giving the VFiles unit file numbers 
        DO i = 1, nFiles                                                      
          IO%UNIT%VFile(i) = IO%UNIT%receiverStart + (i-1)
          IF (IO%UNIT%VFile(i).GT.IO%UNIT%maxThisDom) THEN                     
             logError(*) 'maximum Unit numbers reserved for receiver files is ',IO%UNIT%maxThisDom
             STOP                                                              
          END IF                                                               
        END DO                                                                  
        !                                                              
        !After we got the indices of the elements, we have to create the header for the VFiles and give filenames 
        !for saving picked temporal values of the variables  
        !
        DO i = 1, IO%nRecordPoint
             !
             IF(.NOT.IO%UnstructRecPoint(i)%inside) THEN
                CYCLE    ! If receiver is NOT in the (sub-)domain, continue.
             ENDIF
             !
#ifdef PARALLEL
             WRITE(cmyrank,'(I5.5)') MPI%myrank                                   ! myrank -> cmyrank
             WRITE(ptsoutfile, '(a,a10,i5.5,a1,a5,a4)') TRIM(IO%OutputFile),'-receiver-',i,'-',TRIM(cmyrank),'.dat'!
#else
             WRITE(ptsoutfile, '(a,a10,i5.5,a4)') TRIM(IO%OutputFile),'-receiver-',i,'.dat'
#endif
             logInfo(*) '    ... open file ', TRIM(ptsoutfile),' to save time signal'
             logInfo(*) '<--------------------------------------------------------->'
             logInfo(*) ' '
             !
             index_unitnumber = i
             !                                                                      
             INQUIRE(FILE = ptsoutfile, EXIST = exist)
             IF(exist) THEN
                ! If file exists, then append data
                OPEN(UNIT     = IO%UNIT%VFILE(index_unitnumber)                  , & !
                     FILE     = ptsoutfile                                       , & !
                     FORM     = 'FORMATTED'                                      , & !                         
                     STATUS   = 'OLD'                                            , & !
                     POSITION = 'APPEND'                                         , & !
                     RECL     = 80000                                            , & !
                     IOSTAT = stat                                                 ) !
                IF(stat.NE.0) THEN                                                   !
                   logError(*) 'cannot open ',ptsoutfile          !
                   logError(*) 'Error status: ', stat
                   STOP                                                              !
                END IF                                                               !
                CLOSE( IO%UNIT%VFILE(index_unitnumber) )
             ELSE
                ! If file does not exist, then write header
                OPEN(UNIT     = IO%UNIT%VFile(index_unitnumber)                  , & !
                     FILE     = ptsoutfile                                       , & !
                     FORM     = 'FORMATTED'                                      , & !                         
                     STATUS   = 'NEW'                                            , & !
                     RECL     = 80000                                            , & !
                     IOSTAT = stat                                                 ) !
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
                WRITE(IO%UNIT%VFile(index_unitnumber),'(a,I8.8,a)') & 
                    'TITLE = "Temporal Signal for receiver number ', i, ' "'
                !
                ! Second line of the header (variable names) 
                !               
                VariableList=TRIM('VARIABLES = "Time"') 
                !
                DO h=1,EQN%nVar
                   IF(IO%OutputMask(EQN%Dimension+h))THEN 
                      VariableList_temp=TRIM(VariableList)
                      WRITE(VariableList,'(a,a2,a,a1)')                   &
                         TRIM(VariableList_temp),',"',TRIM(VName(h)),'"'   
                   ENDIF
                ENDDO
                !  
                ! Rotational output
                IF(IO%Rotation.EQ.1) THEN
                    Vel_first = 7
                    Vel_last  = 9
                    j=0
                    DO h=Vel_first,Vel_last
                      j=j+1
                      IF(IO%OutputMask(EQN%Dimension+h))THEN 
                        VariableList_temp=TRIM(VariableList)
                        WRITE(VariableList,'(a,a2,a,a1)')                   &
                           TRIM(VariableList_temp),',"',TRIM(VNameRot(j)),'"'   
                      ENDIF
                    ENDDO
                ENDIF
                ! Moment Tensor output
                IF(IO%Rotation.EQ.2) THEN
                    DO h=1,9
                        VariableList_temp=TRIM(VariableList)
                        WRITE(VariableList,'(a,a2,a,a1)')                   &
                           TRIM(VariableList_temp),',"',TRIM(VNameMTC(h)),'"'   
                    ENDDO
                ENDIF
                ! Curl/Divergence output
                IF(IO%Rotation.EQ.3) THEN
                    DO h=1,4
                        VariableList_temp=TRIM(VariableList)
                        WRITE(VariableList,'(a,a2,a,a1)')                   &
                           TRIM(VariableList_temp),',"',TRIM(VNameCD(h)),'"'
                    ENDDO
                ENDIF
                !
                WRITE(IO%UNIT%VFile(index_unitnumber),'(a)') TRIM(VariableList)          
                !
                ! Third-fifth line of the header (comment, containing x positions) 
                !
                WRITE(IO%UNIT%VFile(index_unitnumber),'(a4,e25.12)') TRIM('# x1'), IO%UnstructRecPoint(i)%X
                WRITE(IO%UNIT%VFile(index_unitnumber),'(a4,e25.12)') TRIM('# x2'), IO%UnstructRecPoint(i)%Y
                WRITE(IO%UNIT%VFile(index_unitnumber),'(a4,e25.12)') TRIM('# x3'), IO%UnstructRecPoint(i)%Z
                !
                CLOSE( IO%UNIT%VFile(index_unitnumber) )
                !
             ENDIF
             !
        ENDDO
        !       
    ENDIF 
    !
    IF(IO%ntotalRecordPoint.NE.0) THEN
        IF(DISC%Galerkin%DGMethod.EQ.3) THEN
            ALLOCATE( IO%LocalPickTime(IO%ntotalRecordPoint) ) 
            IO%LocalPickTime(:) = 0.
        ENDIF
    ENDIF
    !
    !-------------------------------------
    ! Initialization of large output case
    !-------------------------------------
    IF(IO%nTotalRecordPoint.LT.nContinuousOutputRecordPoints) THEN
      IO%PickLarge    = 0  !Standard, opening and closing files at each output level
    ELSE
      IO%PickLarge    = 1    !Intermediate storage of output to reduce IO operations
      if (IO%pickDtType .eq. 2) then
            ! in 20% intervals of the maximum timesteps or every 500 levels
            IO%MaxPickStore = MIN(CEILING(DISC%MaxIteration/IO%pickdt*0.2, 8), 500_8)
      else
            ! in 10% intervals of the total simulate time or every 500 levels
            IO%MaxPickStore = MIN(CEILING(DISC%EndTime/IO%Pickdt*0.1, 8),500_8)
      endif

      ALLOCATE( IO%CurrentPick( IO%nRecordPoint))
      ALLOCATE( IO%TmpTime( IO%nRecordPoint,IO%MaxPickStore))
      IO%CurrentPick(:) = 0.

      j = 0
      DO h = 1,EQN%nVar
        IF(IO%OutputMask(EQN%Dimension+h)) j=j+1
      ENDDO
      ALLOCATE( IO%TmpState(IO%nRecordPoint,IO%MaxPickStore,j))
      
      IF(IO%Rotation.EQ.1) THEN
          Vel_first = 7
          Vel_last  = 9
          j = 0
          DO h = Vel_first,Vel_last
            IF(IO%OutputMask(EQN%Dimension+h)) j=j+1
          ENDDO
        !
        ALLOCATE( IO%TmpState_rot(IO%nRecordPoint,IO%MaxPickStore,j))
      ENDIF

      IF(IO%Rotation.EQ.2) THEN
        ALLOCATE( IO%TmpState_rot(IO%nRecordPoint,IO%MaxPickStore,9))
      ENDIF
      
      IF(IO%Rotation.EQ.3) THEN
        ALLOCATE( IO%TmpState_rot(IO%nRecordPoint,IO%MaxPickStore,4))
      ENDIF
        
    ENDIF   
    !------------------------------------------------
    ! Initialization of maximum values for PGM output
    !------------------------------------------------               
    !
    IF(IO%nPGMRecordPoint.GT.0)THEN
        !Check, which output is required (PGD, PGV, PGA)
        ALLOCATE( IO%PGM(IO%nPGMRecordPoint)      ) !PGV
        !ALLOCATE( IO%PGD_tmp(IO%nPGMRecordPoint,3)) !PGD
        !ALLOCATE( IO%PGA_tmp(IO%nPGMRecordPoint,3)) !PGA
        
        IO%PGM     = 0.
        !IO%PGD_tmp = 0.
        !IO%PGA_tmp = 0.

    ENDIF

  END SUBROUTINE ini_receiver


  SUBROUTINE PGM_output(IO,MPI)
    !--------------------------------------------------------------------------
    ! PGM_output writes the 
    ! Peak Ground Displacement (PGD),
    ! Peak Ground Velocity (PGV),
    ! Peak Ground Acceleration (PGA),
    ! into a file together with the xyz-coordinates of the registration point.
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    ! local Variables
    INTEGER                  :: i, iErr
    INTEGER                  :: allocstat
    !--------------------------------------------------------------------------
    INTENT(INOUT)            :: IO
    !--------------------------------------------------------------------------
    !
    ! output peak ground motion data
    IF(IO%nPGMRecordPoint.GT.0)THEN
        logInfo(*) 'Writing PGM data to file PGM_map.dat ... !'

#ifdef PARALLEL
        !IF(MPI%myrank.EQ.0)THEN
           ALLOCATE(MPI%PGMarray(IO%nPGMRecordPoint,MPI%nCPU),  &
                    STAT = allocStat                              )
           IF (allocStat .NE. 0) THEN
                 logError(*) 'could not allocate',&
                 ' PGMarray for PGM output! '
                 STOP
           END IF
           MPI%PGMarray = 0.
        !ENDIF

        !Collect PGM data from all CPUs in a common array MPI%PGMarray
        CALL MPI_GATHER(IO%PGM,IO%nPGMRecordPoint,MPI%MPI_AUTO_REAL,MPI%PGMarray(:,:),IO%nPGMRecordPoint,MPI%MPI_AUTO_REAL,0,MPI%commWorld,iErr)
        
        !CPU 0 is outputting the collected PGM data
        IF(MPI%myrank.EQ.0)THEN
            MPI%PGMarray(:,1) = MAXVAL(MPI%PGMarray,DIM = 2)
            OPEN(UNIT=999,FILE='PGM_map.dat',FORM='FORMATTED',RECL=500)
            WRITE(999,'(a75)')'x                       y                       z                       PGV'
            DO i = 1, IO%nPGMRecordPoint
                WRITE(999,*) IO%UnstructRecPoint(IO%PGMstartindex+i-1)%X,    &
                             IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Y,    &
                             IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Z,    &
                             MPI%PGMarray(i,1)
            ENDDO
            CLOSE(999)
            DEALLOCATE(MPI%PGMarray)
        ENDIF
#else   
        OPEN(UNIT=999,FILE='PGM_map.dat',FORM='FORMATTED',RECL=500)
        WRITE(999,'(a75)')'x                       y                       z                       PGV'
        DO i = 1, IO%nPGMRecordPoint
           WRITE(999,*) IO%UnstructRecPoint(IO%PGMstartindex+i-1)%X,    &
                        IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Y,    &
                        IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Z,    &
                        IO%PGM(i)
        ENDDO
        CLOSE(999)
#endif  
        
        DEALLOCATE(IO%PGM)
        !DEALLOCATE(IO%PGD_tmp)
        !DEALLOCATE(IO%PGA_tmp)

    ENDIF
    ! 
  END SUBROUTINE PGM_output
  !
  !
#ifdef GENERATEDKERNELS
  SUBROUTINE receiver(i_fullUpdateTime, i_timeStepWidth, i_receiverTime, i_numberOfReceivers, i_receiverIds, EQN,MESH,DISC,MPI,IO,time_op,dt_op)
#else
  SUBROUTINE receiver(EQN,MESH,DISC,MPI,IO,time_op,dt_op)
#endif
    !--------------------------------------------------------------------------
    USE common_receiver_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! argument list declaration
    !--------------------------------------------------------------------------
#ifdef GENERATEDKERNELS
    real*8  :: i_fullUpdateTime
    real*8  :: i_timeStepWidth
    real*8  :: i_receiverTime
    integer :: i_numberOfReceivers
    integer :: i_receiverIds(:)
#endif
    TYPE (tEquations)              :: EQN
    TYPE (tUnstructMesh)           :: MESH 
    TYPE (tDiscretization)         :: DISC
    TYPE (tMPI)                    :: MPI
    TYPE (tInputOutput)            :: IO
    REAL, OPTIONAL                 :: time_op, dt_op
    !--------------------------------------------------------------------------
    ! local Variables
    !--------------------------------------------------------------------------
    REAL                           :: time, dt
    REAL                           :: state( EQN%nVarTotal), state_rot(9)
    REAL                           :: stateToWrite( EQN%nVarTotal ), stateToWrite_rot(9)  ! local variables to avoid fortran runtime warnings when writing output
    REAL                           :: PGD_tot, PGV_tot, PGA_tot
    REAL                           :: TaylorDOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal,0:DISC%Galerkin%nPolyRec)
    REAL                           :: localpicktime
    INTEGER                        :: allocstat, stat
    INTEGER                        :: i, j, k, n, l, h, g, m
    CHARACTER (len=200)            :: ptsoutfile
    CHARACTER (LEN=5)              :: cmyrank
#ifdef GENERATEDKERNELS
    integer :: l_receiver
#endif
    !--------------------------------------------------------------------------
    INTENT(IN)                     :: EQN,MESH                                !
    INTENT(INOUT)                  :: IO,DISC                                 !
    !--------------------------------------------------------------------------
    !
    ! register epik/scorep function receiver
    EPIK_FUNC_REG("receiver")
    SCOREP_USER_FUNC_DEFINE()
    !--------------------------------------------------------------------------
    ! start epik/scorep function receiver
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN("receiver")
    !
    localpicktime = 0
    stateToWrite(:) = 0.0
    stateToWrite_rot(:) = 0.0
#ifdef GENERATEDKERNELS
    do l_receiver =1,i_numberOfReceivers
      j = i_receiverIds(l_receiver)
#else
    ! loop over global receivers
    DO j = 1,IO%ntotalRecordPoint
#endif
        !
        IF(.NOT.IO%UnstructRecpoint(j)%inside) THEN
           CYCLE    ! If receiver is NOT in the (sub-)domain, continue.
        ENDIF
        !
        IF(j.LE.IO%nRecordPoint)THEN
#ifdef PARALLEL
           WRITE(cmyrank,'(I5.5)') MPI%myrank                                   ! myrank -> cmyrank
           WRITE(ptsoutfile, '(a,a10,i5.5,a1,a5,a4)') TRIM(IO%OutputFile),'-receiver-',j,'-',TRIM(cmyrank),'.dat'!
#else
           WRITE(ptsoutfile, '(a,a10,i5.5,a4)') TRIM(IO%OutputFile),'-receiver-',j,'.dat'
#endif
        ENDIF
        !
        SELECT CASE(DISC%Galerkin%DGMethod)
        CASE(3)
          ! LTS
          CALL common_receiver_ck(EQN,MESH,DISC,IO,j,TaylorDof,dt,time,localpicktime)
        CASE DEFAULT
          ! GTS
          CALL common_receiver_ck(EQN,MESH,DISC,IO,j,TaylorDof,dt,time,localpicktime,dt_op,time_op)
        END SELECT
        !
#ifdef GENERATEDKERNELS
        time          = i_fullUpdateTime
        dt            = i_timeStepWidth
        localpicktime = i_receiverTime
#endif
        DO WHILE( (localpicktime.GE.time).AND.(localpicktime.LE.time+dt+1e-10).AND.(localpicktime.LE.DISC%EndTime+1e-10) )
            !
            CALL common_receiver_interp(EQN,MESH,DISC,IO,j,TaylorDof,time,localpicktime,state,state_rot)
            !
            IF(j.LE.IO%nRecordPoint)THEN   !only output subsequent samples for seismogram record points

                IF(IO%PickLarge.EQ.0) THEN
                  !
                  OPEN(UNIT     = IO%UNIT%VFILE(j)                                 , & !
                       FILE     = ptsoutfile                                       , & !
                       FORM     = 'FORMATTED'                                      , & !                         
                       STATUS   = 'OLD'                                            , & !
                       POSITION = 'APPEND'                                         , & !
                       RECL     = 80000                                            , & !
                       IOSTAT = stat                                                 ) !
                  IF(stat.NE.0) THEN                                                   !
                    logError(*) 'cannot open ',ptsoutfile          !
                    logError(*) 'Error status: ', stat
                    STOP                                                              !
                  END IF                                                               !
                  !
                  stateToWrite(1:size(IO%pickmask)) = State(IO%pickmask(:))
                  IF(IO%Rotation.EQ.0) THEN
                    WRITE(IO%UNIT%VFile(j),*) localpicktime, stateToWrite(1:size(IO%pickmask))
                  ELSE
                    stateToWrite_rot(1:size(IO%pickmask_rot)) = State_rot(IO%pickmask_rot(:))
                    WRITE(IO%UNIT%VFile(j),*) localpicktime, stateToWrite(1:size(IO%pickmask)), stateToWrite_rot(1:size(IO%pickmask_rot))
                  ENDIF

                  CLOSE( IO%UNIT%VFILE(j) )

                ELSEIF(IO%PickLarge.EQ.1) THEN
                  !@TODO we might want to fix frontran runtime warnings here as well!
                  IO%CurrentPick(j) = IO%CurrentPick(j) +1
                  IO%TmpTime(j,IO%CurrentPick(j)) = localpicktime
                  IO%TmpState(j,IO%CurrentPick(j),:) = State(IO%pickmask(:))
                  IF(IO%Rotation.NE.0) IO%TmpState_rot(j,IO%CurrentPick(j),:) = State_rot(IO%pickmask_rot(:))
              
                  IF(IO%CurrentPick(j).GE.IO%MaxPickStore.OR.ABS(DISC%EndTime-localpicktime)/DISC%EndTime.LE.1.0e-8) THEN
                    !
                    OPEN(UNIT     = IO%UNIT%VFILE(j)                                 , & !
                         FILE     = ptsoutfile                                       , & !
                         FORM     = 'FORMATTED'                                      , & !                         
                         STATUS   = 'OLD'                                            , & !
                         POSITION = 'APPEND'                                         , & !
                         RECL     = 80000                                            , & !
                         IOSTAT = stat                                                 ) !
                    IF(stat.NE.0) THEN                                                   !
                      logError(*) 'cannot open ',ptsoutfile           !
                      logError(*) 'Error status: ', stat
                      STOP                                                               !
                    END IF                                                               !
                    !
                    IF(IO%Rotation.EQ.0) THEN
                      DO k=1,IO%CurrentPick(j)
                        WRITE(IO%UNIT%VFile(j),*) IO%TmpTime(j,k), IO%TmpState(j,k,:)
                      ENDDO
                    ELSE
                      DO k=1,IO%CurrentPick(j)
                        WRITE(IO%UNIT%VFile(j),*) IO%TmpTime(j,k), IO%TmpState(j,k,:), IO%TmpState_rot(j,k,:)
                      ENDDO
                    ENDIF

                    IO%CurrentPick(j) = 0                

                    CLOSE( IO%UNIT%VFILE(j) )   
                
                  ENDIF               
              
                ENDIF !PickLarge
            
            ELSE  ! check PGM record point 

                k = j-IO%nRecordPoint
                ! First compute PGD as integral of velocities
                !IO%PGD_tmp(k,1) = IO%PGD_tmp(k,1) + state(7)*IO%pickdt
                !IO%PGD_tmp(k,2) = IO%PGD_tmp(k,2) + state(8)*IO%pickdt
                !IO%PGD_tmp(k,3) = IO%PGD_tmp(k,3) + state(9)*IO%pickdt
                !PGD_tot      = SQRT(IO%PGD_tmp(k,1)**2 + IO%PGD_tmp(k,2)**2 + IO%PGD_tmp(k,3)**2)
                !IF( PGD_tot.GT.IO%PGM(k,1) )THEN
                !   IO%PGM(k,1) = PGD_tot
                !ENDIF

                ! Second compute PGV 
                PGV_tot = SQRT(state(7)**2+state(8)**2+state(9)**2)
                IF( PGV_tot.GT.IO%PGM(k) )THEN
                   IO%PGM(k) = PGV_tot
                ENDIF

                ! Third compute PGA as derivative of velocities
                !IO%PGA_tmp(k,1) = (state(7)-IO%PGA_tmp(k,1)) / IO%pickdt
                !IO%PGA_tmp(k,2) = (state(8)-IO%PGA_tmp(k,2)) / IO%pickdt
                !IO%PGA_tmp(k,3) = (state(9)-IO%PGA_tmp(k,3)) / IO%pickdt
                !PGA_tot      = SQRT(IO%PGA_tmp(k,1)**2 + IO%PGA_tmp(k,2)**2 + IO%PGA_tmp(k,3)**2)
                !IF( PGA_tot.GT.IO%PGM(k,3) )THEN
                !   IO%PGM(k,3) = PGA_tot
                !ENDIF
                ! Store current velocity values for derivation of PGA in the next time step 
                !IO%PGA_tmp(k,1) = state(7)
                !IO%PGA_tmp(k,2) = state(8)
                !IO%PGA_tmp(k,3) = state(9)

            ENDIF  ! seismogram or PGM record points

            if (IO%pickDtType .eq. 2) then
                localpicktime = localpicktime + IO%pickdt * dt
            else
                localpicktime = localpicktime + IO%pickdt
            endif

        ENDDO ! END DO WHILE localpicktime
        !
        IF( DISC%Galerkin%DGMethod.EQ.3) THEN
          IO%LocalPickTime(j) = localpicktime
        ENDIF
        !
    ENDDO ! END DO IO%ntotalRecordPoint
    !
    IF( DISC%Galerkin%DGMethod.EQ.1) THEN
      IO%picktime = localpicktime
    ENDIF
    ! end epik/scorep function receiver
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()

  END SUBROUTINE receiver

END MODULE receiver_mod
