!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Michael Dumbser (michael.dumbser AT unitn.it, https://www5.unitn.it/People/it/Web/Persona/PER0029602#INFO)
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2007-2016, SeisSol Group
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

#include <Initializer/preProcessorMacros.fpp>

MODULE ini_SeisSol_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE ini_SeisSol
     MODULE PROCEDURE ini_SeisSol
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC  :: ini_SeisSol
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE ini_SeisSol(time,timestep,pvar,cvar,EQN,IC,MESH,MPI,      &
       SOURCE,DISC,BND,OptionalFields,IO, &
       programTitle) !
    !--------------------------------------------------------------------------
    USE COMMON_readpar_mod
    USE ini_OptionalFields_mod
    USE dg_setup_mod
    USE ini_MODEL_mod
 !   USE DGSponge_mod
    USE calc_deltaT_mod

    use MeshReaderCBinding
#ifdef HDF
    USE HDF5
#endif
    use jacobiNormal_mod
#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc, c_null_char
    use f_ftoc_bind_interoperability
#endif
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                                              !
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    REAL                           :: time,x,y,Variable(8),k1,k2               !
    INTEGER                        :: timestep                                 !
    REAL,POINTER                   :: pvar(:,:)                                !
    REAL,POINTER                   :: cvar(:,:)                                !
    TYPE (tEquations)              :: EQN                                      !
    TYPE (tInitialCondition)       :: IC                                       !
    TYPE (tUnstructMesh)           :: MESH                                     !
    TYPE (tMPI), OPTIONAL          :: MPI                                      !
    TYPE (tSource)                 :: SOURCE                                   !
    TYPE (tDiscretization)         :: DISC                                     !
    TYPE (tUnstructOptionalFields) :: OptionalFields                           !
    TYPE (tInputOutput)            :: IO                                       !
    TYPE (tBoundary)               :: BND                                      !
    CHARACTER(LEN=100)             :: programTitle                             !
    ! local variable declaration                                               !
    INTEGER                        :: iElem, iSide, iVtx, i                    !
    INTEGER                        :: counter
    INTEGER                        :: allocStat,iErr                           !
    INTEGER                        :: dummy(4)                                 !
    INTEGER                        :: nGraphVertex                             !
    INTEGER                        :: LocPoly                                  !
	INTEGER						   :: IntegrationMask(1:9)					   !
    INTEGER, POINTER               :: MetisWeight(:)                           !
    REAL                           :: AnelasticFactor(7,10)                    ! Factor of additional cost due to anelasticity
    REAL                           :: OrderFactor(7)                           ! Factor of additional cost due to order
    CHARACTER(LEN=256)             :: outfile
    CHARACTER(LEN=5)               :: cmyrank
    CHARACTER(LEN=600)             :: Filename, InFileName, OutFileName
    LOGICAL                        :: exists
    REAL, ALLOCATABLE              :: MaterialVals(:)
    REAL, ALLOCATABLE              :: WaveSpeeds(:)
#ifdef HDF
! HDF5 variables
     CHARACTER(LEN=100) :: file_name = "" ! File name
     CHARACTER(LEN=8) :: dsetname   ! Dataset name

     INTEGER(HID_T) :: file_id       ! File identifier
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
     INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
     INTEGER(HID_T) :: plist_id      ! Property list identifier

     INTEGER(HSIZE_T), DIMENSION(2) :: dimsf ! Dataset dimensions.
     INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi
     INTEGER        :: info   ! Number of variables to output
     INTEGER(HSIZE_T), DIMENSION(2) :: count
     INTEGER(HSSIZE_T), DIMENSION(2) :: offset
     REAL (4), ALLOCATABLE :: data (:,:)  ! Data to write
     INTEGER :: rank = 2 ! Dataset rank
#endif
     INTEGER :: error, error_n  ! Error flags
     INTEGER  :: t1, t2, clock_max      ! needed for dwalltime
     REAL     :: clock_rate = 0.0d0         ! needed for dwalltime
! HDF5 variables end
     integer :: l_ruptureFace, l_side, l_elementId, l_localFaceId ! loop counter of the dynamic rupture faces
     real :: l_jInv
     real, dimension(9, 9) :: l_A, l_B
     real, dimension(3, 3) :: l_attenuation

    !--------------------------------------------------------------------------
    INTENT(IN)                     :: programTitle                             !
    INTENT(INOUT)                  :: EQN,DISC,IO                              ! Some values are set in the TypesDef
    INTENT(OUT)                    :: IC,MESH,SOURCE,BND, OptionalFields       !
    INTENT(OUT)                    :: time,timestep             !

    ! register ini_seissol function
    EPIK_FUNC_REG("ini_SeisSol")
    SCOREP_USER_FUNC_DEFINE()
    ! register read_mesh/compute mesh region
    EPIK_USER_REG(r_read_compute_mesh, "read_compute_mesh")
    SCOREP_USER_REGION_DEFINE( r_read_compute_mesh )
    !--------------------------------------------------------------------------

    ! start epik/scorep function ini_SeisSol
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN( "ini_SeisSol" )

    time                     = 0.0                                             ! Initialize
    timestep                 = 0                                               ! Global iteration number, incl. restarts
    DISC%iterationstep       = 0                                               ! local iteration number (0 after restart)
    DISC%LoopCPUTime         = 0.0                                             ! Initialize CPU-time-measurement
    EQN%pi                   = ACOS(-1.)                                       ! Initialize
    MESH%MaxElementsOnVertex = 20                                              ! Initialize
    MESH%nZones              = 1                                               ! Default values
    !                                                                          !
    IO%ContourFile = 'contour'                                                 ! File for plot body contour
    !
    !
    CALL readpar(                                         &                    ! read the parameter file
         EQN          = EQN                             , &                    !
         IC           = IC                              , &                    !
         usMESH       = MESH                            , &                    !
         DISC         = DISC                            , &                    !
         SOURCE       = SOURCE                          , &                    !
         BND          = BND                             , &                    !
         IO           = IO                              , &                    !
         programTitle = programTitle                    , &                    !
         MPI          = MPI                               )

    !
    IF(MPI%nCPU.GT.1) THEN
     DISC%IterationCriterion = 3 ! Associated to local time step
    ENDIF
    !

#ifdef GENERATEDKERNELS
    ! Set ouput and checkpoint parameters
    ! This has to be done before the LTS setup!
    if( io%format .eq. 6 ) then
      call c_interoperability_enableWaveFieldOutput( i_waveFieldInterval = io%outInterval%timeInterval, &
                                                     i_waveFieldFilename = trim(io%OutputFile) // c_null_char )
    endif

    if( io%checkpoint%interval .gt. 0 ) then
        call c_interoperability_enableCheckPointing( i_checkPointInterval = io%checkpoint%interval, &
                                                     i_checkPointFilename = trim(io%checkpoint%filename) // c_null_char, &
                                                     i_checkPointBackend = trim(io%checkpoint%backend) // c_null_char )
    endif
    if( io%printEnergiesinterval .gt. 0 ) then
        call c_interoperability_setPrintEnergiesInterval( i_printEnergiesinterval = io%printEnergiesinterval) 
    endif

#ifdef INTEGRATE_QUANTITIES
	do i = 1,9
		if ( io%IntegrationMask(i) ) then
			IntegrationMask(i) = 1
		else
			IntegrationMask(i) = 0
		end if
	end do

	call c_interoperability_getIntegrationMask( i_integrationMask = IntegrationMask(1:9) )
#endif
#endif

    ! Start mesh reading/computing section
    EPIK_USER_START(r_read_compute_mesh)
    SCOREP_USER_REGION_BEGIN( r_read_compute_mesh, "read_compute_mesh", SCOREP_USER_REGION_TYPE_COMMON )
    call read_mesh_fast(IO,EQN,DISC,MESH,BND,MPI)
    !                                                                          !
    ! output neighbour list
    !OPEN(UNIT=999,FILE='Neighbourhood.dat')
    !   DO iElem = 1,MESH%nElem
    !         WRITE(999,*) MESH%ELEM%SideNeighbor(:,iElem)
    !   ENDDO
    !CLOSE(999)

    logInfo(*) '<--------------------------------------------------------->' !
    logInfo(*) '<           Calling DG Initialization level 1             >' !
    logInfo(*) '<--------------------------------------------------------->' !

    CALL iniGalerkin3D_us_level1_new(                 &              !
        EQN    = EQN                                , &              !
        DISC   = DISC                               , &              !
        MESH   = MESH                               , &              !
        BND    = BND                                , &              !
        IC     = IC                                 , &              !
        SOURCE = SOURCE                             , &              !
        OptionalFields = OptionalFields             , &              !
        MPI    = MPI                                , &              !
        IO     = IO                                   )              !

    ! End mesh reading/computing section
    EPIK_USER_END(r_read_compute_mesh)
    SCOREP_USER_REGION_END( r_read_compute_mesh )

#ifdef GENERATEDKERNELS
    logInfo(*) 'Generated Kernels: Checking boundary conditions'

    do iElem = 1, mesh%nElem
      do iSide = 1,4
        ! write an error message in the case of not supported boundary conditions
        if(       ( mesh%elem%reference( iSide, iElem ) .ne. 0 )& ! no boundary conditions: insides the computational domain
            .and. ( mesh%elem%reference( iSide, iElem ) .ne. 6 )& ! periodic boundary conditions
            .and. ( mesh%elem%reference( iSide, iElem ) .ne. 5 )& ! absorbing boundary conditions
            .and. ( mesh%elem%reference( iSide, iElem ) .ne. 1 )& ! free surface boundary conditions
            .and. ( mesh%elem%reference( iSide, iElem ) .ne. 3 )& ! dynamic rupture boundary conditions
          ) then
          logError(*) 'boundary conditions not supported'
          stop
        endif
      enddo
    enddo
#endif

    !                                                                               !
    DISC%NGP = 1                                                               !
    !                                                                          !
    !                                                                          !
    ALLOCATE(                                          &                       !
         pvar( MESH%nElem , EQN%nVar),                 &                       ! primitive variables
         cvar( MESH%nElem , EQN%nVar),                 &                       ! conservative var.
         STAT=allocStat                                )                       !
    !                                                                          !
    IF (allocStat .NE. 0) THEN                                                 ! Error Handling
       logError(*) 'could not allocate all variables!'      !
       STOP                                                                    ! STOP
    END IF                                                                     !
    !                                                                          !
    pvar(:,:)  = 0.                                                            ! Initialize
    cvar(:,:)  = 0.                                                            ! Initialize
    !                                                                          !
    !                                                                          ! Compute the number of background variables
    CALL ini_OptionalFields(                              &                    ! Allociert: OptionalFields%dt
         OptionalFields = OptionalFields                , &                    !            OptionalFields%RK_k
         SOURCE         = SOURCE                        , &                    !            OptionalFields%Residual
         EQN            = EQN                           , &                    !            OptionalFields%Backgroundvalue
         MESH           = MESH                          , &                    !
         DISC           = DISC                          , &                    !
         IO             = IO                              )                    !
    !                                                                          !
    !
    !T. Ulrich 08.2015 Read 2D basis Function for Vr output
    CALL Read2dGF(DISC,IO)
    !
    IF (EQN%linearized) THEN                                                   !
       CALL ini_MODEL(                                      &                  ! Initialize Local Linearized calculation
            MaterialVal    = OptionalFields%BackgroundValue,&                  ! Initialize Local Linearized calculation
            EQN            = EQN                          , &                  ! Initialize Local Linearized calculation
            MESH           = MESH                         , &                  ! Initialize Local Linearized calculation
            IO             = IO                           , &                  ! Initialize Local Linearized calculation
            DISC           = DISC                         , &                  ! Initialize Local Linearized calculation
            BND            = BND                            )                  ! Initialize Local Linearized calculation
    END IF                                                                     !
    logInfo(*) '<--------------------------------------------------------->' !
    logInfo(*) '<           Calling DG Initialization level 2             >' !
    logInfo(*) '<--------------------------------------------------------->' !
    !                                                                       !
    CALL iniGalerkin3D_us_level2_new(                    &              !
         EQN    = EQN                                , &              !
         DISC   = DISC                               , &              !
         MESH   = MESH                               , &              !
         BND    = BND                                , &              !
         IC     = IC                                 , &              !
         SOURCE = SOURCE                             , &              !
         OptionalFields = OptionalFields             , &              !
         MPI    = MPI                                , &              !
         IO     = IO                                   )              !
    !                                                                       !
    logInfo(*) 'Galerkin module initialized correctly.'    !
    !
    !aheineck: Metisweighs are not used, @TODO we should delete this
#if 0
    IF(IO%METISWeights) THEN ! not yet done for hybrids
        ! Compute Metis weights
        logInfo(*) 'Computing Metis weights ...'

        IF(BND%periodic.NE.0)THEN
           logError(*) 'MetisWeights can only be computed for non-periodic boundary conditions!'
           STOP
        ENDIF

        ALLOCATE(MetisWeight(MESH%nElem))

        !OrderFactor(:)       = (/4.34, 5.36, 8.39, 16.69, 37.08, 83.69, 180.70/)
        OrderFactor(:)       = (/1.00, 1.07, 1.29, 1.60, 2.62, 5.07, 10.29/)

        AnelasticFactor(1,:) = (/1.38, 1.48, 1.62, 1.68, 1.82, 1.88, 1.95, 2.08, 2.18, 2.27/)
        AnelasticFactor(2,:) = (/1.38, 1.48, 1.62, 1.68, 1.82, 1.88, 1.95, 2.08, 2.18, 2.27/)
        AnelasticFactor(3,:) = (/1.49, 1.70, 1.91, 2.07, 2.28, 2.49, 2.66, 2.82, 2.99, 3.20/)
        AnelasticFactor(4,:) = (/1.77, 2.07, 2.43, 2.70, 2.99, 3.34, 3.69, 4.03, 4.47, 4.67/)
        AnelasticFactor(5,:) = (/1.92, 2.30, 2.83, 3.27, 3.84, 4.38, 4.83, 5.40, 5.77, 6.32/)
        AnelasticFactor(6,:) = (/1.93, 2.44, 3.00, 3.55, 4.13, 4.65, 5.06, 5.82, 6.44, 6.87/)
        AnelasticFactor(7,:) = (/1.93, 2.44, 3.00, 3.55, 4.13, 4.65, 5.06, 5.82, 6.44, 6.87/)

        DO iElem = 1,MESH%nElem
           LocPoly = INT(DISC%Galerkin%LocPoly(iElem))
!           IF(EQN%LocAnelastic(iElem).EQ.1)THEN
!                  MetisWeight(iElem) = CEILING( (OrderFactor(LocPoly+1)                            &
!                                              * (AnelasticFactor(LocPoly+1,EQN%nMechanisms))) *100 )
!           ELSE
                  MetisWeight(iElem) = CEILING( OrderFactor(LocPoly+1) *100 )
!           ENDIF

        ENDDO

        !MetisWeight(:) = CEILING( REAL(MetisWeight(:) / minval(MetisWeight(:))) )
        logInfo(*) 'Minimum Metis weight:  ', minval(MetisWeight(:))
        logInfo(*) 'Maximum Metis weight:  ', maxval(MetisWeight(:))

        ! Read Metis graph file and output weighted Metis graph file
        InFileName  = TRIM(IO%MetisFile) // '.dgraph'
        OutFileName = TRIM(IO%MetisFile) // '.dgraph.weighted'

        OPEN(UNIT=998,FILE=TRIM(InFileName))
        OPEN(UNIT=999,FILE=TRIM(OutFileName))
        READ(998,*)  dummy(1), dummy(2)
        WRITE(999,*) dummy(1), dummy(2), 10, 1
        nGraphVertex = dummy(1)
        IF(nGraphVertex.EQ.MESH%nElem)THEN
           DO iElem = 1,MESH%nElem
              dummy(:) = 0
              counter  = 0
              DO iSide = 1,MESH%GlobalElemType
                 IF(MESH%ELEM%SideNeighbor(iSide,iElem).NE.MESH%nElem+1)THEN
                     counter = counter+1
                 ENDIF
              ENDDO
              READ(998,*)  dummy(1:counter)
              WRITE(999,*) MetisWeight(iElem), dummy(1:counter)
           ENDDO
        ELSE
           PRINT *, 'Number of Vertices in Graph File not equal to number of element of the mesh!'
        ENDIF

        CLOSE(998)
        CLOSE(999)

        logInfo(*) 'Metis weights computed successfully!'
        logInfo(*) 'weighted graph file written to ',TRIM(IO%MetisFile) // '.dgraph.weighted'
        logInfo(*) 'Terminating Programm !'
        STOP
    ENDIF
#endif

    ! TODO do we need to call this function when we read a checkpoint?? (SR)
    CALL icGalerkin3D_us_new(EQN, DISC, MESH, IC, SOURCE, IO)
    ! Pre-compute basis functions at Gauss points of non-conforming boundaries
 !   CALL NonConformingGPEvaluation3D(DISC, MESH, IO, BND)
 !   CALL DGSponge(0., EQN, DISC, MESH, BND, IC, SOURCE, IO) ! not yet done for hybrids/unified version
    !

    logInfo(*) 'Galerkin3D module initialized correctly.'     !
    logInfo(*) ' '                                             !
    !                                                                          !
    logInfo(*) '<--------------------------------------------------------->'  !
    logInfo(*) '<     ini_SeisSol finished successfully                   >'  !
    logInfo(*) '<--------------------------------------------------------->'  !
    !                                                                          !

    ! end epik function ini_SeisSol
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()
  END SUBROUTINE ini_SeisSol                                                    !

END MODULE ini_SeisSol_mod
