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

MODULE COMMON_readpar_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  USE COMMON_operators_mod
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  interface
    subroutine getParameterFileC( i_maxlen, o_file ) bind( C, name='getParameterFile' )
        use iso_c_binding
        implicit none
        integer(kind=c_int), value                        :: i_maxlen
        character(kind=c_char), dimension(*), intent(out) :: o_file
    end subroutine
  end interface
  !----------------------------------------------------------------------------
  PUBLIC  :: readpar
  !----------------------------------------------------------------------------

  LOGICAL :: CalledFromStructCode !

CONTAINS

  subroutine getParameterFile(parafile)
    use iso_c_binding
    character(*), intent(out) :: parafile
    integer(c_int) :: i_maxlen
    character(c_char) :: o_file(len(parafile))
    integer :: i
    parafile=''; i_maxlen = len(parafile)
    call getParameterFileC(i_maxlen,o_file)
    do i =1,i_maxlen
       if(o_file(i)==c_null_char)return
       parafile(i:i) = o_file(i)
    end do
  end subroutine getParameterFile

  SUBROUTINE readpar(EQN,IC,usMESH,DISC,SOURCE,BND,IO, &
                     programTitle,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)               :: EQN
    TYPE (tInitialCondition)        :: IC
    TYPE (tUnstructMesh) , OPTIONAL :: usMESH
    TYPE (tDiscretization)          :: DISC
    TYPE (tSource)                  :: SOURCE
    TYPE (tBoundary)                :: BND
    TYPE (tInputOutput)             :: IO
    TYPE (tMPI)          , OPTIONAL :: MPI
    CHARACTER(LEN=100)              :: programTitle
    ! local variables
    INTEGER                         :: actual_version_of_readpar
    CHARACTER(LEN=600)              :: Name
    CHARACTER(LEN=801)              :: Name1
    logical :: existence
    !--------------------------------------------------------------------------
    INTENT(IN)                      :: programTitle
    INTENT(OUT)                     :: IC, BND, DISC, SOURCE
    INTENT(INOUT)                   :: EQN,IO, usMESH
    !--------------------------------------------------------------------------
    PARAMETER(actual_version_of_readpar = 20)
    !--------------------------------------------------------------------------
    !                                                                        !
    IO%Mesh_is_structured       = .FALSE.                                    ! PostProcessing in default position
    SOURCE%Type                 = 0                                          ! switch for source terms in deflt pos.
    !                                                                        !
    IF (PRESENT(usMESH)) THEN                                                !
       CalledFromStructCode = .FALSE.                                        !
       DISC%CalledFromStructCode = .FALSE.                                   !
    ELSE                                                                     !
       logError(*) 'No MESH-Type in Argumentlist of readpar!'
       STOP                                                                  !
    END IF                                                                   !
    !                                                                        !
    call getParameterFile(IO%ParameterFile)

    ! Check if the parameter file exists
    inquire(file=IO%ParameterFile,EXIST=existence)
    if(existence)then
    else
       logError(*) 'You did not specify a valid parameter-file'
       stop
    endif
    !
    logInfo0(*) '<  Parameters read from file: ', TRIM(IO%ParameterFile) ,'              >'
    logInfo0(*) '<                                                         >'
    !                                                                        !
    CALL OpenFile(UnitNr=IO%UNIT%FileIn, name=trim(IO%ParameterFile), create=.FALSE.)
    !                                                                        !
    CALL readpar_header(IO,IC,actual_version_of_readpar,programTitle) !
    !                                                                        !
    CALL readpar_equations(EQN,DISC,SOURCE,IC,IO)                  !
    !                                                                        !
    CALL readpar_ini_condition(EQN,IC,SOURCE,IO)                   !
    !                                                                        !
    CALL readpar_boundaries(EQN,BND,IC,DISC,IO,CalledFromStructCode)    !
    !                                                                        !                                                                        !
    CALL readpar_sourceterm(EQN,SOURCE,IO)                                   !
    !                                                                        !
    CALL readpar_spongelayer(DISC,EQN,SOURCE,IO)
    !                                                                        !
    CALL readpar_mesh(EQN,IC,usMESH,DISC,BND,SOURCE,IO)            ! Unstrukturiertes Gitter definieren
    !                                                                        !
    CALL readpar_discretisation(EQN,usMESH,DISC,SOURCE,IO)         !
    !                                                                        !
    CALL readpar_output(EQN,DISC,IO,CalledFromStructCode)          !
    !                                                                        !
    CALL readpar_abort(DISC,IO)                                    !
    !                                                                        !
    CLOSE(IO%UNIT%FileIn,status='keep')                                      !
    !                                                                        !
    CALL analyse_readpar(EQN,DISC,usMESH,IC,SOURCE,IO,MPI)                   ! Check parameterfile...
    !                                                                        ! and write restart.par
    RETURN                                                                   !
    !                                                                        !
  END SUBROUTINE readpar                                                     !

  SUBROUTINE RaiseErrorNml(FID, NMLname)
    INTEGER                    :: FID
    CHARACTER(LEN=*)         :: NMLname
    CHARACTER(LEN=600)         :: line
    INTENT(IN)                 :: FID, NMLname

    backspace(FID)
    read(FID,fmt='(A)') line
    logError(*) 'invalid line in namelist '//trim(NMLname)//': '//trim(line)
    stop
    RETURN

  END SUBROUTINE

  !============================================================================
  ! H E A D E R
  !============================================================================

  SUBROUTINE readpar_header(IO, IC, actual_version_of_readpar,programTitle)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tInputOutput)        :: IO
    TYPE (tInitialCondition)   :: IC
    INTEGER                    :: actual_version_of_readpar
    INTEGER                    :: WeightsFlag
    CHARACTER(LEN=100)         :: programTitle
    !------------------------------------------------------------------------
    INTENT(IN)                 :: actual_version_of_readpar, programTitle
    INTENT(INOUT)              :: IO
    !------------------------------------------------------------------------
    !
    logInfo(*) '<<<<<<<<<<<<<<<<<<< ',TRIM(programTitle),' file log >>>>>>>>>>>>>>>>>>>>>>'
    logInfo(*) '<                                                         >'
    logInfo(*) '<        Welcome to ',TRIM(programTitle),' file logging system           >'
    logInfo(*) '<                                                         >'
    logInfo(*) '<--------------------------------------------------------->'

    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<              ',TRIM(programTitle),' - V E R S I O N                    >'
    logInfo(*) '<--------------------------------------------------------->'

    ! There should be only one version left, which I call here temporarily 19
    !logInfo(*) 'Parameter file is for ',TRIM(programTitle),' VERSION: 19'

  END SUBROUTINE readpar_header

  !============================================================================
  ! E Q U A T I O N S
  !============================================================================

  SUBROUTINE readpar_equations(EQN,DISC,SOURCE,IC, IO)
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tDiscretization)     :: DISC
    TYPE (tSource)             :: SOURCE
    TYPE (tInitialCondition)   :: IC
    TYPE (tInputOutput)        :: IO
    ! localVariables
    INTEGER                    :: intDummy
    INTEGER                    :: readStat
    !------------------------------------------------------------------------
    INTENT(INOUT)              :: EQN, IC, IO, SOURCE
    !------------------------------------------------------------------------
    LOGICAL                    :: fileExists
    INTEGER                    :: Anisotropy, Anelasticity, Plasticity, pmethod, Adjoint
    REAL                       :: FreqCentral, FreqRatio, Tv
    CHARACTER(LEN=600)         :: MaterialFileName, AdjFileName
    NAMELIST                   /Equations/ Anisotropy, Plasticity, &
                                           Tv, pmethod, &
                                           Adjoint,  &
                                           MaterialFileName, FreqCentral, &
                                           FreqRatio, AdjFileName
    !------------------------------------------------------------------------
    !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<  E Q U A T I O N S                                      >'
    logInfo(*) '<--------------------------------------------------------->'
    !
    EQN%Dimension = 3
    logInfo(*) 'Solve 3-dimensional equations'
    EQN%nvar = 9
    logInfo(*) 'Number of Variables    is ',EQN%nVar
    EQN%EqType = 8                                                ! equationtype
    !                                                             ! (8=seismic wave equations)
    !
    logInfo(*) 'Solving evolution equation for seismic wave propagation. '
    !
    EQN%linearized = .TRUE.
    ! aheineck, @TODO these values are used, but not initialized > Begin
    EQN%Poroelasticity = 0
    EQN%nBackgroundVar = 0
    EQN%Advection = 0
    ! aheineck, @TODO these values are used, but not initialized < End

    ! Setting the default values
    Anisotropy          = 0
#if NUMBER_OF_RELAXATION_MECHANISMS != 0
    Anelasticity        = 1
#else
    Anelasticity        = 0
#endif
    Plasticity          = 0
    Tv                  = 0.03  !standard value from SCEC benchmarks
    pmethod             = 0 !high-order approach as default for plasticity
    Adjoint             = 0
    MaterialFileName    = ''
    !
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Equations)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Equations")
    ENDIF
    !
    SELECT CASE(Anisotropy)
    CASE(0)
      logInfo(*) 'Isotropic material is assumed. '
      EQN%Anisotropy = Anisotropy
      EQN%nBackgroundVar = 3+EQN%nBackgroundVar
      EQN%nNonZeroEV = 3
    CASE(1)
      logInfo(*) 'Full triclinic material is assumed. '
      EQN%Anisotropy = Anisotropy
      EQN%nBackgroundVar = 22+EQN%nBackgroundVar
      EQN%nNonZeroEV = 3
    CASE DEFAULT
      logError(*) 'Choose 0 or 1 as anisotropy assumption. '
      STOP
    END SELECT
    !

#if defined(GENERATEDKERNELS) && defined(USE_PLASTICITY)
    if (Plasticity .eq. 0) then
      logError(*) 'Plasticity is disabled, but this version was compiled with Plasticity.'
      stop
    endif
#endif

    SELECT CASE(Plasticity)
    CASE(0)
      logInfo0(*) 'No plasticity assumed. '
      EQN%Plasticity = Plasticity                                                     !
    CASE(1)
#if defined(GENERATEDKERNELS) && !defined(USE_PLASTICITY)
       logError(*) 'Plasticity is assumed, but this version was not compiled with Plasticity.'
       stop
#else
       logInfo0(*) '(Drucker-Prager) plasticity assumed .'

#if defined(USE_PLASTIC_IP)
       logInfo0(*) 'Integration Points approach used for plasticity.'
#elif defined(USE_PLASTIC_NB)
       logInfo0(*) 'Nodal Basis approach used for plasticity.'

#endif

#endif
        EQN%Plasticity = Plasticity
        !first constant, can be overwritten in ini_model
        EQN%Tv = Tv
        EQN%PlastMethod = pmethod
        SELECT CASE (EQN%PlastMethod) !two different methods for plasticity
        CASE(0)
                 logInfo0(*) 'Plastic relaxation Tv is set to: ', EQN%Tv
                 logInfo0(*) 'High-order points are used for plasticity. '
        CASE(2)
                 logInfo0(*) 'Plastic relaxation Tv is set to: ', EQN%Tv
                 logInfo0(*) 'Average of an element is used for plasticity. '
        CASE DEFAULT
                 logError(*) 'ERROR: choose 0 or 2 as plasticity method'
        stop
        END SELECT
    CASE DEFAULT
      logError(*) 'Choose 0 or 1 as plasticity assumption. '
      STOP
    END SELECT


    EQN%Anelasticity = Anelasticity
    SELECT CASE(Anelasticity)
    CASE(0)
      logInfo(*) 'No attenuation assumed. '
      EQN%nAneMaterialVar = 3
      EQN%nMechanisms    = 0
      EQN%nAneFuncperMech= 0
      EQN%nVarTotal = EQN%nVar
      EQN%nBackgroundVar = 3
    CASE(1)
       logInfo(*) 'Viscoelastic attenuation assumed ... '
       EQN%nAneMaterialVar = 5        ! rho, mu, lambda, Qp, Qs
       EQN%nMechanisms = NUMBER_OF_RELAXATION_MECHANISMS
       IF(EQN%Anisotropy.NE.2)   THEN
         EQN%nAneFuncperMech = 6                                                    !
         logInfo(*) '... using ', EQN%nAneFuncperMech,' anelastic functions per Mechanism.'                                                   !
       ENDIF
       EQN%nVarTotal = EQN%nVar + EQN%nAneFuncperMech * EQN%nMechanisms
       EQN%nBackgroundVar  = 3 + EQN%nMechanisms * 4
    CASE DEFAULT
      logError(*) 'Choose 0 or 1 as anelasticity assumption. '
      STOP
    END SELECT

    DISC%Galerkin%CKMethod = 0
    !
      SELECT CASE(Adjoint)
      CASE(0)
         logInfo(*) 'No adjoint wavefield generated. '
         EQN%Adjoint = Adjoint
      CASE(1)
         logInfo(*) 'Adjoint wavefield simultaneously generated. '
         EQN%Adjoint = Adjoint
      CASE DEFAULT
        logError(*) 'Choose 0, 1 as adjoint wavefield assumption. '
        STOP
      END SELECT

    IF(EQN%Adjoint.EQ.1) THEN
     call readadjoint(IO, DISC, SOURCE, AdjFileName)
    END IF
    !
    inquire(file=MaterialFileName , exist=fileExists)
    if (.NOT. fileExists) then
     logError(*) 'Material file "', trim(MaterialFileName), '" does not exist.'
     STOP
    endif
    !
    EQN%MaterialFileName = MaterialFileName
    EQN%FreqCentral = FreqCentral
    EQN%FreqRatio = FreqRatio
    !
    intDummy = 1                                                  ! coordinate type index
    !                                                             ! (1=cartesian)
    EQN%CartesianCoordSystem = int_to_logical(                  & !
           int        = intDummy                              , & !
           TrueValue  = 1                                     , & !
           FalseValue = 2                                     , & !
          Default    = .TRUE.                                   ) !
                                                                  !
    IF (EQN%CartesianCoordSystem) THEN                            ! Cartesian system
    logInfo(*) 'Cartesian coordinate',&          !
               ' system (X,Y,Z) chosen'          !
    ELSE                                                          !
        logError(*) 'wrong coordinate system chosen! '
        stop
    END IF                                                        !
    !                                                             !
  END SUBROUTINE readpar_equations

    !------------------------------------------------------------------------
     !Reading the Random Field Files
    !------------------------------------------------------------------------
  SUBROUTINE readrffiles(IO, number, RF_Files)
    IMPLICIT NONE
    TYPE (tInputOutput)                               :: IO
    INTENT(INOUT)                                     :: IO
    INTEGER                                           :: number
    INTEGER                                           :: readStat
    CHARACTER(600), DIMENSION(:), ALLOCATABLE         :: RF_Files
    NAMELIST                                         /RFFile/ RF_Files
    !------------------------------------------------------------------------
    ALLOCATE(RF_Files(number))
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = RFFile)      ! Write in namelistfile RF_File(1) = ... and in the next line RF_Files(2) = ...
                                            ! according to the number of Random Fields
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "RFFile")
    ENDIF
  END SUBROUTINE
    !------------------------------------------------------------------------
     !Adjoint set to yes
    !------------------------------------------------------------------------
  SUBROUTINE readadjoint(IO, DISC, SOURCE, AdjFileName)
    IMPLICIT NONE
    TYPE (tInputOutput)                               :: IO
    TYPE (tDiscretization)                            :: DISC
    TYPE (tSource)                                    :: SOURCE
    INTENT(INOUT)                                     :: IO, SOURCE
    INTEGER                                           :: i
    CHARACTER(600)                                    :: AdjFileName
    !------------------------------------------------------------------------
      ! Read-in adjoint inversion parameters
      logInfo(*) 'Beginning adjoint inversion initialization. '
      logInfo(*) 'Inversion parameters read from ', TRIM(AdjFileName)
      CALL OpenFile(                                       &
        UnitNr       = IO%UNIT%other01                , &
        Name         = ADJFileName                    , &
        create       = .FALSE.                          )
      !
      ! Inversion information is read now
      !
      DO i = 1,4
        READ(IO%UNIT%other01,*)               ! Read unimportant comments
      ENDDO

      !KERNEL TYPE CHOICE
      READ(IO%UNIT%other01,*)  DISC%Adjoint%KernelType
      READ(IO%UNIT%other01,*)  DISC%Adjoint%nIter
      READ(IO%UNIT%other01,*)  DISC%Adjoint%IterIni
      READ(IO%UNIT%other01,*)  DISC%Adjoint%IterFin

      DO i = 1,4
        READ(IO%UNIT%other01,*)               ! Read unimportant comments
      ENDDO

      !DATA AND SYNTHETIC PREFIXES
      READ(IO%UNIT%other01,*)       SOURCE%AdjSource%DataFormat
      READ(IO%UNIT%other01,'(a14)') IO%ObsFile
      READ(IO%UNIT%other01,*)       SOURCE%AdjSource%nSamples
      READ(IO%UNIT%other01,*)       SOURCE%AdjSource%Dt
      READ(IO%UNIT%other01,*)       SOURCE%AdjSource%nVar

      DO i = 1,4
        READ(IO%UNIT%other01,*)               ! Read unimportant comments
      ENDDO

      !TF MISFIT PARAMETERS
      READ(IO%UNIT%other01,*)  DISC%Adjoint%MisfitType
      READ(IO%UNIT%other01,*)  DISC%Adjoint%fmin
      READ(IO%UNIT%other01,*)  DISC%Adjoint%fmax
      READ(IO%UNIT%other01,*)  DISC%Adjoint%w0

      DO i = 1,4
        READ(IO%UNIT%other01,*)               ! Read unimportant comments
      ENDDO

      !DATA TAPERING PARAMETERS
      READ(IO%UNIT%other01,*)  SOURCE%AdjSource%TaperType
      READ(IO%UNIT%other01,*)  SOURCE%AdjSource%wind_size
      READ(IO%UNIT%other01,*)  SOURCE%AdjSource%var_chk
      READ(IO%UNIT%other01,*)  SOURCE%AdjSource%Tol

      DO i = 1,4
        READ(IO%UNIT%other01,*)               ! Read unimportant comments
      ENDDO

      !KERNEL SMOOTHENING PARAMETERS
      READ(IO%UNIT%other01,*)  DISC%Adjoint%SpFilterType
      READ(IO%UNIT%other01,*)  DISC%Adjoint%SmoothSize

      CLOSE(IO%UNIT%other01)

  END SUBROUTINE

  !
  !============================================================================
  ! INITIAL CONDITION
  !============================================================================

  SUBROUTINE readpar_ini_condition(EQN,IC,SOURCE,IO)
    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tInitialCondition)   :: IC
    TYPE (tSource)             :: SOURCE
    TYPE (tInputOutput)        :: IO
    ! localVariables
    INTEGER                    :: i, j, k, iLambda, intDummy
    INTEGER                    :: counter, iZones, iVar
    INTEGER                    :: allocstat
    REAL                       :: lambda(3), Re, Im
    COMPLEX                    :: IU
    CHARACTER(LEN=600)         :: cdummy
    !------------------------------------------------------------------------
    INTENT(OUT)                :: IC
    INTENT(IN)                 :: EQN
    INTENT(INOUT)              :: IO, SOURCE
    !------------------------------------------------------------------------
    CHARACTER(Len=600)         :: cICType, IniConditionFile
    REAL                       :: xc(3), amplitude, hwidth(3)
    INTEGER                    :: nZones, variable
    INTEGER                    :: readStat
    NAMELIST                   /IniCondition/ cICType, variable, xc, amplitude, hwidth, &
                                              IniConditionFile, nZones
    !------------------------------------------------------------------------
    !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<  INITIAL CONDITION                                      >'
    logInfo(*) '<--------------------------------------------------------->'

    SOURCE%Type = 0         ! set dummy value, sources are specified later in readpar_sourceterm
                                          ! <------>
    ! Setting the default values = no source acting since amplitude is zero
    cICType = 'Zero'
    variable = 1
    xc(:) = 0.0                 ! x,y,z - coordinate, in inputfile you can choose different values vor x,y,z
    amplitude = 0.0
    hwidth(:) = 5.0e3           ! in inputfile you can choose different values for x,y,z
    !
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = IniCondition)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "IniCondition")
    ENDIF

    ! Renaming all variables in the beginning
     IC%cICType = cICType

     logInfo(*) 'Type of INITIAL CONDITION required: ', TRIM(IC%cICType)
       !
   SELECT CASE(IC%cICType)
   !
   CASE('Zero')
       logInfo(*) 'Zero initial condition'
    CASE('Planarwave')                                                                ! CASE tPlanarwave
       logInfo(*) 'Planarwave initial condition'
    CASE DEFAULT                                                             ! CASE DEFAULT
       logError(*) 'none of the possible'           ,&
            ' initial conditions was chosen'
       logError(*) TRIM(IC%cICType),'|'
       STOP
    END SELECT
    !
    logInfo(*) 'to calculate the initial values.'
    !
  END SUBROUTINE readpar_ini_condition

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------

  subroutine readpar_faultAtPickpoint(EQN,BND,IC,DISC,IO,CalledFromStructCode)
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tBoundary)           :: BND
    TYPE (tInitialCondition)   :: IC
    TYPE (tDiscretization)     :: DISC
    TYPE (tInputOutput)        :: IO
    LOGICAL                    :: CalledFromStructCode
    ! localVariables
    INTEGER                    :: allocStat, OutputMask(5), i
    INTEGER                    :: printtimeinterval
    INTEGER                    :: nOutPoints
    INTEGER                    :: readStat
    REAL, DIMENSION(:), ALLOCATABLE ::X, Y, Z
    CHARACTER(LEN=600)         :: PPFileName
    !------------------------------------------------------------------------
    INTENT(INOUT)              :: EQN, IO, DISC
    INTENT(INOUT)              :: BND
    !------------------------------------------------------------------------
    NAMELIST                   /Pickpoint/ printtimeinterval, OutputMask, nOutPoints, PPFileName
    !------------------------------------------------------------------------
    !
    !Setting default values
    printtimeinterval = 1
    OutputMask(1:3) = 1
    OutputMask(4:5) = 0
    !
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Pickpoint)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Pickpoint")
    ENDIF
    !
     DISC%DynRup%DynRup_out_atPickpoint%printtimeinterval = printtimeinterval   ! read time interval at which output will be written
     DISC%DynRup%DynRup_out_atPickpoint%OutputMask(1:5) =  OutputMask(1:5)      ! read info of desired output 1/ yes, 0/ no
                                                                                ! position: 1/ slip rate 2/ stress 3/ normal velocity
     DISC%DynRup%DynRup_out_atPickpoint%nOutPoints = nOutPoints                 ! 4/ in case of rate and state output friction and state variable
     logInfo(*) '| '
     logInfo(*) 'Record points for DR are allocated'
     logInfo(*) 'Output interval:',DISC%DynRup%DynRup_out_atPickpoint%printtimeinterval,'.'

     ALLOCATE(X(DISC%DynRup%DynRup_out_atPickpoint%nOutPoints))
     ALLOCATE(Y(DISC%DynRup%DynRup_out_atPickpoint%nOutPoints))
     ALLOCATE(Z(DISC%DynRup%DynRup_out_atPickpoint%nOutPoints))

      logInfo(*) ' Pickpoints read from ', TRIM(PPFileName)
      CALL OpenFile(                                 &
            UnitNr       = IO%UNIT%other01         , &
            Name         = PPFileName            , &
            create       = .FALSE.                          )
        DO i = 1, nOutPoints
          READ(IO%UNIT%other01,*) X(i), Y(i), Z(i)

            logInfo(*) 'Read in point :'
            logInfo(*) 'x = ', X(i)
            logInfo(*) 'y = ', Y(i)
            logInfo(*) 'z = ', Z(i)

       END DO
       CLOSE(IO%UNIT%other01)
      ALLOCATE ( DISC%DynRup%DynRup_out_atPickpoint%RecPoint(DISC%DynRup%DynRup_out_atPickpoint%nOutPoints),     &
                STAT = allocStat                             )

      IF (allocStat .NE. 0) THEN
            logError(*) 'could not allocate',&
                 ' all variables! Ie. Unstructured record Points'
            STOP
      END IF
      !
      DISC%DynRup%DynRup_out_atPickpoint%RecPoint(:)%X = X(:)
      DISC%DynRup%DynRup_out_atPickpoint%RecPoint(:)%Y = Y(:)
      DISC%DynRup%DynRup_out_atPickpoint%RecPoint(:)%Z = Z(:)

      logInfo(*) 'In total:',DISC%DynRup%DynRup_out_atPickpoint%nOutPoints,'.'

  end SUBROUTINE readpar_faultAtPickpoint
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------

  subroutine readpar_faultElementwise(EQN,BND,IC,DISC,IO,CalledFromStructCode)
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tBoundary)           :: BND
    TYPE (tInitialCondition)   :: IC
    TYPE (tDiscretization)     :: DISC
    TYPE (tInputOutput)        :: IO
    LOGICAL                    :: CalledFromStructCode
    ! localVariables
    INTEGER                    :: OutputMask(11)
    INTEGER                    :: printtimeinterval
    INTEGER                    :: printIntervalCriterion
    INTEGER                    :: refinement_strategy, refinement
    INTEGER                    :: readStat
    REAL                       :: printtimeinterval_sec
    !-----------------------------------------------------------------------
    INTENT(INOUT)              :: EQN, IO, DISC
    INTENT(INOUT)              :: BND
    NAMELIST                   /Elementwise/ printtimeinterval, OutputMask, refinement_strategy, &
                                                refinement, printIntervalCriterion,printtimeinterval_sec
    !Setting default values
    printtimeinterval = 2
    printtimeinterval_sec = 1d0
    printIntervalCriterion = 1
    OutputMask(:) = 1
    OutputMask(4:11) = 0
    refinement_strategy = 2
    refinement = 2
    !
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Elementwise)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Elementwise")
    ENDIF
    !
    DISC%DynRup%DynRup_out_elementwise%printIntervalCriterion = printIntervalCriterion
    if (printIntervalCriterion.EQ.1) THEN
        DISC%DynRup%DynRup_out_elementwise%printtimeinterval = printtimeinterval   ! read time interval at which output will be written
#ifdef GENERATEDKERNELS
        logError(*) 'The generated kernels version does no longer support printIntervalCriterion = 1 for elementwise fault output'
        stop
#endif
    else
        DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec = printtimeinterval_sec   ! read time interval at which output will be written
    endif

    ! if 2, printtimeinterval is set afterwards, when dt is known
    DISC%DynRup%DynRup_out_elementwise%OutputMask(1:11) =  OutputMask(1:11)      ! read info of desired output 1/ yes, 0/ no
                                                                                     ! position: 1/ slip rate 2/ stress 3/ normal velocity
                                                                                     ! 4/ in case of rate and state output friction and state variable
                                                                                     ! 5/ background values 6/Slip 7/rupture speed 8/final slip 9/peak SR
                                                                                     ! 10/rupture arrival 11/dynamic shear stress arrival


    DISC%DynRup%DynRup_out_elementwise%refinement_strategy = refinement_strategy

    IF (DISC%DynRup%DynRup_out_elementwise%refinement_strategy.NE.2 .AND. &
        DISC%DynRup%DynRup_out_elementwise%refinement_strategy.NE.1 .AND. &
        DISC%DynRup%DynRup_out_elementwise%refinement_strategy.NE.0) THEN
        logError(*) 'Undefined refinement strategy for fault output!'
        STOP
    ENDIF

    DISC%DynRup%DynRup_out_elementwise%refinement = refinement                 ! read info of desired refinement level : default 0

    !Dynamic shear stress arrival output currently only for linear slip weakening friction laws
    IF (OutputMask(11).EQ.1) THEN
        SELECT CASE (EQN%FL)
               CASE(2,6,13,16,17,29,30,103) !LSW friction law cases
                    !use only if RF_output=1
                    IF (OutputMask(10).EQ.1) THEN
                        ! set 'collecting DS time' to 1
                        DISC%DynRup%DS_output_on = 1
                    ELSE
                        DISC%DynRup%DynRup_out_elementwise%OutputMask(10) = 1
                        logInfo(*) 'RF output turned on when DS output is used'
                        ! set 'collecting DS time' to 1
                        DISC%DynRup%DS_output_on = 1
                    ENDIF
               CASE DEFAULT
                    logError(*) 'Dynamic shear stress arrival output only for LSW friction laws.'
        END SELECT
    ENDIF

    IF ((OutputMask(7).EQ.1) .AND. (DISC%DynRup%RFtime_on.EQ.0)) THEN
        ! set 'collecting RF time' to 1
        DISC%DynRup%RFtime_on = 1
    ENDIF


  end SUBROUTINE readpar_faultElementwise

  !============================================================================
  ! B O U N D A R I E S
  !============================================================================

  SUBROUTINE readpar_boundaries(EQN,BND,IC,DISC,IO,CalledFromStructCode)
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tBoundary)           :: BND
    TYPE (tInitialCondition)   :: IC
    TYPE (tDiscretization)     :: DISC
    TYPE (tInputOutput)        :: IO
    LOGICAL                    :: CalledFromStructCode
    ! localVariables
    INTEGER                    :: stat
    INTEGER                    :: allocStat
    INTEGER                    :: BC_fs, BC_nc, BC_dr, BC_if, BC_of, BC_pe
    INTEGER                    :: readStat
    !------------------------------------------------------------------------
    INTENT(INOUT)              :: EQN, IO, DISC
    INTENT(INOUT)              :: BND
    !------------------------------------------------------------------------
    NAMELIST                   /Boundaries/ BC_fs, BC_nc, BC_dr, BC_if, BC_of, BC_pe
    !------------------------------------------------------------------------


    logInfo(*) '<--------------------------------------------------------->'          !
    logInfo(*) '<  B O U N D A R I E S                                    >'          !
    logInfo(*) '<--------------------------------------------------------->'          !
    !

    ! Setting default values
    BC_fs = 0
    BC_nc = 0
    BC_dr = 0
    BC_if = 0
    BC_of = 0
    BC_pe = 0
    !
    READ (IO%UNIT%FileIn, IOSTAT=readStat, nml = Boundaries)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Boundaries")
    ENDIF
    !
      !
      BND%NoBndObjects(:) = 0
      !--------------------------------------------------------------------------------------------!
      ! Free surface boundaries
      !--------------------------------------------------------------------------------------------!
      !                                                                                            ! number of
      BND%NoBndObjects(1) = BC_fs                                                                      !free surface boundaries
      !                                                                                            !
      logInfo(*) ' '
      logInfo(*) 'The number of free      surfaces is',BC_fs
      logInfo(*) '-----------------------------------  '
      !                                                                                            !
      IF(BC_fs.NE.0)THEN                                                                              !
         ALLOCATE (BND%ObjFreeSurface(BC_fs), STAT = allocStat)                                       !
          !                                                                                        !
          IF (allocStat .NE. 0) THEN                                                               ! Error Handler
             logError(*) 'could not allocate surface wall variables!'                              ! Error Handler
             STOP                                                                                  ! Error Handler
          END IF                                                                                   ! Error Handler
      ENDIF
      !----------------------------------------------------------------------------------------!
      ! (Internal) boundaries at which non-conforming mesh is allowed
      !----------------------------------------------------------------------------------------!

                                                                                                   ! number of
                                                                                                   ! free surface boundaries
       BND%NoBndObjects(2) = BC_nc
       !                                                                                        !
       logInfo(*) ' '                                                              !
       logInfo(*) 'The number of non-conforming     is',BC_nc                    !
       logInfo(*) '-----------------------------------  '                     !
      !
      EQN%DR = 0 !By default no dynamic rupture
      DISC%DynRup%OutputPointType = 0 !By default no Output
      !
      !----------------------------------------------------------------------------------------!
      ! (Internal) boundaries at which dynamic rupture is allowed
      !----------------------------------------------------------------------------------------!
      !                                                                                        ! number of
      BND%NoBndObjects(3) = BC_dr                                                                 ! rupture inner boundaries
      !                                                                                        !
       IF(BC_dr.NE.0) EQN%DR = 1
          logInfo(*) ' '
          logInfo(*) 'The number of   rupture surfaces is',BC_dr
          logInfo(*) '-----------------------------------  '

          IF(EQN%DR.EQ.1) THEN
         call readdr(IO, EQN, DISC, BND, IC)
       ENDIF ! EQN%DR.EQ.1
      !--------------------------------------------------------------------------------------!
      ! Inflow
      !--------------------------------------------------------------------------------------!

      BND%NoBndObjects(4) = BC_if                                                               ! number of different Inflow boundary
      !
      logInfo(*) ' '
      logInfo(*) 'The number of inflow surfaces is',BC_if
      logInfo(*) '-----------------------------------  '
      !
      IF(BC_if.NE.0)THEN
         ALLOCATE (BND%ObjInflow(BC_if), STAT = allocStat)
         !
         IF (allocStat .NE. 0) THEN                                                          ! Error Handler
            logError(*) 'could not allocate Inflow variables!'                               ! Error Handler
            STOP                                                                             ! Error Handler
         END IF                                                                              ! Error Handler
      call readinfl(BND, IC, EQN, IO, DISC, BC_if)
      END IF
      !
      !--------------------------------------------------------------------------------------!
      ! Outflow
      !--------------------------------------------------------------------------------------!

      BND%NoBndObjects(5) = BC_of                                                               ! number of Outflow boundary
      !                                                                                      !
      logInfo(*) ' '                                                            !
      logInfo(*) 'The number of outflow   surfaces is',BC_of                  !
      logInfo(*) '-----------------------------------  '                   !
      !                                                                                      !
      IF(BC_of.NE.0)THEN                                                                        !
          ALLOCATE (BND%ObjOutflow(BC_of), STAT = allocStat)                                !
         !                                                                                   !
         IF (allocStat .NE. 0) THEN                                                          ! Error Handler
            logError(*) 'could not allocate Outflow wall variables!'      ! Error Handler
            STOP                                                                             ! Error Handler
         END IF                                                                              ! Error Handler
      END IF                                                                                 !
      !                                                                                      !
      !--------------------------------------------------------------------------------------!
      ! Periodic
      !--------------------------------------------------------------------------------------!

      BND%NoBndObjects(6) = BC_pe                                                               ! number of different Periodic boundary
      !                                                                                      !
      logInfo(*) ' '                                                            !
      logInfo(*) 'The number of connected surfaces is',BC_pe                  !
      logInfo(*) '-----------------------------------  '                   !


  END SUBROUTINE readpar_boundaries                                                          !

    !------------------------------------------------------------------------
     !Dynamic Rupture
    !------------------------------------------------------------------------
  SUBROUTINE readdr(IO, EQN, DISC, BND, IC)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    TYPE (tEquations)                      :: EQN
    TYPE (tDiscretization)                 :: DISC
    TYPE (tBoundary)                       :: BND
    TYPE (tInitialCondition)               :: IC
    INTENT(INOUT)                          :: IO, EQN, DISC, BND
    INTEGER                                :: FL, BackgroundType, Nucleation, inst_healing, RF_output_on, DS_output_on, &
                                              OutputPointType, magnitude_output_on,  energy_rate_output_on, read_fault_file,refPointMethod, SlipRateOutputType
    INTEGER                                :: readStat
    LOGICAL                                :: fileExists

    CHARACTER(600)                         :: FileName_BackgroundStress
    CHARACTER(LEN=600)                     :: ModelFileName
    REAL                                   :: RS_sv0, XRef, YRef, ZRef, GPwise,  &
                                              Mu_SNuc_ini, H_Length, RS_f0, &
                                              RS_sr0, RS_b, RS_iniSlipRate1, &
                                              RS_iniSlipRate2, v_star, L, t_0, Mu_W, &
                                              NucRS_sv0, r_s, energy_rate_printtimeinterval

    !------------------------------------------------------------------------
    NAMELIST                              /DynamicRupture/ FL, BackgroundType, &
                                                RS_sv0, XRef, YRef, ZRef,refPointMethod, FileName_BackgroundStress, &
                                                GPwise, inst_healing, &
                                                Mu_SNuc_ini, H_Length, RS_f0, &
                                                RS_sr0, RS_b, RS_iniSlipRate1, &
                                                RS_iniSlipRate2, v_star, L, t_0, Mu_W, &
                                                NucRS_sv0, r_s, RF_output_on, DS_output_on, &
                                                OutputPointType, magnitude_output_on, energy_rate_output_on, energy_rate_printtimeinterval,  &
                                                SlipRateOutputType, ModelFileName
    !------------------------------------------------------------------------

    ! Setting default values
    BackgroundType = 0
    FL = 0
    RF_output_on = 0
    DS_output_on = 0
    magnitude_output_on = 0
    energy_rate_output_on = 0
    energy_rate_printtimeinterval = 1
    OutputPointType = 3
    SlipRateOutputType = 0
    RS_sv0 = 0
    XRef = 0
    YRef = 0
    ZRef = 0
    refPointMethod=0
    GPwise = 1 !1=GPwise and 0=elementwise
    inst_healing = 0
    Mu_SNuc_ini=1.0
    H_Length = 0
    RS_f0 = 0
    RS_sr0 = 0
    RS_b = 0
    RS_iniSlipRate1 = 0
    RS_iniSlipRate2 = 0
    v_star = 0
    t_0 = 0
    L = 0
    Mu_W = 0
    NucRS_sv0 = 0
    r_s = 0
    ModelFileName = ''

    !FileName_BackgroundStress = 'tpv16_input_file.txt'

           ! Read-in dynamic rupture parameters
           READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = DynamicRupture)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "DynamicRupture")
    ENDIF
           logInfo(*) 'Beginning dynamic rupture initialization. '
           
    inquire(file=ModelFileName , exist=fileExists)
    if (.NOT. fileExists) then
     logError(*) 'Dynamic rupture model file "', trim(ModelFileName), '" does not exist.'
     STOP
    endif
    !
    DISC%DynRup%ModelFileName = ModelFileName

           !FRICTION LAW CHOICE
           EQN%FL = FL
           EQN%GPwise = GPwise
           EQN%refPointMethod = refPointMethod
           EQN%XRef = XRef
           EQN%YRef = YRef
           EQN%ZRef = ZRef

           IF (EQN%GPwise .EQ.1) THEN
               logInfo0(*) 'GPwise initialization. '
           ELSE
               logInfo0(*) 'elementwise initialization. '
           ENDIF

           !BACKGROUND VALUES
           DISC%DynRup%BackgroundType = BackgroundType
           SELECT CASE(DISC%DynRup%BackgroundType)
           CASE(0)
             EQN%RS_sv0 = RS_sv0
           CASE DEFAULT
             logError(*) 'Unknown Stress Background Type: ',DISC%DynRup%BackgroundType
             STOP
           END SELECT

           !FRICTION SETTINGS
           SELECT CASE(EQN%FL)
           CASE(0)
             CONTINUE
           CASE(2)
             DISC%DynRup%inst_healing = inst_healing ! instantaneous healing switch (1: on, 0: off)
           CASE(6) ! bimaterial with LSW
             DISC%DynRup%v_star = v_star
             DISC%DynRup%L = L
           CASE(16) ! SCEC TPV 16/17
             ! all parameters are defined in input file of INITIAL VALUES
             DISC%DynRup%inst_healing = 0
             DISC%DynRup%t_0      = t_0       ! forced rupture decay time
             CONTINUE
           CASE(3,4,7,101,103)
             DISC%DynRup%RS_f0 = RS_f0    ! mu_0, reference friction coefficient
             DISC%DynRup%RS_sr0 = RS_sr0  ! V0, reference velocity scale
             DISC%DynRup%RS_b = RS_b    ! b, evolution effect
             DISC%DynRup%Mu_W = Mu_W    ! mu_w, weakening friction coefficient
             DISC%DynRup%RS_iniSlipRate1 = RS_iniSlipRate1! V_ini1, initial sliding velocity
             DISC%DynRup%RS_iniSlipRate2 = RS_iniSlipRate2! V_ini2, initial sliding velocity
             DISC%DynRup%t_0      = t_0       ! forced rupture decay time
           CASE DEFAULT
             logError(*) 'Unknown friction law ',EQN%FL
             STOP
           END SELECT

           !OUTPUT
           ! rupture front (RF) output in extra files: on = 1, off = 0
           DISC%DynRup%RF_output_on = RF_output_on

           IF (DISC%DynRup%RF_output_on.EQ.1) THEN
              ! set 'collecting RF time' to 1
              DISC%DynRup%RFtime_on = 1
              logInfo0(*) 'RF output in extra files on'
           ELSE
              DISC%DynRup%RFtime_on = 0
           ENDIF
           !
           ! dynamic stress output on = 1, off = 0
           DISC%DynRup%DS_output_on = DS_output_on
           IF ((DISC%DynRup%DS_output_on.EQ.1).AND. (DISC%DynRup%RF_output_on.EQ.1)) THEN
               logInfo0(*) 'DS output on'
           ELSE IF ((DISC%DynRup%DS_output_on.EQ.1).AND. (DISC%DynRup%RF_output_on.EQ.0)) THEN
               logInfo0(*) 'DS output on. For ouput in files, RF_output is turned on.'
               DISC%DynRup%RF_output_on = 1
               DISC%DynRup%RFtime_on = 1
           ENDIF


           ! magnitude output on = 1, off = 0
           DISC%DynRup%magnitude_output_on = magnitude_output_on

           ! moment rate and frictional energy rate output on=1, off=0
           DISC%DynRup%energy_rate_output_on = energy_rate_output_on
           DISC%DynRup%energy_rate_printtimeinterval = energy_rate_printtimeinterval

           !
           DISC%DynRup%OutputPointType = OutputPointType
           DISC%DynRup%SlipRateOutputType = SlipRateOutputType

           !
           if (DISC%DynRup%OutputPointType .eq. 0) then
                logInfo0(*) 'Disabling fault output'
           elseif(DISC%DynRup%OutputPointType.EQ.3) THEN
                ! in case of OutputPointType 3, read in receiver locations:
                ! DISC%DynRup%DynRup_out_atPickpoint%nOutPoints is for option 3 the number of pickpoints
                call readpar_faultAtPickpoint(EQN,BND,IC,DISC,IO,CalledFromStructCode)
           ELSEIF(DISC%DynRup%OutputPointType.EQ.4) THEN
                ! elementwise output -> 2 dimensional fault output
                call readpar_faultElementwise(EQN,BND,IC,DISC,IO,CalledFromStructCode)
           ELSEIF(DISC%DynRup%OutputPointType.EQ.5) THEN
                ! ALICE: TO BE DONE
                ! fault receiver + 2 dimensional fault output
                call readpar_faultElementwise(EQN,BND,IC,DISC,IO,CalledFromStructCode)
                call readpar_faultAtPickpoint(EQN,BND,IC,DISC,IO,CalledFromStructCode)
           ELSE
               logError(*) 'Unkown fault output type (e.g.3,4,5)',DISC%DynRup%OutputPointType
               STOP
           ENDIF ! DISC%DynRup%OutputPointType
  !
  END SUBROUTINE
    !------------------------------------------------------------------------
     !Inflow Boundaries
    !------------------------------------------------------------------------
  SUBROUTINE readinfl(BND, IC, EQN, IO, DISC, n4)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    TYPE (tEquations)                      :: EQN
    TYPE (tDiscretization)                 :: DISC
    TYPE (tBoundary)                       :: BND
    TYPE (tInitialCondition)               :: IC
    INTENT(INOUT)                          :: IO, EQN, DISC, BND
    INTEGER                                :: n4, i, j, iVar, readstat, startComment
    INTEGER                                :: nsteps, nPW, iObject, isteps
    REAL                                   :: nx,ny,nz,length
    CHARACTER(LEN=600)                     :: char_option, PWFileName
    INTEGER                                :: setvar
    REAL,DIMENSION(:),ALLOCATABLE          :: varfield, u0_in
    NAMELIST                               /InflowBound/ setvar, char_option, &
                                                         PWFileName
    !------------------------------------------------------------------------
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = InflowBound)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "InflowBound")
    ENDIF

      DO i=1,n4
         j = 1
         DO WHILE(char_option(j:j).NE.'!'.AND.j.LE.LEN(char_option)-1)
            j=j+1
         END DO
         IF (char_option(j:j).EQ.'!') THEN
            startComment= j
         ELSE
         startComment=LEN(char_option)
         END IF

              ALLOCATE(u0_in(EQN%nVar))
              call readuin(IO, EQN, u0_in)
              ALLOCATE(BND%ObjInflow(i)%u0_in(EQN%nVar) )
              BND%ObjInflow(i)%u0_in(:) = u0_in(:)

         READ(char_option,*,IOSTAT=readStat) BND%ObjInflow(i)%u0_in(:)

         IF (readStat.EQ. 0) THEN
            logInfo(*) 'Inflow conditions specified are constant Data'
            BND%ObjInflow(i)%InflowType   = 0
            logInfo(*) 'Inflow values for object n = ', i                   !
            DO iVar = 1, EQN%nVar
                logInfo(*) 'Variable ', iVar, ' = ',  BND%ObjInflow(i)%u0_in(iVar)
            ENDDO
            logInfo(*) ' '
         ELSE
            logInfo(*) 'Inflow conditions specified are not constant Data'

            SELECT CASE(TRIM(char_option(1:startComment-1)))
            CASE('Char_Gauss_Puls')
                logInfo(*) 'Planarwave Gausspulse inflow conditions '
                IF(IC%cICType.NE.'Char_Gauss_Puls') THEN
                  logError(*) 'Inflow boundary condition Char_Gauss_Puls only available with '
                  logError(*) 'corresponding Char_Gauss_Puls initial condition '
                  stop
                ENDIF
                BND%ObjInflow(i)%InflowType   = 1

            CASE('Char_Ricker_Puls')
                logInfo(*) 'Planarwave Rickerpulse inflow conditions '
                IF(IC%cICType.NE.'Char_Ricker_Puls') THEN
                  logError(*) 'Inflow boundary condition Char_Ricker_Puls only available with '
                  logError(*) 'corresponding Char_Ricker_Puls initial condition '
                ENDIF
                BND%ObjInflow(i)%InflowType   = 2

            CASE('Custom_PlaneWave_File')
                logInfo(*) 'Custom Planarwave from file inflow conditions. '
                logInfo(*) '(only for homogeneous elastic isotropic materials!) '
                BND%ObjInflow(i)%InflowType   = 10
                iObject=1 ! So far limited to a single PlaneWave!
                ALLOCATE(BND%ObjInflow(iObject)%PW)
                ALLOCATE(BND%ObjInflow(iObject)%u0_in(EQN%nVar) )
                ALLOCATE(BND%ObjInflow(iObject)%PW%n_vec(3))
                ALLOCATE(BND%ObjInflow(iObject)%PW%p_vec(3))
                ALLOCATE(BND%ObjInflow(iObject)%PW%t_vec(3))
                ALLOCATE(BND%ObjInflow(iObject)%PW%FI_Point(3))
                BND%ObjInflow(iObject)%PW%SetVar = SetVar
                ALLOCATE(BND%ObjInflow(iObject)%PW%varfield(BND%ObjInflow(iObject)%PW%SetVar))

                call readvarfield(IO, SetVar, varfield)

                BND%ObjInflow(iObject)%PW%varfield(:) = varfield (:)
                BND%ObjInflow(iObject)%u0_in(:) = u0_in(:)

                logInfo(*) 'Wave time histories read from file: ', TRIM(PWFileName)
                CALL OpenFile(                                        &
                      UnitNr       = IO%UNIT%other01                , &
                      Name         = PWFileName                     , &
                      create       = .FALSE.                          )
                logInfo(*) 'Reading inflow wave time-history file ...  '
                READ(IO%UNIT%other01,'(i10,a)') nsteps                        ! Number of timesteps included
                BND%ObjInflow(iObject)%PW%TimeHistnSteps = nsteps
                ALLOCATE(BND%ObjInflow(iObject)%PW%TimeHist(2,nsteps))
                ! Read normal vector

                READ(IO%UNIT%other01,*) BND%ObjInflow(iObject)%PW%n_vec(1:3)  ! Normal vector to wave propagation
                nx = BND%ObjInflow(iObject)%PW%n_vec(1)
                ny = BND%ObjInflow(iObject)%PW%n_vec(2)
                nz = BND%ObjInflow(iObject)%PW%n_vec(3)
                length = SQRT(nx*nx+ny*ny+nz*nz)
                BND%ObjInflow(iObject)%PW%n_vec(1) = nx/length
                BND%ObjInflow(iObject)%PW%n_vec(2) = ny/length
                BND%ObjInflow(iObject)%PW%n_vec(3) = nz/length
                ! Read tangential vector
                READ(IO%UNIT%other01,*) BND%ObjInflow(iObject)%PW%p_vec(1:3)  ! Polarization vector of wave
                nx = BND%ObjInflow(iObject)%PW%p_vec(1)
                ny = BND%ObjInflow(iObject)%PW%p_vec(2)
                nz = BND%ObjInflow(iObject)%PW%p_vec(3)
                length = SQRT(nx*nx+ny*ny+nz*nz)
                BND%ObjInflow(iObject)%PW%p_vec(1) = nx/length
                BND%ObjInflow(iObject)%PW%p_vec(2) = ny/length
                BND%ObjInflow(iObject)%PW%p_vec(3) = nz/length
                ! compute orthogonal vector to normal and tangential using the crossproduct
                BND%ObjInflow(iObject)%PW%t_vec(:)= BND%ObjInflow(iObject)%PW%n_vec(:) .x. &
                                                    BND%ObjInflow(iObject)%PW%p_vec(:)
                nx = BND%ObjInflow(iObject)%PW%t_vec(1)
                ny = BND%ObjInflow(iObject)%PW%t_vec(2)
                nz = BND%ObjInflow(iObject)%PW%t_vec(3)
                length = SQRT(nx*nx+ny*ny+nz*nz)
                BND%ObjInflow(iObject)%PW%t_vec(1) = nx/length
                BND%ObjInflow(iObject)%PW%t_vec(2) = ny/length
                BND%ObjInflow(iObject)%PW%t_vec(3) = nz/length

                READ(IO%UNIT%other01,*) BND%ObjInflow(iObject)%PW%FI_Point(1:3)  ! Point of first incidence of wave in the mesh
                DO isteps=1,nsteps
                  READ(IO%UNIT%other01,*) BND%ObjInflow(iObject)%PW%TimeHist(1:2,isteps)
                ENDDO
                CLOSE(IO%UNIT%other01)
                logInfo(*) 'Time histories successfully read from '
                logInfo(*) '  ',BND%ObjInflow(iObject)%PW%TimeHist(1,1),' seconds to '
                logInfo(*) '  ',BND%ObjInflow(iObject)%PW%TimeHist(1,nsteps),' seconds!'
           CASE DEFAULT
               logError(*) 'Inflow conditions specified are unknown!'
               logError(*) TRIM(char_option(1:startComment-1)),'|'
               STOP
            END SELECT
         END IF
      ENDDO

  END SUBROUTINE
    !------------------------------------------------------------------------
     !Reading the varfield
    !------------------------------------------------------------------------
 SUBROUTINE readvarfield(IO, number, varfield)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    INTENT(INOUT)                          :: IO
    INTEGER                                :: number
    INTEGER                                :: readStat
    REAL, DIMENSION(:), ALLOCATABLE        :: varfield
    NAMELIST                               /InflowBoundPWFile/ varfield
    !-----------------------------------------------------------------------
    ALLOCATE(varfield(number))
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = InflowBoundPWFile) ! Write in namelistfile varfield(1) = ... and in the next line varfield(2) = ...
                                                  ! and the same for u0_in
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "InflowBoundPWFile")
    ENDIF
  END SUBROUTINE
    !------------------------------------------------------------------------
     !Reading u0_in
    !------------------------------------------------------------------------
 SUBROUTINE readuin(IO, EQN, u0_in)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    TYPE (tEquations)                      :: EQN
    INTENT(INOUT)                          :: IO, EQN
    REAL, DIMENSION(:), ALLOCATABLE        :: u0_in
    INTEGER                                :: readStat
    NAMELIST                               /InflowBounduin/u0_in
    !-----------------------------------------------------------------------
    ALLOCATE(u0_in(EQN%nVar))

    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = InflowBounduin) ! Write in namelistfile u0_in(1) = ... and in the next line u0_in(2) = ...
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "InflowBounduin")
    ENDIF

  END SUBROUTINE


   !
   !============================================================================
   ! S O U R C E   T E R M
   !============================================================================

  SUBROUTINE readpar_sourceterm(EQN,SOURCE,IO)
   !------------------------------------------------------------------------
   IMPLICIT NONE
   !------------------------------------------------------------------------
   ! argument list declaration
   TYPE (tEquations)                :: EQN
   TYPE (tSource)                   :: SOURCE
   TYPE (tInputOutput)              :: IO
   ! local variable declaration
   INTEGER                          :: i,j,k,iVar,iDummy1
   INTEGER                          :: startComment
   REAL                             :: pi
   REAL, POINTER                    :: dummy(:)
   CHARACTER(LEN=15)                :: char_dummy
   REAL, POINTER                    :: Dist(:,:), loc_src(:,:)
   INTEGER,POINTER                  :: Sorted(:),Group(:),Members(:,:)
   INTEGER                          :: count, count2,nsrc
   !-------------------------------------------------------------------------
   INTENT(INOUT)                    :: SOURCE, EQN
   INTENT(IN)                       :: IO
   !-------------------------------------------------------------------------
   INTEGER                          :: Type, RType, nDirac, nRicker, nPulseSource
   REAL,DIMENSION(:),ALLOCATABLE    :: U0, l1, TimePosition, &
                                       Intensity, EqnNr, Element, Delay, a1, f, &
                                       l2, T, t0, Width, A0
   REAL,DIMENSION(:),ALLOCATABLE    :: SpacePositionx, SpacePositiony, SpacePositionz
   CHARACTER(Len=600)               :: FileName
   INTEGER                          :: readStat
   NAMELIST                        /SourceType/ Type, Rtype, nDirac, nPulseSource, FileName, nRicker

   !------------------------------------------------------------------------
   !
    logInfo(*) '<--------------------------------------------------------->'          !
    logInfo(*) '<  S O U R C E    T E R M S                               >'          !
    logInfo(*) '<--------------------------------------------------------->'          !
    !
    SOURCE%CauchyKovalewski = .FALSE.
    !
    ! Setting default values
    Type = 0
    !
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = SourceType)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "SourceType")
    ENDIF
    SOURCE%Type = Type
   SELECT CASE(SOURCE%Type)                                                 !

    CASE(0)
                                                                    ! No Source Term
       logInfo(*) 'No source specified.'                    !
       !
    CASE(1)
       logInfo(*) 'Source for convergence study of varying coefficient PDE specified.'
       ALLOCATE(SOURCE%CS%U0(3))
       ALLOCATE(SOURCE%CS%k1(3))
       ALLOCATE(SOURCE%CS%l1(3))
       call readsource110(IO, U0, l1)
       SOURCE%CS%U0(:) = U0(:)   ! perturbation amplitudes for A, B, and C
       SOURCE%CS%l1(:) = l1(:)   ! wavelengths
       SOURCE%CS%k1 = 2*EQN%Pi/SOURCE%CS%l1(:)
       !
    CASE(10)
       logInfo(*) 'Source for convergence study of varying coefficient PDE specified.'
       ALLOCATE(SOURCE%CS%k1(3) )
       ALLOCATE(SOURCE%CS%l1(3))
       call readsource110(IO, U0, l1)
       SOURCE%CS%l1(:) = l1(:)   ! wavelengths
       SOURCE%CS%k1 = 2*EQN%Pi/SOURCE%CS%l1(:)
       !
    CASE(15)
       logInfo(*) 'Dirac sources in space and time chosen. '
       SOURCE%Dirac%nDirac = nDirac
       logInfo(*) 'Number of Dirac sources: ', SOURCE%Dirac%nDirac

       ALLOCATE( SOURCE%Dirac%SpacePosition(3,SOURCE%Dirac%nDirac), &
                 SOURCE%Dirac%TimePosition(SOURCE%Dirac%nDirac), &
                 SOURCE%Dirac%Intensity(SOURCE%Dirac%nDirac), &
                 SOURCE%Dirac%EqnNr(SOURCE%Dirac%nDirac), &
                 SOURCE%Dirac%Element(SOURCE%Dirac%nDirac) )

       call readsource15(IO, nDirac, SpacePositionx, SpacePositiony, SpacePositionz, TimePosition, Intensity, EqnNr)
          ! copy on data structure:
          SOURCE%Dirac%SpacePosition(1,:) = SpacePositionx(:)
          SOURCE%Dirac%SpacePosition(2,:) = SpacePositiony(:)
          SOURCE%Dirac%SpacePosition(3,:) = SpacePositionz(:)
          SOURCE%Dirac%TimePosition(:) = TimePosition(:)
          SOURCE%Dirac%Intensity(:) = Intensity(:)
          SOURCE%Dirac%EqnNr(:) = EqnNr(:)


    CASE(16,18)

         SOURCE%CauchyKovalewski = .TRUE.
         SOURCE%Ricker%nRicker = nRicker
       IF(SOURCE%Type.EQ.16) THEN
         logInfo(*) 'Dirac sources in space and Ricker wavelet in time chosen. '
         logInfo(*) 'Number of Ricker sources: ', SOURCE%Ricker%nRicker
       ELSEIF(SOURCE%Type.EQ.18) THEN
         logInfo(*) 'Dirac sources in space and Gaussian wavelet in time chosen. '
         logInfo(*) 'Number of Gaussian sources: ', SOURCE%Ricker%nRicker
       ENDIF

       ALLOCATE( SOURCE%Ricker%SpacePosition(3,SOURCE%Ricker%nRicker), &
                 SOURCE%Ricker%Delay(SOURCE%Ricker%nRicker), &
                 SOURCE%Ricker%a1(SOURCE%Ricker%nRicker), &
                 SOURCE%Ricker%f(SOURCE%Ricker%nRicker), &
                 SOURCE%Ricker%EqnNr(SOURCE%Ricker%nRicker), &
                 SOURCE%Ricker%Element(SOURCE%Ricker%nRicker))

      call readsource1618(IO, nRicker, SpacePositionx, SpacePositiony, SpacePositionz, Delay, a1, f, EqnNr)

         ! copy on data structure:
         SOURCE%Ricker%SpacePosition(1,:) = SpacePositionx(:)
         SOURCE%Ricker%SpacePosition(2,:) = SpacePositiony(:)
         SOURCE%Ricker%SpacePosition(3,:) = SpacePositionz(:)
         SOURCE%Ricker%Delay(:) = Delay(:)
         SOURCE%Ricker%a1(:) = a1(:)
         SOURCE%Ricker%f(:)  = f(:)                    ! F takes the role of "sigma" for the gaussian case!
         SOURCE%Ricker%EqnNr(:) = EqnNr(:)

    CASE(17)

       logInfo(*) 'Sourceterm for convergence studies chosen. '
       logInfo(*) 'Sourceterm is a space-time periodic function consisting of sine-waves.  '
       !
       ALLOCATE(SOURCE%CS%U0(EQN%nVar))
       ALLOCATE(SOURCE%CS%l1(EQN%nVar))
       ALLOCATE(SOURCE%CS%l2(EQN%nVar))
       ALLOCATE(SOURCE%CS%k1(EQN%nVar))
       ALLOCATE(SOURCE%CS%k2(EQN%nVar))
       ALLOCATE(SOURCE%CS%k3(EQN%nVar))
       ALLOCATE(SOURCE%CS%T(EQN%nVar))
       ALLOCATE(SOURCE%CS%omega(EQN%nVar))
       !
         SOURCE%CauchyKovalewski = .TRUE.
       call readsource17(EQN, IO, U0, l1, l2, T)
          ! copy on data structure:
          SOURCE%CS%U0(:) = U0(:)
          SOURCE%CS%l1(:) = l1(:)         ! Read wavelengths (periods), compute wave numbers / frequencies
          SOURCE%CS%l2(:) = l2(:)
          !SOURCE%CS%k3(iVar)
          SOURCE%CS%T(:)  = T(:)
          SOURCE%CS%k1(:)    = 2.*EQN%Pi/SOURCE%CS%l1(:)
          SOURCE%CS%k2(:)    = 2.*EQN%Pi/SOURCE%CS%l2(:)
          SOURCE%CS%omega(:) = 2.*EQN%Pi/SOURCE%CS%T(:)
    CASE(19)

         SOURCE%CauchyKovalewski = .TRUE.
          ! Time Gauss pulse source

        logInfo(*) 'Time gauss pulse source chosen. '

        SOURCE%TimeGP%nPulseSource = nPulseSource
        ALLOCATE(SOURCE%TimeGP%SpacePosition(3,SOURCE%TimeGP%nPulseSource))
        ALLOCATE(SOURCE%TimeGP%t0(SOURCE%TimeGP%nPulseSource))
        ALLOCATE(SOURCE%TimeGP%Width(SOURCE%TimeGP%nPulseSource))
        ALLOCATE(SOURCE%TimeGP%A0(SOURCE%TimeGP%nPulseSource))
        ALLOCATE(SOURCE%TimeGP%Element(SOURCE%TimeGP%nPulseSource))

        call readsource19(IO, nPulseSource, EqnNr, SpacePositionx, SpacePositiony, SpacePositionz, t0, Width, A0)
            ! copy on data structure:
            SOURCE%TimeGP%SpacePosition(1,:) = SpacePositionx(:)
            SOURCE%TimeGP%SpacePosition(2,:) = SpacePositiony(:)
            SOURCE%TimeGP%SpacePosition(3,:) = SpacePositionz(:)
            SOURCE%TimeGP%t0(:) = t0(:)
            SOURCE%TimeGP%Width(:) = Width(:)
            SOURCE%TimeGP%A0(:) = A0(:)
            SOURCE%TimeGP%EqnNr(:) = EqnNr(:)

        DO i = 1, SOURCE%TimeGP%nPulseSource
            logInfo('("Source ",I4 )') i
            logInfo('("   x0  = ",E12.5,"   y0 = ",E12.5,"  z0 = ",E12.5)') SOURCE%TimeGP%SpacePosition(:,i)
            logInfo('("   t0  = ",E12.5)') SOURCE%TimeGP%t0(i)
            logInfo('("   tau = ",E12.5)') SOURCE%TimeGP%Width(i)
            logInfo('("   M0  = ",E12.5)') SOURCE%TimeGP%A0(i)
            logInfo('("   Affecetd var  = ",I4)') SOURCE%TimeGP%EqnNr(i)

        ENDDO

    CASE(20) !Single Force with individual slip rate history for each subfault

       logInfo(*) 'Single Force chosen. '
       SOURCE%FSRMFileName = FileName
       logInfo(*) 'Source term read from ', TRIM(SOURCE%FSRMFileName)
       CALL OpenFile(                                       &
            UnitNr       = IO%UNIT%other01                , &
            Name         = SOURCE%FSRMFileName            , &
            create       = .FALSE.                          )
       logInfo(*) 'Reading single force file ...  '
       !
       ! LEGEND
       !
       ! Number of Sources
       ! 3
       ! Single Force on Variable Nr.
       ! 7
       ! 8
       ! 9
       ! x y z
       ! 10. 0. 20.
       ! 10. 0. 20.
       ! 10. 0. 20.
       ! Source Time Functions
       ! 0.001 3                    ! all STF must have same sampling rate
       ! Samples
       ! 0.                         ! example of 3 delta functions acting on components 7,8,9
       ! 1.3
       ! 0.
       ! 0.
       ! 1.3
       ! 0.
       ! 0.
       ! 1.3
       ! 0.
       ! END LEGEND
       !
       ! Single Force information is read now
       !
       READ(IO%UNIT%other01,*)                                           ! Read comment
       READ(IO%UNIT%other01,*) SOURCE%Ricker%nRicker
       logInfo(*) 'Number of single forces: ', SOURCE%Ricker%nRicker
       ALLOCATE( SOURCE%Ricker%EqnNr(SOURCE%Ricker%nRicker) )
       ALLOCATE ( dummy(3) )
       READ(IO%UNIT%other01,*)                                           ! Read comment
       DO i = 1,SOURCE%Ricker%nRicker
           READ(IO%UNIT%other01,*) dummy(1)                             ! Read variable
           SOURCE%Ricker%EqnNr(i) = dummy(1)
       ENDDO
       READ(IO%UNIT%other01,*)                                           ! Read comment
       ALLOCATE ( SOURCE%RP%SpacePosition(3,SOURCE%Ricker%nRicker), &
                  SOURCE%RP%Element(SOURCE%Ricker%nRicker)    )
       DO i = 1,SOURCE%Ricker%nRicker
           READ(IO%UNIT%other01,*) dummy                                 ! Read coordinate data
           SOURCE%RP%SpacePosition(1:3,i) = dummy(1:3)
       ENDDO

       DEALLOCATE ( dummy )

       READ(IO%UNIT%other01,*)                                           ! Read comment
       READ(IO%UNIT%other01,*) SOURCE%RP%t_samp, SOURCE%RP%nsteps        ! Read time sampling of rupture functions

       SOURCE%RP%T_max = SOURCE%RP%t_samp * (SOURCE%RP%nsteps-1)         ! Equispaced time sampling for all subfaults assumed

       ALLOCATE ( SOURCE%RP%TimeHist(SOURCE%RP%nsteps, SOURCE%Ricker%nRicker)  )

       READ(IO%UNIT%other01,*)                                           ! Read comment
       DO i = 1,SOURCE%Ricker%nRicker
           DO j = 1,SOURCE%RP%nsteps
               READ(IO%UNIT%other01,*) SOURCE%RP%TimeHist(j,i)           !
           ENDDO
       ENDDO
       !
       CLOSE(IO%UNIT%other01)

    CASE(30)        ! Read in Finite Source Rupture Model
                    ! for input formats see: Martin Mai (http://www.seismo.ethz.ch/srcmod/Events.html)

       pi        = ACOS(-1.0)

       logInfo(*) 'Finite Source Rupture Model chosen. '
       SOURCE%FSRMFileName = FileName
       SOURCE%RP%Type = RType
       logInfo(*) 'Sourceterm read from ', TRIM(SOURCE%FSRMFileName)
       CALL OpenFile(                                       &
            UnitNr       = IO%UNIT%other01                , &
            Name         = SOURCE%FSRMFileName            , &
            create       = .FALSE.                          )
       logInfo(*) 'Reading rupture model file of Type  ',SOURCE%RP%Type
       !
       ! Rupture Plane Geometrie and Subfault information is read now
       !
       DO i = 1,5
          READ(IO%UNIT%other01,*)               ! Read unimportant comments
       ENDDO
       READ(IO%UNIT%other01,'(a15,f8.3,a11,f8.3,a11,f7.2)')    char_dummy, SOURCE%Hypocenter%latitude,   &
                                                               char_dummy, SOURCE%Hypocenter%longitude,  &
                                                               char_dummy, SOURCE%Hypocenter%depth
       READ(IO%UNIT%other01,'(a15,f8.3,a11,f8.3,a39)')         char_dummy, SOURCE%RP%length,             &
                                                               char_dummy, SOURCE%RP%width,              &
                                                               char_dummy
       READ(IO%UNIT%other01,'(a15,f14.0,a5,f14.0,a6,f10.0,a6,f5.2,a3)') char_dummy, SOURCE%RP%strk,      &
                                                               char_dummy, SOURCE%RP%dip,                &
                                                               char_dummy, SOURCE%RP%rake_ave,           &
                                                               char_dummy, SOURCE%RP%Htop,               &
                                                               char_dummy
       READ(IO%UNIT%other01,'(a15,f8.3,a12,f7.3,a12,f5.1,a20)')char_dummy, SOURCE%RP%HypX,               &
                                                               char_dummy, SOURCE%RP%HypZ,               &
                                                               char_dummy, SOURCE%RP%RiseTime_ave,       &
                                                               char_dummy
       DO i = 1,3
          READ(IO%UNIT%other01,*)               ! Read unimportant comments
       ENDDO
       READ(IO%UNIT%other01,'(a16,i16,a5,i13,a32)')            char_dummy, SOURCE%RP%nxRP,               &
                                                               char_dummy, SOURCE%RP%nzRP,               &
                                                               char_dummy
       READ(IO%UNIT%other01,'(a16,f8.3,a13,f7.3,a42)')         char_dummy, SOURCE%RP%dxRP,               &
                                                               char_dummy, SOURCE%RP%dzRP,               &
                                                               char_dummy
       READ(IO%UNIT%other01,'(a16,i16,a5,i13,a39)')            char_dummy, SOURCE%RP%nTWindow,           &
                                                               char_dummy, SOURCE%RP%nSegments,          &
                                                               char_dummy
       READ(IO%UNIT%other01,'(a16,f7.3,a14,f7.3,a39)')         char_dummy, SOURCE%RP%TWindowLen,         &
                                                               char_dummy, SOURCE%RP%WindowShift,        &
                                                               char_dummy
       DO i = 1,9
          READ(IO%UNIT%other01,*)            ! Read unimportant comments
       ENDDO
       READ(IO%UNIT%other01,'(a9,i4)')    char_dummy, iDummy1
       READ(IO%UNIT%other01,*)               ! Read unimportant comments
       READ(IO%UNIT%other01,*)               ! Read unimportant comments
       DO i = 1,iDummy1
          READ(IO%UNIT%other01,*)            ! Read velocity field structure and comments
       ENDDO
       DO i = 1,6
          READ(IO%UNIT%other01,*)            ! Read unimportant comments
       ENDDO

       CLOSE(IO%UNIT%other01)

       ALLOCATE( SOURCE%RP%SegStrk(SOURCE%RP%nSegments),           &
                 SOURCE%RP%SegDip(SOURCE%RP%nSegments),            &
                 SOURCE%RP%SegLen(SOURCE%RP%nSegments),            &
                 SOURCE%RP%SegWid(SOURCE%RP%nSegments),            &
                 SOURCE%RP%nSbfs(SOURCE%RP%nSegments),             &
                 SOURCE%RP%Element(SOURCE%RP%nxRP*SOURCE%RP%nzRP), &
                 SOURCE%RP%SpacePosition(3,SOURCE%RP%nxRP*SOURCE%RP%nzRP) )

       logInfo(*) 'Finite Source Rupture Geometry read. '
       logInfo(*) 'Strike           :  ', SOURCE%RP%strk
       logInfo(*) 'Dip              :  ', SOURCE%RP%dip
       logInfo(*) 'Rake(average)    :  ', SOURCE%RP%rake_ave
       logInfo(*) 'RiseTime(average):  ', SOURCE%RP%RiseTime_ave
       !
       ! Different Rupture Model Formats have to be distinguished
       !
       SELECT CASE(SOURCE%RP%Type)

       !CASE(1) !############################DOES NOT WORK YET #####################

       CASE(2)
          !
          ! Total slip and rake available for each subfault and each time window
          ! Only one Rupture Plane Segment
          ! (example: s1994NORTHRhart)
          !
          READ(IO%UNIT%other01,'(a11,i7)')    char_dummy, SOURCE%RP%nSbfs(1)
          DO i = 1,5
               READ(IO%UNIT%other01,*)            ! Read unimportant comments
          ENDDO

          logInfo(*) 'Reading        ',SOURCE%RP%nSbfs(1),' subfaults.'

          ALLOCATE ( dummy(5+SOURCE%RP%nTWindow) )
          ALLOCATE ( SOURCE%RP%Slip(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )
          ALLOCATE ( SOURCE%RP%Sliprate(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )
          ALLOCATE ( SOURCE%RP%Rake(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )

          DO j = 1,SOURCE%RP%nSbfs(1)
               READ(IO%UNIT%other01,*) dummy
               SOURCE%RP%Slip(j,:) = dummy(6:5+SOURCE%RP%nTWindow)
               SOURCE%RP%Rake(j,:) = dummy(5)
          ENDDO

       !CASE(3) !############################ DOES NOT WORK YET #####################
          !
          ! Strike-slip and Dip-slip available for each subfault
          ! Rupture Plane Segments with Subfault coordinates and slip values
          ! of all time windows are read now
          ! (example: s1995KOBEJAseki)
          !

          !ALLOCATE ( dummy(4+2*SOURCE%RP%nTWindow) )
          !ALLOCATE ( SOURCE%RP%Slip(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )
          !ALLOCATE ( SOURCE%RP%Sliprate(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )

          !k = 0            !subfault counter

          !DO i = 1,SOURCE%RP%nSegments

          !    READ(IO%UNIT%other01,*)           ! Read unimportant comment
          !    READ(IO%UNIT%other01,'(a26,f7.0,a15,f7.1,a3)')      char_dummy, SOURCE%RP%SegStrk(i),     &
          !                                                        char_dummy, SOURCE%RP%SegDip(i),      &
          !                                                        char_dummy
          !    READ(IO%UNIT%other01,'(a26,f8.3,a14,f8.3,a3)')      char_dummy, SOURCE%RP%SegLen(i),      &
          !                                                        char_dummy, SOURCE%RP%SegWid(i),      &
          !                                                        char_dummy
          !    READ(IO%UNIT%other01,'(a11,i7)')                    char_dummy, SOURCE%RP%nSbfs(i)
          !    READ(IO%UNIT%other01,*)           ! Read unimportant comment
          !    READ(IO%UNIT%other01,*)           ! Read unimportant comment
          !    READ(IO%UNIT%other01,*)           ! Read unimportant comment
          !
          !    DO j = 1,SOURCE%RP%nSbfs(i)
          !         k = k+1
          !         READ(IO%UNIT%other01,*) SOURCE%RP%Slip(k,:)
          !    ENDDO
          !
          !ENDDO

       CASE(3,4)
          !
          ! Rise time and onset times available for each subfault
          ! Strike-slip and Dip-slip  available for each subfault
          ! Only one Rupture Plane Segment with Subfault coordinates
          ! Only works for one time window
          ! (example: s2005BLINDTkaes)
          !
          READ(IO%UNIT%other01,'(a11,i7)')    char_dummy, SOURCE%RP%nSbfs(1)
          DO i = 1,5
               READ(IO%UNIT%other01,*)          ! Read unimportant comments
          ENDDO

          logInfo(*) 'Reading        ',SOURCE%RP%nSbfs(1),' subfaults.'

          ALLOCATE ( dummy(6) )
          ALLOCATE ( SOURCE%RP%SSlip(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )
          ALLOCATE ( SOURCE%RP%DSlip(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )
          ALLOCATE ( SOURCE%RP%Sliprate(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )
          ALLOCATE ( SOURCE%RP%Rake(SOURCE%RP%nxRP*SOURCE%RP%nzRP,SOURCE%RP%nTWindow) )
          ALLOCATE ( SOURCE%RP%Tonset(SOURCE%RP%nxRP*SOURCE%RP%nzRP) )
          ALLOCATE ( SOURCE%RP%TRise(SOURCE%RP%nxRP*SOURCE%RP%nzRP) )

          DO j = 1,SOURCE%RP%nSbfs(1)
               READ(IO%UNIT%other01,*) dummy
               SOURCE%RP%Tonset(j)  = dummy(3)
               SOURCE%RP%TRise(j)   = dummy(4)
               SOURCE%RP%SSlip(j,1) = dummy(5)
               SOURCE%RP%DSlip(j,1) = dummy(6)
               SOURCE%RP%Rake(j,1)  = ATAN2(-dummy(6),dummy(5))
               IF(SOURCE%RP%Rake(j,1).LT.0.)THEN
                    SOURCE%RP%Rake(j,1) = (2*pi+SOURCE%RP%Rake(j,1))/pi*180.
               ELSE
                    SOURCE%RP%Rake(j,1) = SOURCE%RP%Rake(j,1)/pi*180.
               ENDIF
          ENDDO

       CASE DEFAULT
          logError(*)  'The format type of the Finite Source Rupture Model is unknown! '
          STOP                                                                                     ! STOP

       END SELECT

    CASE(31,32,40,41)

       ALLOCATE ( SOURCE%RP%nSbfs(1) )

       logInfo(*) 'Finite Source Rupture Model (FREE-FORMAT) chosen. '
       SOURCE%FSRMFileName = FileName
       logInfo(*) 'Source term read from ', TRIM(SOURCE%FSRMFileName)
       CALL OpenFile(                                       &
            UnitNr       = IO%UNIT%other01                , &
            Name         = SOURCE%FSRMFileName            , &
            create       = .FALSE.                          )
       logInfo(*) 'Reading rupture model file ...  '
       !
       ! Moment Tensor, Rupture Plane Geometrie and Subfault information is read now
       !
       READ(IO%UNIT%other01,*)                                           ! Read comment
       SOURCE%RP%MomentTensor(:,:) = 0.                                  !
       READ(IO%UNIT%other01,*) SOURCE%RP%MomentTensor(1,:)               ! Read Moment Tensor
       READ(IO%UNIT%other01,*) SOURCE%RP%MomentTensor(2,:)               ! Read Moment Tensor
       READ(IO%UNIT%other01,*) SOURCE%RP%MomentTensor(3,:)               ! Read Moment Tensor
       READ(IO%UNIT%other01,*)                                           ! Read comment
       READ(IO%UNIT%other01,*) SOURCE%RP%nSbfs(1)                        ! Read number of subfaults
       READ(IO%UNIT%other01,*)                                           ! Read comment
       ALLOCATE ( dummy(8) )
       ALLOCATE ( SOURCE%RP%SpacePosition(3,SOURCE%RP%nSbfs(1)), &
                  SOURCE%RP%Strks(1,SOURCE%RP%nSbfs(1)),         &
                  SOURCE%RP%Dips(1,SOURCE%RP%nSbfs(1)),          &
                  SOURCE%RP%Rake(1,SOURCE%RP%nSbfs(1)),          &
                  SOURCE%RP%Tonset(SOURCE%RP%nSbfs(1)),          &
                  SOURCE%RP%Area(SOURCE%RP%nSbfs(1)),          &
                  SOURCE%RP%Element(SOURCE%RP%nSbfs(1))  )
       !
       DO i = 1,SOURCE%RP%nSbfs(1)
           READ(IO%UNIT%other01,*) dummy                                 ! Read subfault data
           SOURCE%RP%SpacePosition(1:3,i) = dummy(1:3)
           SOURCE%RP%Strks(1,i)           = dummy(4)
           SOURCE%RP%Dips(1,i)            = dummy(5)
           SOURCE%RP%Rake(1,i)            = dummy(6)
           SOURCE%RP%Area(i)              = dummy(7)
           SOURCE%RP%Tonset(i)            = dummy(8)
       ENDDO

       CLOSE(IO%UNIT%other01)

       DEALLOCATE ( dummy )
       !
       !

    case(42) ! Netcdf rupture format
      logInfo(*) 'Netcdf rupture format chosen.'
      SOURCE%NRFFileName = FileName
#ifndef USE_NETCDF
      logError(*) 'NRF sources require netcdf support.'
      stop
#endif

    CASE(50) !Finite sources with individual slip rate history for each subfault

       ALLOCATE ( SOURCE%RP%nSbfs(1) )

       logInfo(*) 'Finite Source Rupture Model (FREE-FORMAT) chosen. '
       SOURCE%FSRMFileName = FileName
       logInfo(*) 'Source term read from ', TRIM(SOURCE%FSRMFileName)
       CALL OpenFile(                                       &
            UnitNr       = IO%UNIT%other01                , &
            Name         = SOURCE%FSRMFileName            , &
            create       = .FALSE.                          )
       logInfo(*) 'Reading rupture model file ...  '
       !
       ! Moment Tensor, Rupture Plane Geometrie and Subfault information is read now
       !
       READ(IO%UNIT%other01,*)                                           ! Read comment
       SOURCE%RP%MomentTensor(:,:) = 0.                                  !
       READ(IO%UNIT%other01,*) SOURCE%RP%MomentTensor(1,:)               ! Read Moment Tensor
       READ(IO%UNIT%other01,*) SOURCE%RP%MomentTensor(2,:)               ! Read Moment Tensor
       READ(IO%UNIT%other01,*) SOURCE%RP%MomentTensor(3,:)               ! Read Moment Tensor
       READ(IO%UNIT%other01,*)                                           ! Read comment
       READ(IO%UNIT%other01,*) SOURCE%RP%nSbfs(1)                        ! Read number of subfaults
       READ(IO%UNIT%other01,*)                                           ! Read comment
       ALLOCATE ( dummy(8) )
       ALLOCATE ( SOURCE%RP%SpacePosition(3,SOURCE%RP%nSbfs(1)), &
                  SOURCE%RP%Strks(1,SOURCE%RP%nSbfs(1)),         &
                  SOURCE%RP%Dips(1,SOURCE%RP%nSbfs(1)),          &
                  SOURCE%RP%Rake(1,SOURCE%RP%nSbfs(1)),          &
                  SOURCE%RP%Tonset(SOURCE%RP%nSbfs(1)),          &
                  SOURCE%RP%Area(SOURCE%RP%nSbfs(1)),            &
                  SOURCE%RP%Element(SOURCE%RP%nSbfs(1))  )

       DO i = 1,SOURCE%RP%nSbfs(1)
           READ(IO%UNIT%other01,*) dummy                                 ! Read subfault data
           SOURCE%RP%SpacePosition(1:3,i) = dummy(1:3)
           SOURCE%RP%Strks(1,i)           = dummy(4)
           SOURCE%RP%Dips(1,i)            = dummy(5)
           SOURCE%RP%Rake(1,i)            = dummy(6)
           SOURCE%RP%Area(i)              = dummy(7)
           SOURCE%RP%Tonset(i)            = dummy(8)
       ENDDO

       DEALLOCATE ( dummy )

       READ(IO%UNIT%other01,*)                                           ! Read comment
       READ(IO%UNIT%other01,*) SOURCE%RP%t_samp, SOURCE%RP%nsteps        ! Read time sampling of rupture functions

       SOURCE%RP%T_max = SOURCE%RP%t_samp * (SOURCE%RP%nsteps-1)         ! Equispaced time sampling for all subfaults assumed

       ALLOCATE ( SOURCE%RP%TimeHist(SOURCE%RP%nsteps,SOURCE%RP%nSbfs(1))  )

       READ(IO%UNIT%other01,*)                                           ! Read comment
       DO i = 1,SOURCE%RP%nSbfs(1)
           DO j = 1,SOURCE%RP%nsteps
               READ(IO%UNIT%other01,*) SOURCE%RP%TimeHist(j,i)           !
           ENDDO
       ENDDO

       CLOSE(IO%UNIT%other01)
       !
    CASE DEFAULT                                                                                   !
       logError(*)  'The sourctype specified (', SOURCE%Type, ') is unknown! '                  !
       STOP                                                                                        ! STOP
    END SELECT                                                                                     !

                                                                                                   !
    IF(EQN%Adjoint.EQ.1) THEN  !verify compliance of sources to adjoint simulations and sort them by location
      !
      SOURCE%Ricker%nRicker = nRicker
      SOURCE%TimeGP%nPulseSource = nPulseSource
      IF(SOURCE%Ricker%nRicker.EQ.0.AND.SOURCE%TimeGP%nPulseSource.EQ.0) THEN
        logError(*)  'Adjoint simulations require point sources of type 16, 18 or 19! '                  !
        STOP
      ENDIF
      !
      count=0
      !
      IF(SOURCE%Ricker%nRicker.NE.0) THEN
        nsrc = SOURCE%Ricker%nRicker
        ALLOCATE(loc_src(3,nsrc))
        loc_src(:,1:nsrc) = SOURCE%Ricker%SpacePosition(:,1:nsrc)
      ELSEIF(SOURCE%TimeGP%nPulseSource.NE.0) THEN
        nsrc = SOURCE%TimeGP%nPulseSource
        ALLOCATE(loc_src(3,nsrc))
        loc_src(:,1:nsrc) = SOURCE%TimeGP%SpacePosition(:,1:nsrc)
      ENDIF
      !
      !Group sources at the same location
      IF(nsrc.GE.1) THEN
        ALLOCATE(Dist(nsrc,nsrc),Sorted(nsrc),Group(nsrc),Members(nsrc,nsrc))
        DO i=1,nsrc-1
          DO j=i+1,nsrc
            Dist(i,j)=SQRT((loc_src(1,i)-loc_src(1,j))**2 + (loc_src(2,i)-loc_src(2,j))**2 + (loc_src(3,i)-loc_src(3,j))**2 )
          ENDDO
        ENDDO
        Sorted  = 0
        Group   = 0
        Members = 0
        count=0
        DO i=1,nsrc-1
          IF(Sorted(i).EQ.0) THEN
            Sorted(i)= 1
            count=count+1
            Group(count)=1
            Members(count,1)=i
            count2=1
            DO j=i+1,nsrc
              IF(Dist(i,j).LE.1.0e-20) THEN
                Sorted(j)=1
                count2=count2+1
                Group(count)=count2
                Members(count,count2)=j
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IF(Sorted(nsrc).EQ.0) THEN
          Sorted(i)= 1
          count=count+1
          Group(nsrc)=1
          Members(count,1)=nsrc
        ENDIF
        !
        SOURCE%Adjsource%nFwdSourceLoc   = count
        SOURCE%Adjsource%MaxFwdSrcPerLoc = MAXVAL(Group(:))
        ALLOCATE(SOURCE%Adjsource%FwdSourceLoc(count,SOURCE%Adjsource%MaxFwdSrcPerLoc))
        !
        SOURCE%Adjsource%FwdSourceLoc = 0
        !
        DO i=1,count
          DO j=1,SOURCE%Adjsource%MaxFwdSrcPerLoc
            !IF(Members(i,j).EQ.j) THEN
              SOURCE%Adjsource%FwdSourceLoc(i,j)= Members(i,j)
            !ENDIF
          ENDDO
        ENDDO

      ENDIF
    ENDIF

  END SUBROUTINE readpar_sourceterm

 SUBROUTINE readsource110(IO, U0, l1)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    INTENT(IN)                             :: IO
    REAL                                   :: U0(3), l1(3)
    INTEGER                                :: readStat
    NAMELIST                               /Source110/ U0, l1
    !-----------------------------------------------------------------------

    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Source110) ! Write in namelistfile U0(1) = ... and in the next line U0(2) = ...
                                                  ! and the same for l1
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Source110")
    ENDIF
  END SUBROUTINE

 SUBROUTINE readsource15(IO, nDirac, SpacePositionx, SpacePositiony, SpacePositionz, TimePosition, Intensity, EqnNr)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    INTENT(IN)                             :: IO
    INTEGER                                :: nDirac
    INTEGER                                :: readStat
    REAL, DIMENSION(:), ALLOCATABLE        :: SpacePositionx, SpacePositiony, SpacePositionz, TimePosition, &
                                              Intensity, EqnNr
    NAMELIST                               /Source15/SpacePositionx, SpacePositiony, SpacePositionz, &
                                                     TimePosition, Intensity, EqnNr
    !-----------------------------------------------------------------------
ALLOCATE( SpacePositionx(nDirac), &
          SpacePositiony(nDirac), &
          SpacePositionz(nDirac), &
          TimePosition(nDirac),    &
          Intensity(nDirac),       &
          EqnNr(nDirac))

    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Source15) ! Write in namelistfile SpacePositionx(1) = ... and in the next line SpacePositionx(2) = ...
                                                  ! and the same for SpacePositiony, SpacePositionz,...
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Source15")
    ENDIF
  END SUBROUTINE

 SUBROUTINE readsource1618(IO, nRicker, SpacePositionx, SpacePositiony, SpacePositionz, Delay, a1, f, EqnNr)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    INTENT(IN)                             :: IO
    INTEGER                                :: nRicker
    INTEGER                                :: readStat
    REAL, DIMENSION(:), ALLOCATABLE        :: SpacePositionx, SpacePositiony, SpacePositionz, Delay, &
                                              a1, f, EqnNr
    NAMELIST                               /Source1618/ SpacePositionx, SpacePositiony, SpacePositionz, Delay, &
                                                        a1, f, EqnNr
    !-----------------------------------------------------------------------
     ALLOCATE( SpacePositionx(nRicker), &
               SpacePositiony(nRicker), &
               SpacePositionz(nRicker), &
               Delay(nRicker),           &
               a1(nRicker),              &
               f(nRicker),               &
               EqnNr(nRicker))

    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Source1618) ! Write in namelistfile SpacePositionx(1) = ... and in the next line SpacePositionx(2) = ...
                                                  ! and the same for SpacePositiony, ...
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Source1618")
    ENDIF
  END SUBROUTINE

 SUBROUTINE readsource17(EQN, IO, U0, l1, l2, T)
    IMPLICIT NONE
    TYPE (tEquations)                      :: EQN
    TYPE (tInputOutput)                    :: IO
    INTENT(INOUT)                          :: EQN
    INTENT(IN)                             :: IO
    REAL, DIMENSION(:), ALLOCATABLE        :: U0, l1, l2, T
    INTEGER                                :: readStat
    NAMELIST                               /Source17/ U0, l1, l2, T
    !-----------------------------------------------------------------------
    ALLOCATE(U0(EQN%nVar), &
             l1(EQN%nVar), &
             l2(EQN%nVar), &
             T(EQN%nVar))

    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Source17) ! Write in namelistfile U0(1) = ... and in the next line U0(2) = ...
                                                  ! and the same for l1, l2, ...
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Source17")
    ENDIF
  END SUBROUTINE

 SUBROUTINE readsource19(IO, nPulseSource, EqnNr, SpacePositionx, SpacePositiony, SpacePositionz, t0, Width, A0)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    INTENT(IN)                          :: IO
    INTEGER                                :: nPulseSource
    INTEGER                                :: readStat
    REAL, DIMENSION(:), ALLOCATABLE        :: EqnNr, SpacePositionx, SpacePositiony, SpacePositionz, t0, Width, A0
    NAMELIST                               /Source19/ SpacePositionx, SpacePositiony, SpacePositionz, t0, Width, A0
    !----------------------------------------------------------------------
    ALLOCATE(EqnNr(nPulseSource), &
             SpacePositionx(nPulseSource), &
             SpacePositiony(nPulseSource), &
             SpacePositionz(nPulseSource), &
             t0(nPulseSource), &
             Width(nPulseSource), &
             A0(nPulseSource))

    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Source19) ! Write in namelistfile EqnNr(1) = ... and in the next line EqnNr(2) = ...
                                                  ! and the same for Spacepositionx, ...
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Source19")
    ENDIF
  END SUBROUTINE
  !
  !============================================================================
  ! S P O N G E   L A Y E R
  !============================================================================

  SUBROUTINE readpar_spongelayer(DISC,EQN,SOURCE,IO)
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tSource)             :: Source
    TYPE (tDiscretization)     :: DISC
    TYPE (tInputOutput)        :: IO
    INTEGER                    :: iSponge
    INTEGER                    :: intDummy
    INTEGER                    :: readStat
    !--------------------------------------------------------------------------
    INTENT(INOUT)              :: SOURCE, DISC,EQN
    INTENT(IN)                 :: IO
    !--------------------------------------------------------------------------
    INTEGER                          :: enabled, nDGSponge
    REAL,DIMENSION(:),ALLOCATABLE    :: SpongeDelta, SpongePower, SigmaMax
    REAL                             :: DGSpongeTol, PMLDelta, Refl_Coeff, &
                                        PMLPrefactor, PMLFrequency
    NAMELIST                        /SpongeLayer/ enabled, DGSpongeTol, nDGSponge, &
                                                  intDummy, PMLDelta, Refl_Coeff, &
                                                  PMLPrefactor, PMLFrequency
    !------------------------------------------------------------------------
    !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<  S P O N G E    L A Y E R                               >'
    logInfo(*) '<--------------------------------------------------------->'
    !
    !Setting default values
    enabled = 0
    !
    READ (IO%UNIT%FileIn, IOSTAT=readStat, nml = SpongeLayer)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "SpongeLayer")
    ENDIF
    SOURCE%Sponge%enabled = enabled

    SELECT CASE(SOURCE%Sponge%enabled)
    !
    CASE(0)
       logInfo(*)  'No sponge layer used '
       !
    CASE(1)
       logInfo(*)  'A sponge layer for general domains is used '
       DISC%Galerkin%DGSpongeTol = DGSpongeTol
       DISC%Galerkin%nDGSponge = nDGSponge
       logInfo(*)  'Number of sponge boundaries specified: ' , DISC%Galerkin%nDGSponge
       ALLOCATE(DISC%Galerkin%DGSponge(DISC%Galerkin%nDGSponge) )
       call readSponges(IO, nDGSponge, SpongeDelta, SpongePower, SigmaMax)
       DO iSponge = 1, DISC%Galerkin%nDGSponge
          DISC%Galerkin%DGSponge(iSponge)%SpongeObject = intDummy/100
          DISC%Galerkin%DGSponge(iSponge)%SpongeBND    = MOD(intDummy,100)
          DISC%Galerkin%DGSponge(iSponge)%SpongeDelta = SpongeDelta(iSponge)
          DISC%Galerkin%DGSponge(iSponge)%SpongePower = SpongePower(iSponge)
          DISC%Galerkin%DGSponge(iSponge)%SigmaMax = SigmaMax(iSponge)
          logInfo(*) 'Sponge Nr : ', iSponge
          logInfo(*) 'Object, Type, Width, Power: ',         &
                                  DISC%Galerkin%DGSponge(iSponge)%SpongeObject, &
                                  DISC%Galerkin%DGSponge(iSponge)%SpongeBND,    &
                                  DISC%Galerkin%DGSponge(iSponge)%SpongeDelta,  &
                                  DISC%Galerkin%DGSponge(iSponge)%SpongePower
       ENDDO
    CASE(3)
        logInfo(*)  'A PML layer for general domains is used '
        DISC%Galerkin%PMLayer%PMLDelta = PMLDelta
        DISC%Galerkin%PMLayer%Refl_Coeff = Refl_Coeff
        DISC%Galerkin%PMLayer%PMLPrefactor = PMLPrefactor
        DISC%Galerkin%PMLayer%PMLFrequency = PMLFrequency
        ! Set number of variables for CPML
        EQN%nVar = 13
        EQN%nVarTotal = 13
        !
    CASE default
       !WRITE(IO%UNIT%errOut,*)  '|   The option is not valid! 0 '        , &
       !     'disables, 1 enables the sponge layer, 2 takes PML and 3 CPML '
       !STOP
       !
    END SELECT
    !
  END SUBROUTINE readpar_spongelayer

  SUBROUTINE readSponges(IO, nDGSponge, SpongeDelta, SpongePower, SigmaMax)
    IMPLICIT NONE
    TYPE (tInputOutput)                    :: IO
    INTENT(IN)                             :: IO
    INTEGER                                :: nDGSponge
    INTEGER                                :: readStat
    REAL, DIMENSION(:), ALLOCATABLE        :: SpongeDelta, SpongePower, SigmaMax
    NAMELIST                               /Sponges/ SpongeDelta, SpongePower, SigmaMax
    !----------------------------------------------------------------------
       ALLOCATE(SpongeDelta(nDGSponge), &
                SpongePower(nDGSponge), &
                SigmaMax(nDGSponge))

    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Sponges) ! Write in namelistfile SpongeDelta(1) = ... and in the next line SpongeDelta(2) = ...
                                                  ! and the same for SpongePower, ...
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Sponges")
    ENDIF
   END SUBROUTINE
  !
  !============================================================================
  ! M E S H
  !============================================================================

  SUBROUTINE readpar_mesh(EQN,IC,MESH,DISC,BND,SOURCE,IO)
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tInitialCondition)   :: IC
    TYPE (tUnstructMesh)       :: MESH
    Type (tDiscretization)     :: DISC
    TYPE (tBoundary)           :: BND
    TYPE (tSource)             :: SOURCE
    TYPE (tInputOutput)        :: IO
    REAL                       :: periodic_direction(3)
    INTEGER                    :: j ,k
    INTEGER                    :: i, stat
    INTEGER                    :: readStat
    CHARACTER(LEN=600)          :: Name
    LOGICAL                    :: file_exits
    !------------------------------------------------------------------------
    INTENT(INOUT)              :: MESH,BND,SOURCE,IO
    INTENT(IN)                 :: IC
    !------------------------------------------------------------------------
    INTEGER                          :: periodic
    REAL                             :: ScalingMatrixX(3), ScalingMatrixY(3), ScalingMatrixZ(3), &
                                        displacement(3)
    CHARACTER(LEN=600)               :: MeshFile, meshgenerator
    NAMELIST                         /MeshNml/ MeshFile, meshgenerator, periodic, &
                                            periodic_direction, displacement, ScalingMatrixX, &
                                            ScalingMatrixY, ScalingMatrixZ
    !------------------------------------------------------------------------
    !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<  M E S H                                                >'
    logInfo(*) '<--------------------------------------------------------->'
    ! Put a redundant dimension information also into the MESH data structure
    !
    MESH%Dimension = EQN%Dimension
    !
    MESH%MESHversionID = 0.0

    ! Setting default values
    MeshFile = 'LOH1'
    meshgenerator = 'Gambit3D-fast'
    displacement(:) = 0.
    ScalingMatrixX(:) = 0.0
    ScalingMatrixX(1) = 1.0
    ScalingMatrixY(:) = 0.0
    ScalingMatrixY(2) = 1.0
    ScalingMatrixZ(:) = 0.0
    ScalingMatrixZ(3) = 1.0
    periodic = 0
    periodic_direction(:) = 0
    !
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = MeshNml)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "MeshNml")
    ENDIF

    IO%MeshFile = MeshFile                               ! mesh input (mesh file name, no_file)

    Name = TRIM(IO%MeshFile) // '.met'
    IO%MetisFile = Name(1:600)

    MESH%iniSquareMesh = .FALSE.
    MESH%iniDiscMesh   = .FALSE.

    IO%meshgenerator = trim(meshgenerator)

       EQN%HexaDimension = 3
       SELECT CASE(IO%meshgenerator)
       CASE('Gambit3D-fast','Netcdf','PUML')
          if (IO%meshgenerator .eq. 'Netcdf') then
            logInfo0(*) 'Read a netCDF mesh ...'
            Name = trim(IO%MeshFile) // '.nc'
          elseif (IO%meshgenerator .eq. 'PUML') then
#if defined(USE_METIS) && defined(USE_HDF) && defined(USE_MPI)
            Name = trim(IO%MeshFile)
            logInfo0(*) 'Read a PUML mesh file'
#else
            logError(*) 'PUML requires METIS and HDF5'
#endif
          else
            logInfo0(*) 'Read a Gambit 3-D neutral mesh ... '
            Name = TRIM(IO%MeshFile)//'.neu'
          endif

          IO%MeshFile=Name(1:600)

          inquire( file=IO%MeshFile , exist=file_exits )
          if ( .NOT.file_exits ) then
             logError(*) 'mesh file ',IO%MeshFile,'does not exists'
             STOP
          endif

          !
          logInfo(*) 'Mesh is READ from file      ',    IO%MeshFile
          !
          BND%periodic = periodic
          !
          IF(BND%periodic.EQ.0) THEN
             BND%DirPeriodic(:) = .FALSE.
             logInfo(*) 'No periodic boundary conditions specified.    '
          ELSE
             WHERE(periodic_direction(:).EQ.1)
                BND%DirPeriodic(:) = .TRUE.
             ELSEWHERE
                BND%DirPeriodic(:) = .FALSE.
             ENDWHERE
             BND%Periodic = SUM(periodic_direction)
          ENDIF
          !
          IF(BND%DirPeriodic(1)) THEN
            logInfo(*) 'Periodic boundary in x-direction. '
          ENDIF
          IF(BND%DirPeriodic(2)) THEN
            logInfo(*) 'Periodic boundary in y-direction. '
          ENDIF
          IF(BND%DirPeriodic(3)) THEN
            logInfo(*) 'Periodic boundary in z-direction. '
          ENDIF

       CASE DEFAULT
          logError(*) 'Meshgenerator ', TRIM(IO%meshgenerator), ' is unknown!'
          STOP
       END SELECT
    ! specify element type (3-d = tetrahedrons)

      IF(IO%meshgenerator.eq.'Gambit3D-fast' .or. IO%meshgenerator.eq.'Netcdf' .or. IO%meshgenerator.eq.'PUML') THEN
          MESH%GlobalElemType = 4
          MESH%GlobalSideType = 3
          MESH%GlobalVrtxType = 4
          MESH%nVertexMax = 4
          MESH%nSideMax = 4
       ELSE
          logError(*) 'Wrong definition of meshgenerator.'
          STOP
       ENDIF

       SELECT CASE (MESH%GlobalElemType)
       CASE(4)
          logInfo(*) 'Mesh consits of TETRAHEDRAL elements.'
          logInfo(*) 'Mesh type is', MESH%GlobalElemType
       CASE DEFAULT
          logError(*) 'MESH%GlobalElemType must be {4}, {6} or {7} '
          STOP
       END SELECT
          !logInfo(*) 'Mesh consits of TETRAHEDRAL elements.' ! Is obvious as only this is possible
          !logInfo(*) 'Mesh type is', MESH%GlobalElemType
    !
    IF(MESH%MESHversionID.eq.0.0)THEN
         ALLOCATE(MESH%Displacement(1:EQN%Dimension), MESH%ScalingMatrix(1:EQN%Dimension,EQN%Dimension))
         MESH%Displacement(:) = displacement(:)
         logInfo(*) 'Displacement of original mesh'
         logInfo(*) '   ', MESH%Displacement(:)
         !
         MESH%ScalingMatrix(1,:) = ScalingMatrixX(:)
         MESH%ScalingMatrix(2,:) = ScalingMatrixY(:)
         MESH%ScalingMatrix(3,:) = ScalingMatrixZ(:)

         !
         logInfo(*) 'Scaling and Rotation of original mesh'
         logInfo(*) ScalingMatrixX
         logInfo(*) ScalingMatrixY
         logInfo(*) ScalingMatrixZ
    ENDIF
    !
  END SUBROUTINE readpar_mesh

  !============================================================================
  ! D I S C R E T I S A T I O N
  !============================================================================

  SUBROUTINE readpar_discretisation(EQN,MESH,DISC,SOURCE,IO)
    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tUnstructMesh)       :: MESH
    TYPE (tDiscretization)     :: DISC
    TYPE (tSource)             :: SOURCE
    TYPE (tInputOutput)        :: IO
    ! localVariables
    INTEGER                    :: intDummy, stat, i
    INTEGER                    :: readStat
    CHARACTER(LEN=5)           :: cInput
    CHARACTER(LEN=300)         :: cDummy

    !------------------------------------------------------------------------
    INTENT(IN   )              :: SOURCE
    INTENT(INOUT)              :: IO, EQN
    !------------------------------------------------------------------------
    INTEGER                          :: DGFineOut1D, DGMethod, ClusteredLTS, CKMethod, &
                                        FluxMethod, IterationCriterion, nPoly, nPolyRec, &
                                        StencilSecurityFactor, LimiterSecurityFactor, &
                                        Order, Material, nPolyMap
    REAL                             :: CFL, FixTimeStep
    NAMELIST                         /Discretization/ DGFineOut1D, DGMethod, ClusteredLTS, &
                                                      CKMethod, FluxMethod, IterationCriterion, &
                                                      nPoly, nPolyRec, &
                                                      LimiterSecurityFactor, Order, Material, &
                                                      nPolyMap, CFL, FixTimeStep
    !------------------------------------------------------------------------
    !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<  D I S C R E T I S A T I O N                            >'
    logInfo(*) '<-------------------------------------------------------- >'
    logInfo(*) 'Discontinuous Galerkin technique is used. '
    DISC%Galerkin%ZoneOrderFlag = 0 ! aheineck, this is used but never set, but we need to init it

    ! Setting default values
    DGFineOut1D = 0
    CKMethod = 0
    FluxMethod = 0
    DGMethod = 1
    ! 0: read from file, 1: GTS, 2-n: multi-rate
    ClusteredLTS = 1
    CFL = 0.5
    nPolyMap = 0                                                               !                                                                  !
    Material = 1
    FixTimeStep = 5000
    !                                                              ! DGM :
    READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Discretization)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Discretization")
    ENDIF
    DISC%Galerkin%DGFineOut1D = DGFineOut1D                        ! No. of red-refinements
    !                                                              ! for 2-D fine output
    IF(DISC%Galerkin%DGFineOut1D.GT.0) THEN
      logInfo(*) 'Fine output for DG method required. '
    ELSE
      logInfo(*) 'Fine output for DG method not required. '
    ENDIF
    !
    ! =========== NEW ordering of DGMethod and Order ===============================
    !
    disc%galerkin%clusteredLts = ClusteredLts
    select case( disc%galerkin%clusteredLts )
    case(0)
      logError(*) 'TODO: Using clustered LTS with clustering provided file input'
    case(1)
      logInfo(*) 'Using GTS'
    case default
      logInfo(*) 'Using multi-rate clustered LTS:', disc%galerkin%clusteredLts
    endselect

    DISC%Galerkin%DGMethod = DGMethod
    DISC%Galerkin%CKMethod = CKMethod    ! Default: standard CK procedure (0)
    !
    SELECT CASE(DISC%Galerkin%CKMethod)
    CASE(0)
     logInfo(*) 'Using standard CK procedure. '
    CASE(1)
     logInfo(*) 'Using the space-time DG approach. '
    END SELECT
    !
    DISC%Galerkin%FluxMethod = FluxMethod
    !
    SELECT CASE(DISC%Galerkin%FluxMethod)
     CASE(0)
      logInfo(*) 'Using Godunov flux. '
     CASE(1)
      logInfo(*) 'Using Rusanov flux. '
     CASE DEFAULT
      logError(*) 'Flux case not defined ! '
      STOP
     ENDSELECT
    !
    SELECT CASE(DISC%Galerkin%DGMethod)
    CASE(1)
           logInfo(*) 'ADER-DG with global time stepping is used.'
    CASE(3)
           logInfo(*) 'ADER-DG with local timestepping is used.'
           DISC%IterationCriterion = IterationCriterion
           SELECT CASE(DISC%IterationCriterion)
           CASE(1)
               logInfo(*) 'One iteration is defined as one cycle. '
           CASE(2)
               logInfo(*) 'One iteration is defined by the update of all elements. '
           END SELECT
    CASE DEFAULT
         logError(*) 'Wrong DGmethod. Must be 1 or 3.!'
         STOP
    END SELECT
       !
    SELECT CASE(DISC%Galerkin%DGMethod)
    CASE(1,3)
           DISC%SpaceOrder = Order
#if GENERATEDKERNELS
           if (DISC%SpaceOrder .ne. CONVERGENCE_ORDER) then
                logWarning0(*) 'Ignoring min space order from parameter file, using', CONVERGENCE_ORDER
           endif
           DISC%SpaceOrder = CONVERGENCE_ORDER
#endif
           DISC%Galerkin%nMinPoly = DISC%SpaceOrder - 1

#if GENERATEDKERNELS
           if (DISC%SpaceOrder .ne. CONVERGENCE_ORDER) then
                logWarning0(*) 'Ignoring space order from parameter file, using', CONVERGENCE_ORDER
           endif
           DISC%SpaceOrder = CONVERGENCE_ORDER
#endif
           DISC%Galerkin%nPoly    = DISC%SpaceOrder - 1

             ! The choice for p-adaptivity is not possible anymore
             DISC%Galerkin%pAdaptivity = 0
             logInfo(*) 'No p-Adaptivity used. '
             logInfo(*) 'Basis functions degree:',DISC%Galerkin%nPoly

           DISC%Galerkin%nPolyRec = DISC%Galerkin%nPoly
           DISC%Galerkin%nPolyMat       = 0
           DISC%Galerkin%nDegFrMat      = 1
           DISC%Galerkin%nPolyMatOrig = Material
           DISC%Galerkin%nPolyMatOrig = DISC%Galerkin%nPolyMatOrig - 1
           logInfo(*) 'Material basis functions degree:',DISC%Galerkin%nPolyMatOrig
           DISC%Galerkin%nPolyMap = nPolyMap
           DISC%Galerkin%nPolyMat  = DISC%Galerkin%nPolyMatOrig + DISC%Galerkin%nPolyMap
           DISC%Galerkin%nDegFrMat = (DISC%Galerkin%nPolyMat+1)*(DISC%Galerkin%nPolyMat+2)*(DISC%Galerkin%nPolyMat+3)/6
           IF(DISC%Galerkin%nPolyMat.GT.DISC%Galerkin%nPoly) THEN
             logError(*) 'nPolyMat larger than nPoly. '
             STOP
           ENDIF

           IF(MESH%GlobalElemType.EQ.6) THEN
                  READ(IO%UNIT%FileIn,*) DISC%Galerkin%nPolyMatOrig
                  READ(IO%UNIT%FileIn,*) DISC%Galerkin%nPolyMap
                  DISC%Galerkin%nPolyMat  = DISC%Galerkin%nPolyMatOrig + DISC%Galerkin%nPolyMap
                  DISC%Galerkin%nDegFrMat = (DISC%Galerkin%nPolyMat+1)*(DISC%Galerkin%nPolyMat+2)*(DISC%Galerkin%nPolyMat+3)/6
                    IF(DISC%Galerkin%nPolyMat.GT.DISC%Galerkin%nPoly) THEN
                         logError(*) 'nPolyMat larger than nPoly. '
                         STOP
                    ENDIF
           ENDIF

    END SELECT
    !
    DISC%CFL = CFL                               ! minimum Courant number
    logInfo(*) 'The minimum COURANT number:    ', DISC%CFL
    !
        DISC%FixTimeStep = FixTimeStep
        logInfo(*) 'Specified dt_fix            : ', DISC%FixTimeStep
        logInfo(*) 'Actual timestep is min of dt_CFL and dt_fix. '
    !
  END SUBROUTINE readpar_discretisation

  !===========================================================================
  ! O U T P U T
  !============================================================================

    SUBROUTINE readpar_output(EQN,DISC,IO,CalledFromStructCode)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------
      TYPE (tEquations)             :: EQN
      TYPE (tDiscretization)        :: DISC
      TYPE (tInputOutput)           :: IO
      LOGICAL                       :: CalledFromStructCode
      ! local Variables
      INTEGER                       :: i,n, NPTS
      INTEGER                       :: allocstat
      INTEGER                       :: iOutputMask(29)
      INTEGER                       :: idimensionMask(3)
      INTEGER                       :: readStat
      CHARACTER(LEN=620)            :: Name
      REAL,DIMENSION(:),ALLOCATABLE :: X, Y, Z
      !------------------------------------------------------------------------
      INTENT(INOUT)              :: EQN,IO
      !------------------------------------------------------------------------
      INTEGER                          :: Rotation, Format, printIntervalCriterion, &
                                          pickDtType, nRecordPoint, PGMFlag, FaultOutputFlag, &
                                          iOutputMaskMaterial(1:3), nRecordPoints, Refinement, energy_output_on, IntegrationMask(1:9), SurfaceOutput, SurfaceOutputRefinement
      REAL                             :: TimeInterval, pickdt, pickdt_energy, Interval, checkPointInterval, OutputRegionBounds(1:6), SurfaceOutputInterval
      REAL                             :: printEnergiesInterval
      CHARACTER(LEN=600)               :: OutputFile, RFileName, PGMFile, checkPointFile
      !> The checkpoint back-end is specified via a string.
      !!
      !! If none is specified, checkpoints are disabled. To use the HDF5, MPI-IO or SIONlib
      !! back-ends you need to compile SeisSol with HDF5, MPI or SIONlib respectively.
      !!
      !! @allowed_values 'posix', 'hdf5', 'mpio', 'mpio_async', 'sionlib', 'none'
      !! @warning When using an asynchronous back-end (mpio_async), you might lose
      !!  2 * checkPointInterval of your computation.
      !! @more_info https://github.com/SeisSol/SeisSol/wiki/Parameter-File
      !!
      character(LEN=64)                :: checkPointBackend
      character(LEN=64)                :: xdmfWriterBackend
      NAMELIST                         /Output/ OutputFile, Rotation, iOutputMask, iOutputMaskMaterial, &
                                                Format, Interval, TimeInterval, printIntervalCriterion, Refinement, &
                                                pickdt, pickDtType, RFileName, PGMFlag, &
                                                PGMFile, FaultOutputFlag, nRecordPoints, &
                                                checkPointInterval, checkPointFile, checkPointBackend, energy_output_on, pickdt_energy, OutputRegionBounds, IntegrationMask, &
                                                SurfaceOutput, SurfaceOutputRefinement, SurfaceOutputInterval, xdmfWriterBackend, printEnergiesInterval
    !------------------------------------------------------------------------
    !
      logInfo(*) '<--------------------------------------------------------->'
      logInfo(*) '<  O U T P U T                                            >'
      logInfo(*) '<--------------------------------------------------------->'

      ! Setting default values
      OutputFile = 'data'
      iOutputMaskMaterial(:) =  0
	  IntegrationMask(:) = 0
      Rotation = 0
      Format = 10
      Refinement = 0
      pickdt = 0.1
      pickDtType = 1
      nRecordPoints = 0
      energy_output_on = 0
      pickdt_energy = 1.0
      OutputRegionBounds(:) = 0.0
!      RFileName = 'RecordPoints'
      pickDtType = 1
      PGMFlag = 0
      FaultOutputFlag = 0
      checkPointInterval = 0
      printEnergiesInterval = 0
      checkPointBackend = 'none'
#ifdef USE_HDF
      xdmfWriterBackend = 'hdf5'
#else
      xdmfWriterBackend = 'posix'
#endif
      SurfaceOutput = 0
      SurfaceOutputRefinement = 0
      SurfaceOutputInterval = 1.0e99
      !
      READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = Output)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "Output")
    ENDIF
      IO%OutputFile = OutputFile                                                   ! read output field file

      IO%OutputFile  = TRIM(IO%OutputFile)

      IO%SurfaceOutput = SurfaceOutput
      IO%SurfaceOutputRefinement = SurfaceOutputRefinement
      IO%SurfaceOutputInterval = SurfaceOutputInterval

      logInfo(*) 'Data OUTPUT is written to files '
      logInfo(*) '  ' ,IO%OutputFile

      IO%nOutputMask = 60
      ALLOCATE(IO%OutputMask(1:IO%nOutputMask), IO%TitleMask(1:IO%nOutputMask),  &
               STAT=allocStat                                            )
      IF (allocStat .NE. 0) THEN
         logError(*) 'could not allocate IO%OutputMask in readpar!'
         STOP
      END IF
      !
        IO%Rotation = Rotation
        IF(IO%Rotation.GT.0.AND.DISC%SpaceOrder.EQ.1) THEN
          logError(*) 'Space derivatives of polynomials of degree 0 cannot be computed!'
          logError(*) '   Rotations or Seismic Moment Tensor Contributions cannot be outputted!'
          logError(*) '   Increase the polynomial order or choose not to output rotational rates.'
          STOP
        ENDIF
        IF(IO%Rotation.EQ.1) THEN
          logInfo(*) 'Outputting rotational seismograms in addition to translational '                               !
        ElSEIF(IO%Rotation.EQ.2) THEN
          logInfo(*) 'Outputting moment tensor seismograms in addition to translational '                               !
        ElSEIF(IO%Rotation.EQ.3) THEN
          logInfo(*) 'Outputting curl and divergence seismograms in addition to translational '                               !
        ENDIF
      !
      IO%OutputMask = .FALSE.                                                                      !
         IO%OutputMask(:)      = .FALSE.                                                    !
         IO%OutputMask(1:3)    = .TRUE.                                                     ! x-y-z Coordinates
         IO%OutputMask(4:12)   = iOutputMask(1:9)                                           ! State vector

         IF(EQN%Anisotropy.EQ.0.AND.EQN%Poroelasticity.EQ.0.AND.EQN%Plasticity.EQ.0) THEN                           ! Isotropic material
            IO%OutputMask(13:15)  = iOutputMaskMaterial(1:3)                                      ! Constants for Jacobians
         ENDIF

         IF(EQN%Anisotropy.EQ.1.AND.EQN%Poroelasticity.EQ.0) THEN                           ! Anisotropic (triclinic) material
            IO%OutputMask(13:23)  = iOutputMask(1:11)                                       ! Constants for Jacobians
            IO%OutputMask(24:34)  = iOutputMask(1:11)                                       ! Constants for Jacobians
         ENDIF

         IF(EQN%Plasticity.EQ.1) THEN                                                       ! Plastic material properties
            IO%OutputMask(13:19)  = iOutputMask(10:16)                                      ! plastic strain output
         ENDIF

         IF(IO%Rotation.EQ.1) THEN
          ALLOCATE(IO%RotationMask(3),STAT=allocStat )                                      !
           IF (allocStat .NE. 0) THEN                                                       !
             logError(*) 'could not allocate IO%RotationMask in readpar!'!
             STOP                                                                           !
           END IF
           IO%RotationMask = .FALSE.
           IF(IO%OutputMask(10)) IO%RotationMask(1) = .TRUE.
           IF(IO%OutputMask(11)) IO%RotationMask(2) = .TRUE.
           IF(IO%OutputMask(12)) IO%RotationMask(3) = .TRUE.
         ENDIF
         !Seismic Moment Tensor Contribution
         IF(IO%Rotation.EQ.2) THEN
          ALLOCATE(IO%RotationMask(9),STAT=allocStat )                                      !
           IF (allocStat .NE. 0) THEN                                                       !
             logError(*) 'could not allocate IO%RotationMask in readpar!'!
             STOP                                                                           !
           END IF
           IO%RotationMask(1:9) = .TRUE.
         ENDIF
         !Curl/Divergence
         IF(IO%Rotation.EQ.3) THEN
          ALLOCATE(IO%RotationMask(4),STAT=allocStat )                                      !
           IF (allocStat .NE. 0) THEN                                                       !
             logError(*) 'could not allocate IO%RotationMask in readpar!'!
             STOP                                                                           !
           END IF
           IO%RotationMask(1:4) = .TRUE.
         ENDIF

      ALLOCATE(IO%OutputRegionBounds(6),STAT=allocStat )                                      !
       IF (allocStat .NE. 0) THEN                                                       !
         logError(*) 'could not allocate IO%OutputRegionBounds in readpar!'!
         STOP                                                                           !
       END IF
      IO%OutputRegionBounds(1:6) = OutputRegionBounds(1:6)
      ! Check if all are non-zero and then check if the min and max are not the same
      IF ((OutputRegionBounds(1).NE.0.0).OR.(OutputRegionBounds(2).NE.0.0).OR.&
          (OutputRegionBounds(3).NE.0.0).OR.(OutputRegionBounds(4).NE.0.0).OR.&
          (OutputRegionBounds(5).NE.0.0).OR.(OutputRegionBounds(6).NE.0.0)) THEN

          IF (OutputRegionBounds(2)-OutputRegionBounds(1) <= 0.0) THEN
              logError(*) 'Please make sure the x bounds are correct'
              STOP
          ENDIF
          IF (OutputRegionBounds(4)-OutputRegionBounds(3) <= 0.0) THEN
              logError(*) 'Please make sure the y bounds are correct'
              STOP
          ENDIF
          IF (OutputRegionBounds(6)-OutputRegionBounds(5) <= 0.0) THEN
              logError(*) 'Please make sure the z bounds are correct'
              STOP
          ENDIF
      END IF

	  ALLOCATE(IO%IntegrationMask(9),STAT=allocstat )                        !
      IF (allocStat .NE. 0) THEN                                             !
        logError(*) 'could not allocate IO%IntegrationMask in readpar!'      !
        STOP                                                                 !
      END IF
      IO%IntegrationMask(1:9) = IntegrationMask(1:9)

      IF(DISC%Galerkin%pAdaptivity.GT.0) THEN
        IO%OutputMask(59) = .TRUE.
      ENDIF
      IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        IO%OutputMask(60) = .TRUE.
      ENDIF
      !
      !                                                                        !
      IO%Format = Format                                                       ! Plot format
      !                                                                        !
      SELECT CASE(IO%Format)
      case(6)
         logInfo0(*) 'Output data is in XDMF format (new implementation)'
      case(10)
         logInfo0(*) 'Output data is disabled'
      CASE DEFAULT
         logError(*) 'print_format must be {6,10}'
         STOP
      END SELECT

      IO%TitleMask( 1) = TRIM(' "x"')
      IO%TitleMask( 2) = TRIM(' "y"')
      IF(EQN%EQType.EQ.8) THEN
                IO%TitleMask( 3) = TRIM(' "z"')
                IO%TitleMask( 4) = TRIM(' "sigma_xx"')
                IO%TitleMask( 5) = TRIM(' "sigma_yy"')
                IO%TitleMask( 6) = TRIM(' "sigma_zz"')
                IO%TitleMask( 7) = TRIM(' "sigma_xy"')
                IO%TitleMask( 8) = TRIM(' "sigma_yz"')
                IO%TitleMask( 9) = TRIM(' "sigma_xz"')
                IO%TitleMask(10) = TRIM(' "u"')
                IO%TitleMask(11) = TRIM(' "v"')
                IO%TitleMask(12) = TRIM(' "w"')

                IF(EQN%Anisotropy.EQ.0.AND.EQN%Poroelasticity.EQ.0.AND.EQN%Plasticity.EQ.0) THEN
                    IO%TitleMask(13) = TRIM(' "rho0"')
                    IO%TitleMask(14) = TRIM(' "mu"')
                    IO%TitleMask(15) = TRIM(' "lambda"')
                ENDIF
                IF(EQN%Anisotropy.EQ.1.AND.EQN%Poroelasticity.EQ.0) THEN
                    IO%TitleMask(13) = TRIM(' "rho0"')
                    IO%TitleMask(14) = TRIM(' "c11"')
                    IO%TitleMask(15) = TRIM(' "c12"')
                    IO%TitleMask(16) = TRIM(' "c13"')
                    IO%TitleMask(17) = TRIM(' "c14"')
                    IO%TitleMask(18) = TRIM(' "c15"')
                    IO%TitleMask(19) = TRIM(' "c16"')
                    IO%TitleMask(20) = TRIM(' "c22"')
                    IO%TitleMask(21) = TRIM(' "c23"')
                    IO%TitleMask(22) = TRIM(' "c24"')
                    IO%TitleMask(23) = TRIM(' "c25"')
                    IO%TitleMask(24) = TRIM(' "c26"')
                    IO%TitleMask(25) = TRIM(' "c33"')
                    IO%TitleMask(26) = TRIM(' "c34"')
                    IO%TitleMask(27) = TRIM(' "c35"')
                    IO%TitleMask(28) = TRIM(' "c36"')
                    IO%TitleMask(29) = TRIM(' "c44"')
                    IO%TitleMask(30) = TRIM(' "c45"')
                    IO%TitleMask(31) = TRIM(' "c46"')
                    IO%TitleMask(32) = TRIM(' "c55"')
                    IO%TitleMask(33) = TRIM(' "c56"')
                    IO%TitleMask(34) = TRIM(' "c66"')
                ENDIF
                IF(EQN%Plasticity.EQ.1) THEN !plastic strain output
                    IO%TitleMask(13) = TRIM(' "eps_p_xx"')
                    IO%TitleMask(14) = TRIM(' "eps_p_yy"')
                    IO%TitleMask(15) = TRIM(' "eps_p_zz"')
                    IO%TitleMask(16) = TRIM(' "eps_p_xy"')
                    IO%TitleMask(17) = TRIM(' "eps_p_yz"')
                    IO%TitleMask(18) = TRIM(' "eps_p_xz"')
                    IO%TitleMask(19) = TRIM(' "eta_p"')

                ENDIF
      ENDIF
      !
      !
!       IF(DISC%Galerkin%pAdaptivity.GT.0) THEN
!         IO%TitleMask(59) = TRIM(' "N"')
!       ENDIF
      !
      IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        IO%TitleMask(60) = TRIM(' "t"')
      ENDIF
      !
      IO%Title='VARIABLES = '
      IO%nrPlotVar = 0
      logInfo(*) 'Variables plotted: '
      logInfo(*) ' '
      DO i=1,IO%nOutputMask
         IF(IO%OutputMask(i)) THEN
            Name = TRIM(IO%Title) // TRIM(IO%TitleMask(i))
            IO%Title     = Name(1:600)
            IO%nrPlotVar = IO%nrPlotVar + 1
            logInfo(*) '  - ', TRIM(IO%TitleMask(i))
         ENDIF
      ENDDO

      logInfo(*) ' '
      !

      IO%outInterval%printIntervalCriterion = printIntervalCriterion
      !
      IF (IO%outInterval%printIntervalCriterion.EQ.1.AND.DISC%Galerkin%DGMethod.EQ.3) THEN
        logError(*) 'specifying IO%outInterval%printIntervalCriterion: '
        logError(*) 'When local time stepping is used, only Criterion 2 can be used! '
        STOP
      END IF
      IF (      IO%outInterval%printIntervalCriterion .EQ. 1 &                 !
           .OR. IO%outInterval%printIntervalCriterion .EQ. 3 ) THEN
          IO%outInterval%Interval = Interval
         logInfo0(*) 'Output data are generated '          , & !
              'every ', IO%outInterval%Interval, '. timestep'                  !
#ifdef GENERATEDKERNELS
         logError(*) 'Time step-wise output only with classic version'
         stop
#else
         logWarning0(*) 'Time step-wise output is deprecated! Your parameter file is not compatible with GK version!'
#endif
      END IF                                                                   !
      IF (      IO%outInterval%printIntervalCriterion .EQ. 2 &                 !
           .OR. IO%outInterval%printIntervalCriterion .EQ. 3 ) THEN
         IF( IO%Format .eq. 10) THEN
           ! we don't want output, so avoid confusing time stepping by setting
           ! plot interval to "infinity"
           IO%outInterval%TimeInterval = 1E99
           logInfo0(*) 'No output (FORMAT=10) specified, delta T set to: ', IO%OutInterval%TimeInterval
         ELSE
           IO%outInterval%TimeInterval = TimeInterval          !
           logInfo0(*) 'Output data are generated at delta T= ', IO%OutInterval%TimeInterval
         ENDIF
      END IF                                                                   !
      !                                                                        !
      !
      ! Initiate the total number of point to record with zero
      IO%ntotalRecordPoint = 0
      IO%nRecordPoint      = 0         ! points, where whole time series       is recorded
      IO%nPGMRecordPoint   = 0         ! points, where only Peak Ground Motion is recorded
      !
      ! Read pickdt and pickDtType
         IO%pickdt = pickdt
         IO%pickDtType = pickDtType
         IF (DISC%Galerkin%DGMethod .ne. 1 .and. IO%pickDtType .ne. 1) THEN
            logError(*) 'Pickpoint sampling every x timestep can only be used with global timesteping'
            STOP
         ENDIF
         SELECT CASE (IO%pickDtType)
         CASE (1)
         CASE (2)
#ifdef GENERATEDKERNELS
            logError(*) 'Time step-wise output only with classic version'
            stop
#else
            logWarning0(*) 'Time step-wise output is deprecated! Your parameter file is not compatible with GK version!'
#endif
         CASE DEFAULT
            logError(*) 'PickDtType must be 1 = pickdt or 2 = pickdt*dt'
            STOP
         ENDSELECT

       ! energy output on = 1, off =0
       IO%energy_output_on = energy_output_on

       IF(IO%energy_output_on .EQ. 1) THEN
#ifdef GENERATEDKERNELS
            logWarning0(*) 'Energy output currently only working with classic version. Turning it off.'
            IO%energy_output_on = 0
#else
       !own timestep for energy output but output-type is the same as for receivers
       !default type is 1 (time interval-wise)
       !default time interval = 0.1
       IO%pickdt_energy = pickdt_energy
       logInfo0(*) 'current energy dt is', IO%pickdt_energy
#endif
       ENDIF
     IO%printEnergiesInterval = printEnergiesInterval

     IO%nRecordPoint = nRecordPoints  ! number of points to pick temporal signal
     logInfo(*) 'Number of Record Points = ', IO%nRecordPoint
     ALLOCATE( X(IO%nRecordPoint), Y(IO%nRecordPoint), Z(IO%nRecordPoint) )
      ! Read the single record points
      IF (nRecordPoints .GT. 0) THEN
         logInfo(*) 'Record Points read from ', TRIM(RFileName)
         CALL OpenFile(                                 &
               UnitNr       = IO%UNIT%other01         , &
               Name         = RFileName            , &
               create       = .FALSE.                          )

         DO i = 1,nRecordPoints
            READ(IO%UNIT%other01,*) X(i), Y(i), Z(i)

               logInfo0(*) 'in point :'                             !
               logInfo0(*) 'x = ', X(i)         !
               logInfo0(*) 'y = ', Y(i)         !
               logInfo0(*) 'z = ', Z(i)         !

         ENDDO
         !
         CLOSE(IO%UNIT%other01)
      ELSE
         logInfo(*) 'No single record points required. '
      END IF
        ! Allocate the record points for the Unstructured Mesh

      logInfo(*) ' '
      logInfo(*) 'Unstructured record points are allocated'
      logInfo(*) 'Local monitoring time stepping '
      logInfo(*) 'is required in ', IO%nRecordPoint,'points:'
      ! Update total number of record points
      IO%ntotalRecordPoint = IO%nRecordPoint
      logInfo(*) 'Allocating ',IO%ntotalRecordPoint, ' unstructured record points'

      ALLOCATE (                                             &
              IO%UnstructRecPoint(IO%ntotalRecordPoint),     &
              STAT = allocStat                             )

      IF (allocStat .NE. 0) THEN
            logError(*) 'could not allocate',&
                 ' all variables! Ie. Unstructured record Points'
            STOP
      END IF
      !
      IO%UnstructRecPoint(:)%X = X(:)
      IO%UnstructRecPoint(:)%Y = Y(:)
      IO%UnstructRecPoint(:)%Z = Z(:)


      ! Points, where Peak Ground Motion is measured
          IO%PGMLocationsFlag = PGMFlag
          SELECT CASE(IO%PGMLocationsFlag)

            CASE(0)

              logInfo(*) 'No Peak Ground Motion output required ! '

            CASE(1)

              IO%PGMLocationsFile = PGMFile                       !
              logInfo(*) ' '
              logInfo(*) 'Peak Ground Motion (PGM) locations read from file : ', TRIM(IO%PGMLocationsFile)
              CALL OpenFile(                                       &
                   UnitNr       = IO%UNIT%other01                , &
                   Name         = IO%PGMLocationsFile            , &
                   create       = .FALSE.                          )
              READ(IO%UNIT%other01,'(i10)') IO%nPGMRecordPoint                         ! Number of Peak Ground Motion Locations
              logInfo(*) 'Reading ',IO%nPGMRecordPoint,' PGM locations ... '
              logInfo(*) ' '
              ! Update total number of record points
              IO%ntotalRecordPoint = IO%ntotalRecordPoint + IO%nPGMRecordPoint

              ! Enlarge IO%UnstructRecPoint to add PGM record points
              ALLOCATE(IO%tmpRecPoint(IO%nRecordPoint))
              DO i = 1, IO%nRecordPoint
                   IO%tmpRecPoint(i)%X = IO%UnstructRecPoint(i)%X
                   IO%tmpRecPoint(i)%Y = IO%UnstructRecPoint(i)%Y
                   IO%tmpRecPoint(i)%Z = IO%UnstructRecPoint(i)%Z
              ENDDO
              DEALLOCATE(IO%UnstructRecPoint)
              ALLOCATE (IO%UnstructRecPoint(IO%ntotalRecordPoint), &
                   STAT = allocStat                                )
              IF (allocStat.NE.0) THEN
                   logError(*) 'could not allocate all PGM locations in IO%UnstructRecPoint !'
                   STOP
              END IF
              DO i = 1, IO%nRecordPoint
                   IO%UnstructRecPoint(i)%X = IO%tmpRecPoint(i)%X
                   IO%UnstructRecPoint(i)%Y = IO%tmpRecPoint(i)%Y
                   IO%UnstructRecPoint(i)%Z = IO%tmpRecPoint(i)%Z
              ENDDO
              DEALLOCATE(IO%tmpRecPoint)
              DO i = IO%nRecordPoint+1, IO%ntotalRecordPoint
                  READ(IO%UNIT%other01,*)                          &
                       IO%UnstructRecPoint(i)%X,                   &
                       IO%UnstructRecPoint(i)%Y,                   &
                       IO%UnstructRecPoint(i)%Z
              ENDDO

              IO%PGMstartindex = IO%nRecordPoint+1

              CLOSE(IO%UNIT%other01)

            CASE DEFAULT

              logError(*) 'Peak Ground Motion Flag in  O U T P U T  must be set to 0 or 1 ! '
              STOP

          END SELECT
      !
      IF(EQN%DR.NE.0) THEN
          IO%FaultOutputFlag = FaultOutputFlag

        SELECT CASE(IO%FaultOutputFlag)

          CASE(0)

             logInfo(*) 'No Fault Output required ! '

          CASE(1)

             logInfo(*) 'Fault variables will be outputted ! '

          CASE DEFAULT

             logError(*) 'Fault Output Flag in  O U T P U T  must be set to 0 or 1 ! '
             STOP

        END SELECT

      ENDIF

      ! xdmf writer backend
      io%xdmfWriterBackend = xdmfWriterBackend
      select case (io%xdmfWriterBackend)
        case ("posix")
          logInfo0(*) 'Use POSIX XdmfWriter backend'
        case ("hdf5")
#ifndef USE_HDF
          logError(*) 'This version does not support HDF5 backends'
          stop
#endif
          logInfo0(*) 'Use HDF5 XdmfWriter backend'
        case default
          logError(*) 'Unknown XdmfWriter backend ', io%xdmfWriterBackend
          stop
      end select

      ! Check point config
      if (checkPointInterval .lt. 0) then
        logError(*) 'The interval for checkpoints cannot be negative'
        stop
      endif
      io%checkpoint%interval = checkPointInterval
      io%checkpoint%filename = checkPointFile
      io%checkpoint%backend = checkPointBackend

      IO%Refinement = Refinement
      SELECT CASE(Refinement)
         CASE(0)

            logInfo0(*) 'Refinement is disabled'

         CASE(1)

             logInfo0(*) 'Refinement strategy is Face Extraction :  4 subcells per cell'

         CASE(2)

             logInfo0(*) 'Refinement strategy is Equal Face Area : 8 subcells per cell'

         CASE(3)

             logInfo0(*) 'Refinement strategy is Equal Face Area and Face Extraction : 32 subcells per cell'

         CASE DEFAULT

             logError(*) 'Refinement strategy is N O T supported'
             STOP

      END SELECT

      select case (io%checkpoint%backend)
        case ("posix")
            logInfo0(*) 'Using POSIX checkpoint backend'
        case ("hdf5")
#ifndef USE_HDF
            logError(*) 'This version does not support HDF5 checkpoints'
            stop
#endif
            logInfo0(*) 'Using HDF5 checkpoint backend'
        case ("mpio")
#ifndef USE_MPI
            logError(*) 'This version does not support MPI-IO checkpoints'
            stop
#endif
            logInfo0(*) 'Using MPI-IO checkpoint backend'
        case ("mpio_async")
#ifndef USE_MPI
            logError(*) 'This version does not support MPI-IO checkpoints'
            stop
#endif
            logInfo0(*) 'Using async MPI-IO checkpoint backend'
        case ("sionlib")
#ifndef USE_SIONLIB
            logError(*) 'This version does not support SIONlib checkpoints'
            stop
#endif
            logInfo0(*) 'Using SIONlib checkpoint backend'
        case ("none")
            io%checkpoint%interval = 0
        case default
            logInfo0(*) 'Unknown checkpoint backend:', io%checkpoint%backend
            logInfo0(*) 'Disabling checkpoints'
            io%checkpoint%interval = 0
      end select

  END SUBROUTINE readpar_output

  !============================================================================
  ! A B O R T
  !============================================================================

  SUBROUTINE readpar_abort(DISC,IO)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tDiscretization)     :: DISC
    TYPE (tInputOutput)        :: IO
    !------------------------------------------------------------------------
    INTENT(INOUT)              :: DISC,IO
    !------------------------------------------------------------------------
    INTEGER                          :: MaxIteration
    INTEGER                          :: readStat
    REAL                             :: EndTime, MaxTolerance, MaxTolCriterion, WallTime_h, Delay_h
    NAMELIST                         /AbortCriteria/ EndTime, MaxIteration, MaxTolerance, &
                                                      MaxTolCriterion, WallTime_h, Delay_h
    !------------------------------------------------------------------------
    !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<  A B O R T - C R I T E R I A                            >'
    logInfo(*) '<--------------------------------------------------------->'

    ! Setting default values
    EndTime = 15.0
    MaxIteration = 10000000
    MaxTolerance = 1.0E-07
    MaxTolCriterion = 0
    WallTime_h = 1e20
    Delay_h = 0.

   READ(IO%UNIT%FileIn, IOSTAT=readStat, nml = AbortCriteria)
    IF (readStat.NE.0) THEN
        CALL RaiseErrorNml(IO%UNIT%FileIn, "AbortCriteria")
    ENDIF

    DISC%EndTime =  EndTime                                         ! time required

    logInfo(*) 'Maximum computed TIME allowed:',    DISC%EndTime
    !
    DISC%MaxIteration =  MaxIteration                               ! nr of the end iteration
    !
#ifdef GENERATEDKERNELS
    if (DISC%MaxIteration .lt. 10000000) then
      logError(*) 'GK version does not support MaxIteration!'
      stop
    endif
#else
    if (DISC%MaxIteration .lt. 10000000) then
      logWarning(*) 'MaxIteration is deprecated! Your parameter file is not compatible with GK version!'
    endif
#endif
    logInfo(*) 'Maximum ITERATION number allowed:', DISC%MaxIteration
    !
    !
        IO%WallTime_h = WallTime_h
        IO%Delay_h = Delay_h

        IO%WallTime_s = 3600.*IO%WallTime_h
        IO%Delay_s    = 3600.*IO%Delay_h
    !
  END SUBROUTINE readpar_abort

  !============================================================================
  ! A N A L Y S E
  !============================================================================
  ! Checks the correct setting of the .par file
  SUBROUTINE analyse_readpar(EQN,DISC,MESH,IC,SOURCE,IO,MPI)
  !SUBROUTINE analyse_readpar_unstruct(EQN,DISC,MESH,IC,SOURCE,IO,MPI)
    !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------------
    TYPE (tEquations)          :: EQN
    TYPE (tDiscretization)     :: DISC
    TYPE (tUnstructMesh)       :: MESH
    TYPE (tInitialCondition)   :: IC
    TYPE (tSource)             :: SOURCE
    TYPE (tInputOutput)        :: IO
    TYPE (tMPI), OPTIONAL      :: MPI

    ! local variables
    CHARACTER(LEN=80)          :: Filerestart
    CHARACTER(LEN=256)         :: e(1000)
    INTEGER                    :: i
    INTEGER                    :: status
    !------------------------------------------------------------------------
    INTENT(IN)                 :: EQN, DISC, MESH, SOURCE, IO
    !------------------------------------------------------------------------

! Generated kernels sanity check
#ifdef GENERATEDKERNELS
    if (NUMBER_OF_QUANTITIES .NE. EQN%nVarTotal) then
      logError(*) 'Generated kernels: The number of quantities defined by the parameter file (', EQN%nVarTotal, ') does not the number of quantities this version was compiled for (', NUMBER_OF_QUANTITIES, ').'
      stop
    end if
#endif

    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '<  END OF PARAMETER FILE                                  >'
    logInfo(*) '<-------------------------------------------------------- >'
    logInfo(*) ' '
    !
  END SUBROUTINE analyse_readpar
  !END SUBROUTINE analyse_readpar_unstruct

END MODULE COMMON_readpar_mod
