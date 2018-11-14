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

#include <Initializer/preProcessorMacros.fpp>

MODULE data_output_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  USE COMMON_operators_mod
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE data_output
     MODULE PROCEDURE data_output
  END INTERFACE

  INTERFACE GalerkinFineOutput
     MODULE PROCEDURE GalerkinFineOutput
  END INTERFACE

  INTERFACE GalerkinFineOutput_Adj
     MODULE PROCEDURE GalerkinFineOutput_Adj
  END INTERFACE
  INTERFACE KernelsOutput
     MODULE PROCEDURE KernelsOutput
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC :: data_output, &
            GalerkinFineOutput_Adj, &
            KernelsOutput
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE data_output(dt,time,timestep,EQN,MESH,DISC, &
       SOURCE,BND,IO,MPI,OptionalFields,ANALYSE)
    !--------------------------------------------------------------------------
    
    USE plot_fields_mod        , ONLY : plot_fields
    use WaveFieldWriter
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tEquations)             :: EQN
    TYPE (tUnstructMesh)          :: MESH 
    TYPE (tDiscretization)        :: DISC
    TYPE (tSource)                :: SOURCE
    TYPE (tBoundary)              :: BND
    TYPE (tInputOutput)           :: IO
    TYPE (tMPI)                   :: MPI
    TYPE (tUnstructOptionalFields):: OptionalFields
    TYPE (tAnalyse)               :: ANALYSE
    REAL                          :: dt(    MESH%nElem           )
    INTEGER                       :: timestep
    REAL                          :: time
    ! local Variables
    CHARACTER (LEN=350)           :: outfile, outfile_t, outfile_q ! auxiliary output mesh
    CHARACTER (LEN=350)           :: filename                      ! auxiliary output mesh
    CHARACTER (LEN=10)            :: citer                         ! CHAR. for plot file 
    CHARACTER (LEN=350)           :: DataFilename
    CHARACTER (LEN=5)             :: cmyrank
    !--------------------------------------------------------------------------
    INTENT(IN)                    :: dt,time,timestep
    INTENT(IN)                    :: EQN,DISC,SOURCE,BND                                    
    INTENT(INOUT)                 :: IO, OptionalFields, MESH

    ! register epik/scorep function data_output
    EPIK_FUNC_REG("data_output")
    SCOREP_USER_FUNC_DEFINE()
    !--------------------------------------------------------------------------
    !                                                                          !
    ! start epik function data_output
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN("data_output")

    if (IO%Format .eq. 6) then
        !call WaveFieldWriterWriteStep(time, disc, mesh, mpi)
    else
        if (IO%Format .eq. 10) then
            ! Output is disabled

            ! end epik/scorep function data_output
            EPIK_FUNC_END()
            SCOREP_USER_FUNC_END()
            return
        endif                                                                    !
        !                                                                          !
        IF (time.LE.DISC%EndTime) THEN                                             ! time and time step
           !                                                                       ! information in
           !                                                                       ! calculation
           logInfo(*)'-----------------------------------------------------------'
           !                                                                       !
           logInfo(*)'    PLOT wavefield for current Time   :',time                ! Unsteady calc
           logInfo(*)'                        iteration #   :',timestep            ! Unsteady calc
           logInfo(*)'                  current time step   :',dt(1)               ! Unsteady calc
           !                                                                       !
        ELSE                                                                       ! last time step: only
           !                                                                       ! time information, no
           !                                                                       ! further time steps
           !                                                                       ! are calculated
           logInfo(*)'-----------------------------------------------------------'
           !                                                                       !
           logInfo(*)'    PLOT wavefield for current Time   :',time                ! Unsteady calc
           logInfo(*)'                        iteration #   :', timestep           ! Unsteady calc
           !                                                                       !
        ENDIF                                                                      !
        !                                                                          !
        WRITE(citer,'(I9.9)') timestep                                             ! timestep -> citer
        !                                                                          !
#ifdef PARALLEL
        WRITE(cmyrank,'(I5.5)') MPI%myrank                                         ! myrank -> cmyrank
        if (IO%Format .eq. 5) then
            ! HDF5
            outfile = TRIM(IO%OutputFile )//'-'//TRIM(citer)//'.'//TRIM(cmyrank)
        else
            outfile = TRIM(IO%OutputFile )//'-'//TRIM(cmyrank)//'/timestep-'//TRIM(citer)     ! mshfile Name
        endif
          CALL plot_fields(                                                         & !
                fileop     = outfile                                              , & !
                time       = time                                                 , & !
                timestep   = timestep                                             , & !
                EQN        = EQN                                                  , & !
                DISC       = DISC                                                 , & !
                MESH       = MESH                                                 , & !
                IO         = IO                                                   , & !
                BND        = BND                                                  , & !
                OptionalFields = OptionalFields                                   , & !
                MPI        = MPI                                                    )

#else
          outfile = TRIM(IO%OutputFile )//'-'//TRIM(citer)                      ! mshfile Name
          CALL plot_fields(                                                         & !
                fileop     = outfile                                              , & !
                time       = time                                                 , & !
                timestep   = timestep                                             , & !
                EQN        = EQN                                                  , & !
                DISC       = DISC                                                 , & !
                MESH       = MESH                                                 , & !
                IO         = IO                                                   , & !
                BND        = BND                                                  , & !
                OptionalFields = OptionalFields                                   , & !
                MPI        = MPI                                                    )
#endif
        !                                                                          !
        IF(DISC%Galerkin%DGFineOut1D.GT.0) THEN
           CALL GalerkinFineOutput(time,timestep,EQN,MESH,DISC,IO,MPI)
        ENDIF
    endif
    !                                                                          !
    ! end epik/scorep function data_output
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()
  END SUBROUTINE data_output                                                   !

  SUBROUTINE GalerkinFineOutput(time,timestep,EQN,MESH,DISC,IO,MPI)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tEquations)             :: EQN 
    TYPE (tUnstructMesh)          :: MESH 
    TYPE (tDiscretization)        :: DISC
    TYPE (tInputOutput)           :: IO
    TYPE (tMPI)                   :: MPI
    REAL                          :: time
    INTEGER                       :: timestep
    ! local Variables
    INTEGER                       :: i,ielem,iNode, iVar, iDegFr
    INTEGER                       :: iPoint
    INTEGER                       :: stat
    REAL                          :: xGP,yGP,xi,eta
    REAL                          :: state(EQN%nVar)
    CHARACTER (LEN=350)           :: Filename
    CHARACTER (LEN=10)            :: citer            ! CHAR. for plot file 
    CHARACTER (LEN=5)             :: cmyrank
    !--------------------------------------------------------------------------
    INTENT(IN)                    :: EQN,MESH,DISC,IO,time,timestep
    !--------------------------------------------------------------------------
    !
    WRITE(citer,'(I9.9)') timestep                                             ! timestep -> citer
#ifdef PARALLEL
    WRITE(cmyrank,'(I5.5)') MPI%myrank                                         ! myrank -> cmyrank
    Filename = TRIM(IO%OutputFile )//'.GF.'//TRIM(citer)//'.'//TRIM(cmyrank)//'.'//'dat'
#else
    Filename = TRIM(IO%OutputFile )//'.GF.'//TRIM(citer)//'.dat'
#endif
    !
    OPEN(UNIT   = IO%UNIT%other01                   , &                        ! open file
         FILE   = TRIM(Filename)                  ,   &                        ! open file
         STATUS = 'UNKNOWN'                         , &                        ! open file
         RECL   = 700                               , &                        ! open file
         IOSTAT = stat                                )                        ! open file
    IF (stat.NE.0) THEN                                                        ! Error Handler
       logError(*) 'cannot open file in GalerkinFineOutput!'                   ! Error Handler
       STOP                                                                    ! Error Handler
    END IF                                                                     ! Error Handler
    !
    ! Write the header                                                         
    WRITE(IO%UNIT%other01,*) '# Galerkin Fine Output 4.0 ' 
    ! Write polynomial order and the number of degrees of freedom per element
    WRITE(IO%UNIT%other01,*) '# Version ID '
    WRITE(IO%UNIT%other01,*) 5
    WRITE(IO%UNIT%other01,*) '# Number of dimensions '
    WRITE(IO%UNIT%other01,*) EQN%Dimension
    WRITE(IO%UNIT%other01,*) '# Time / Timestep '
    WRITE(IO%UNIT%other01,*) time, timestep   
    WRITE(IO%UNIT%other01,*) '# Number of variables '
    WRITE(IO%UNIT%other01,*) EQN%nVar
    WRITE(IO%UNIT%other01,*) '# Order of basis / Number of DOF per element'
    WRITE(IO%UNIT%other01,*) DISC%Galerkin%nPolyRec, DISC%Galerkin%nDegFrRec
    WRITE(IO%UNIT%other01,*) '# Type of elements / Number of elements / Number of nodes '
    WRITE(IO%UNIT%other01,*) MESH%GlobalElemType, MESH%nElem, MESH%nNode
    WRITE(IO%UNIT%other01,*) '# The nodes of the mesh: '                
    DO iNode = 1, MESH%nNode
       WRITE(IO%UNIT%other01,*) MESH%VRTX%xyNode(:,iNode)
    ENDDO
    WRITE(IO%UNIT%other01,*) '# The connectivity of the mesh: '                
    DO iElem = 1, MESH%nElem
       WRITE(IO%UNIT%other01,*) MESH%ELEM%Vertex(:,iElem)
    ENDDO
    WRITE(IO%UNIT%other01,*) '# The boundary conditions of the mesh: '                
    DO iElem = 1, MESH%nElem
       WRITE(IO%UNIT%other01,*) MESH%ELEM%Reference(1:MESH%nVertexMax,iElem)
    ENDDO
    WRITE(IO%UNIT%other01,*) '# The distribution of the polynomial degrees '                
    IF(DISC%Galerkin%pAdaptivity.GT.0) THEN
      DO iElem = 1, MESH%nElem
         WRITE(IO%UNIT%other01,*) DISC%Galerkin%LocPoly(iElem)
      ENDDO
    ELSE
      DO iElem = 1, MESH%nElem
         WRITE(IO%UNIT%other01,*) DISC%Galerkin%nPolyRec
      ENDDO
    ENDIF
    WRITE(IO%UNIT%other01,*) '# The degrees of freedom: '                
    DO ielem = 1,MESH%nElem                                                  
        DO iDegFr = 1, DISC%Galerkin%nDegFr                                    
           WRITE(IO%UNIT%other01,*) DISC%Galerkin%dgvar(iDegFr,1:EQN%nVar,iElem,1)
        ENDDO                                                                  
    ENDDO                                                    
    !
    CLOSE(IO%UNIT%other01)
    !
  END SUBROUTINE GalerkinFineOutput

  SUBROUTINE GalerkinFineOutput_Adj(time,timestep,EQN,MESH,DISC,IO,MPI)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tEquations)             :: EQN 
    TYPE (tUnstructMesh)          :: MESH 
    TYPE (tDiscretization)        :: DISC
    TYPE (tInputOutput)           :: IO
    TYPE (tMPI)                   :: MPI
    REAL                          :: time
    INTEGER                       :: timestep
    ! local Variables
    INTEGER                       :: i,ielem,iNode, iVar, iDegFr
    INTEGER                       :: iPoint
    INTEGER                       :: stat
    REAL                          :: xGP,yGP,xi,eta
    REAL                          :: state(EQN%nVar)
    CHARACTER (LEN=350)           :: Filename
    CHARACTER (LEN=10)            :: citer            ! CHAR. for plot file 
    CHARACTER (LEN=5)             :: cmyrank
    !--------------------------------------------------------------------------
    INTENT(IN)                    :: EQN,MESH,DISC,IO,time,timestep
    !--------------------------------------------------------------------------
    !
    WRITE(citer,'(I9.9)') timestep                                             ! timestep -> citer
#ifdef PARALLEL
    WRITE(cmyrank,'(I5.5)') MPI%myrank                                         ! myrank -> cmyrank
    Filename = TRIM(IO%OutputFile )//'_adj.GF.'//TRIM(citer)//'.'//TRIM(cmyrank)//'.'//'dat'
#else
    Filename = TRIM(IO%OutputFile )//'_adj.GF.'//TRIM(citer)//'.dat'
#endif
    !
    OPEN(UNIT   = IO%UNIT%other01                   , &                        ! open file
         FILE   = TRIM(Filename)                  ,   &                        ! open file
         STATUS = 'UNKNOWN'                         , &                        ! open file
         RECL   = 700                               , &                        ! open file
         IOSTAT = stat                                )                        ! open file
    IF (stat.NE.0) THEN                                                        ! Error Handler
       logError(*) 'cannot open file in GalerkinFineOutput!'                   ! Error Handler
       STOP                                                                    ! Error Handler
    END IF                                                                     ! Error Handler
    !
    ! Write the header                                                         
    WRITE(IO%UNIT%other01,*) '# Galerkin Fine Output 4.0 '

    ! Write polynomial order and the number of degrees of freedom per element
    WRITE(IO%UNIT%other01,*) '# Version ID '
    WRITE(IO%UNIT%other01,*) 5
    WRITE(IO%UNIT%other01,*) '# Number of dimensions '
    WRITE(IO%UNIT%other01,*) EQN%Dimension
    WRITE(IO%UNIT%other01,*) '# Time / Timestep '
    WRITE(IO%UNIT%other01,*) time, timestep   
    WRITE(IO%UNIT%other01,*) '# Number of variables '
    WRITE(IO%UNIT%other01,*) EQN%nVar+2          ! We additionally output the kernels computed
    WRITE(IO%UNIT%other01,*) '# Order of basis / Number of DOF per element'
    WRITE(IO%UNIT%other01,*) DISC%Galerkin%nPolyRec, DISC%Galerkin%nDegFrRec
    WRITE(IO%UNIT%other01,*) '# Type of elements / Number of elements / Number of nodes '
    WRITE(IO%UNIT%other01,*) MESH%GlobalElemType, MESH%nElem, MESH%nNode
    WRITE(IO%UNIT%other01,*) '# The nodes of the mesh: '                
    DO iNode = 1, MESH%nNode
       WRITE(IO%UNIT%other01,*) MESH%VRTX%xyNode(:,iNode)
    ENDDO
    WRITE(IO%UNIT%other01,*) '# The connectivity of the mesh: '                
    DO iElem = 1, MESH%nElem
       WRITE(IO%UNIT%other01,*) MESH%ELEM%Vertex(:,ielem)
    ENDDO
    WRITE(IO%UNIT%other01,*) '# The boundary conditions of the mesh: '                
    DO iElem = 1, MESH%nElem
       WRITE(IO%UNIT%other01,*) MESH%ELEM%Reference(1:MESH%nVertexMax,iElem)
    ENDDO

    WRITE(IO%UNIT%other01,*) '# The distribution of the polynomial degrees '                
    IF(DISC%Galerkin%pAdaptivity.GT.0) THEN
      DO ielem = 1, MESH%nElem
         WRITE(IO%UNIT%other01,*) DISC%Galerkin%LocPoly(iElem)
      ENDDO
    ELSE
      DO ielem = 1, MESH%nElem
         WRITE(IO%UNIT%other01,*) DISC%Galerkin%nPolyRec
      ENDDO
    ENDIF

    WRITE(IO%UNIT%other01,*) '# The degrees of freedom: '                
        DO ielem = 1,MESH%nElem                                                  
         DO iDegFr = 1, DISC%Galerkin%nDegFr                                    
           WRITE(IO%UNIT%other01,*) DISC%Adjoint%dgvar_adj(iDegFr,1:EQN%nVar,ielem,1),DISC%Adjoint%kernel(iDegFr,ielem,1),DISC%Adjoint%kernel(iDegFr,ielem,2)
         ENDDO                                                                  
        ENDDO                    
    !
    CLOSE(IO%UNIT%other01)
    !
  END SUBROUTINE GalerkinFineOutput_Adj
 
  SUBROUTINE KernelsOutput(EQN,MESH,DISC,IO,MPI)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tEquations)             :: EQN 
    TYPE (tUnstructMesh)          :: MESH 
    TYPE (tDiscretization)        :: DISC
    TYPE (tInputOutput)           :: IO
    TYPE (tMPI)                   :: MPI
    ! local Variables
    INTEGER                       :: i,ielem,iNode, iVar, iDegFr
    CHARACTER (LEN=350)           :: Filename
    CHARACTER (LEN=10)            :: citer            ! CHAR. for plot file 
    CHARACTER (LEN=5)             :: cmyrank
    !--------------------------------------------------------------------------
    INTENT(IN)                    :: EQN,MESH,DISC,IO
    !--------------------------------------------------------------------------  

#ifdef PARALLEL
    WRITE(cmyrank,'(I5.5)') MPI%myrank                                         ! myrank -> cmyrank
    Filename = TRIM(IO%OutputFile )//'_kernels.'//TRIM(cmyrank)//'.'//'dat'
#else
    Filename = TRIM(IO%OutputFile )//'_kernels.dat'
#endif
    !
    ! The extra 2 value accounts for the double precision
    OPEN(         UNIT=IO%UNIT%other01,                &
                  FILE=TRIM(Filename),                 &
                  FORM='unformatted',                  &
                  ACCESS='direct',                     &
                  RECL=2*DISC%Galerkin%nDegFr*MESH%nElem*DISC%Adjoint%nKernels, &
                  STATUS='unknown'                    )
    WRITE(IO%UNIT%other01, REC=1) DISC%Adjoint%Kernel
    CLOSE(IO%UNIT%other01)   

  END SUBROUTINE KernelsOutput
                                                                                                     !
END MODULE data_output_mod
