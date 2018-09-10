!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2014-2016, SeisSol Group
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
!! routine outputs the rupture front arrival time for each GP node for an iFace

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE output_rupturefront_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  PUBLIC  :: output_rupturefront
  !---------------------------------------------------------------------------!  
  INTERFACE output_rupturefront
     MODULE PROCEDURE output_rupturefront
  END INTERFACE

CONTAINS

    SUBROUTINE output_rupturefront(DISC,MESH,MPI,IO, BND)
    !< routine outputs the rupture front arrival time for each GP node for an iFace
    !-------------------------------------------------------------------------!
    USE DGBasis_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization), target   :: DISC                                              !
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    TYPE (tBoundary)                :: BND
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: iFace, iBndGP, iElem, iSide
    INTEGER                         :: stat, UNIT_RF, lines
    INTEGER                         :: iLocalNeighborSide,iNeighbor
    INTEGER                         :: MPIIndex, iObject
#ifdef OMP
    INTEGER                         :: TID,omp_get_thread_num
    CHARACTER (LEN=2)               :: c_TID
#endif
    REAL                            :: xV(1:4),yV(1:4),zV(1:4)
    REAL                            :: chi,tau,xi,eta,zeta
    REAL                            :: xGP,yGP,zGP
    LOGICAL                         :: exist
    CHARACTER (LEN=5)               :: cmyrank
    CHARACTER (len=200)             :: RF_FILE
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, MESH, MPI, IO
    !-------------------------------------------------------------------------!
    IF ( ALL(DISC%DynRup%rupture_time .EQ. 0.0 ) ) THEN
        RETURN !no output if no element ruptured
    ENDIF

    ! generate unique name out of iElem,iBndGP and MPI rank
#ifdef PARALLEL
#ifdef OMP
    ! hybrid MPI + openMP case
    TID = omp_get_thread_num()
    WRITE(c_TID,'(I2.2)') TID
    WRITE(cmyrank,'(I5.5)') MPI%myrank                           ! myrank -> cmyrank
    WRITE(RF_FILE, '(a,a4,a5,a5,a2,a4)') TRIM(IO%OutputFile),'-RF-',TRIM(cmyrank),'-TID-',TRIM(c_TID),'.dat'
    UNIT_RF = 99875+MPI%myrank+TID
#else
    ! pure MPI case
    WRITE(cmyrank,'(I5.5)') MPI%myrank                           ! myrank -> cmyrank
    WRITE(RF_FILE, '(a,a4,a5,a4)') TRIM(IO%OutputFile),'-RF-',TRIM(cmyrank),'.dat'
    UNIT_RF = 99875+MPI%myrank
#endif
#else
    WRITE(RF_FILE, '(a,a4,a4)') TRIM(IO%OutputFile),'-RF-','.dat'
    UNIT_RF = 99875
#endif
    !
    INQUIRE(FILE = RF_FILE, EXIST = exist)
    IF(exist) THEN
        ! If file exists, then append data
        OPEN(UNIT     = UNIT_RF                                          , & !
             FILE     = RF_FILE                                          , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'OLD'                                            , & !
             POSITION = 'APPEND'                                         , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        IF(stat.NE.0) THEN                                              !
           logError(*) 'cannot open ',RF_FILE         !
           logError(*) 'Error status: ', stat                !
           STOP                                                          !
        ENDIF
    ELSE
        ! open file
        OPEN(UNIT   = UNIT_RF                                            , & !
             FILE     = RF_FILE                                          , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'NEW'                                            , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        IF(stat.NE.0) THEN                                              !
           logError(*) 'cannot open ',RF_FILE         !
           logError(*) 'Error status: ', stat                !
           STOP                                                          !
        ENDIF
        !
        ! add header with information of total nr of lines:
        lines = MESH%Fault%nSide * DISC%Galerkin%nBndGP
        WRITE(UNIT_RF,*) lines
    ENDIF

    ! !loop over the positive fault sides, for each GP
    DO iFace=1,mesh%fault%nSide

    iElem               = MESH%Fault%Face(iFace,1,1)          ! Remark:
    iSide               = MESH%Fault%Face(iFace,2,1)          ! iElem denotes "+" side
    DO iBndGP=1,DISC%Galerkin%nBndGP

        IF (MESH%Fault%Face(iFace,1,1) == 0) THEN
          ! iElem is in the neighbor domain
          ! The neighbor element belongs to a different MPI domain
            iNeighbor           = MESH%Fault%Face(iFace,1,2)          ! iNeighbor denotes "-" side
            iLocalNeighborSide  = MESH%Fault%Face(iFace,2,2)
            iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
            MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)

            xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
            yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
            zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
        ELSE
           ! get vertices
            xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
            yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
            zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
        ENDIF

    ! get vertices of complete tet
    !xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
    !yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
    !zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
    !
    ! Transformation of boundary GP's into XYZ coordinate system
    chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
    tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
    CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
    CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)


    ! Write output for positive elements
    IF (MESH%FAULT%Face(iFace,1,1).NE.0) THEN
       IF (DISC%DynRup%rupture_time(iBndGP,iFace).GT. 0.0) THEN
           IF (DISC%DynRup%DS_output_on .EQ. 1) THEN
               WRITE(UNIT_RF,*) xGP, yGP, zGP, DISC%DynRup%rupture_time(iBndGP,iFace), DISC%DynRup%dynStress_time(iBndGP,iFace)
           ELSE
               WRITE(UNIT_RF,*) xGP, yGP, zGP, DISC%DynRup%rupture_time(iBndGP,iFace)
           ENDIF
       ENDIF
    ENDIF

    ENDDO
    ENDDO
    CLOSE( UNIT_RF )

  END SUBROUTINE output_rupturefront

END MODULE output_rupturefront_mod
