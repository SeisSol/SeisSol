!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Cristobal E. Castro (ccastro AT uc.cl https://sites.google.com/a/uc.cl/cristobal-e-castro/)
!!
!! @section LICENSE
!! Copyright (c) 2008, SeisSol Group
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

MODULE ini_MODEL_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE ini_MODEL
     MODULE PROCEDURE ini_MODEL
  END INTERFACE
  INTERFACE ini_ATTENUATION
     MODULE PROCEDURE ini_ATTENUATION
  END INTERFACE
  INTERFACE ModelDefinition_new
     MODULE PROCEDURE ModelDefinition_new
  END INTERFACE
  INTERFACE generate_FacetList
     MODULE PROCEDURE generate_FacetList
  END INTERFACE
  INTERFACE calc_BaryDist
     MODULE PROCEDURE calc_BaryDist
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: ini_MODEL,                &
             ini_ATTENUATION,          &
             generate_FacetList,       &
             ModelDefinition_new,      &
             SolveRiemannProblem
  !----------------------------------------------------------------------------
  integer,  parameter, private :: rk=kind(1.)
  real(rk), parameter, private :: fx_large=1.d30

CONTAINS


  SUBROUTINE ini_MODEL(MaterialVal,OptionalFields,EQN,MESH,IC,IO,DISC,BND)
    !--------------------------------------------------------------------------

    USE COMMON_operators_mod, ONLY: OpenFile, XYinTriangle
    USE QuadPoints_mod
    USE DGBasis_mod
    USE TrilinearInterpolation_mod
    USE read_backgroundstress_mod
    USE ini_model_DR_mod
    use VelocityFieldReader
    use modules

    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tEquations)               :: EQN
    TYPE (tUnstructMesh)            :: MESH
    TYPE (tInputOutput)             :: IO
    TYPE (tInitialCondition)        :: IC
    TYPE (tDiscretization)          :: DISC
    TYPE (tBoundary)                :: BND
    TYPE (tUnstructOptionalFields)  :: OptionalFields
    REAL                            :: MaterialVal(MESH%nElem,EQN%nBackgroundVar)
    ! Local variable declaration
    LOGICAL                         :: FoundLayer, DIVZERO=.FALSE., nodewise=.FALSE.
    CHARACTER(LEN=600)              :: Name
    INTEGER                         :: i, j, k, ielem, iVertex, iLayer, iMech, iRFFlag, iMaterial, nVar, mu_ct, la_ct
    INTEGER                         :: iDegFr, iIntGP
    INTEGER                         :: NX, NY, NZ, zoneNum, nPert, counter
    INTEGER, POINTER                :: PertMaterial(:)
    INTEGER, POINTER                :: ElemID(:)
    REAL                            :: Mat_Vert(MESH%nVertexMax,EQN%nBackgroundVar), Mat_Sum(EQN%nBackgroundVar), MaterialTmp(EQN%nAneMaterialVar)
    REAL                            :: x, y, z, r, rBary, cp, cs, rho, mu, lambda, depths(4), radius, radius_ref
    REAL, POINTER                   :: xrf(:), yrf(:), zrf(:), pertrf(:,:,:,:)
    REAL, POINTER                   :: PerturbationVar(:,:), P(:,:)
    REAL, POINTER                   :: posx(:), posy(:), posz(:), pert(:)
    REAL                            :: Material_INF(EQN%nAneMaterialVar-3)
    REAL                            :: xrf_max,xrf_min,yrf_max,yrf_min,zrf_max,zrf_min
    REAL                            :: posx_max,posx_min,posy_max,posy_min,posz_max,posz_min
    REAL                            :: pert_max, pert_min
    REAL, allocatable, dimension(:) :: BaryDist
    REAL                            :: omega
    REAL                            :: circ
    REAL                            :: QPLocVal,QSLocVal
    REAL                            :: Theta(EQN%nMechanisms,3)
    REAL                            :: w_freq(EQN%nMechanisms)
    REAL                            :: material(1:3), MaterialVal_k
    REAL                            :: ZoneIns, ZoneTrans, xG, yG, X2, LocX(3), LocY(3), tmp
    REAL                            :: BedrockVelModel(10,4)
    INTEGER                         :: eType
    INTEGER         :: nTens3GP
    REAL,POINTER    :: Tens3GaussP(:,:)
    REAL,POINTER    :: Tens3GaussW(:)
    REAL,POINTER    :: Tens3BaseFunc(:,:)

    ! ------------------------------------------------------------------------!

    INTEGER                         :: nGaussPoints_quad                      ! Number of Gauss points for integral computation on quads
    REAL, POINTER                   :: GaussPoint_quad(:,:) => NULL()         ! Points
    REAL, POINTER                   :: GaussWeight_quad(:)  => NULL()         ! Weights
    INTEGER                         :: nGaussPoints_tri                       ! Number of Gauss points for integral computation on tri
    REAL, POINTER                   :: GaussPoint_tri(:,:) => NULL()          ! Points
    REAL, POINTER                   :: GaussWeight_tri(:)  => NULL()          ! Weights
    INTEGER                         :: nGaussPoints                           ! Number of Gauss points for integral computation
    REAL, POINTER                   :: GaussPoint(:,:) => NULL()              ! Pointer to the corresponding Points
    REAL, POINTER                   :: GaussWeight(:)  => NULL()              ! Pointer to the corresponding Weights

    REAL, POINTER                   :: BaseFuncVal_quad(:,:) => NULL()        ! Basis functions for quads at Gauss points (iPhi,iIntGp)
    REAL, POINTER                   :: BaseFuncVal_tri(:,:)  => NULL()        ! Basis functions for tria at Gauss points (iPhi,iIntGp)
    REAL, POINTER                   :: BaseFuncVal(:,:) => NULL()             ! Pointer to the corresponding set of basis functions

    REAL, POINTER                   :: iMassMatrix(:,:) => NULL()             ! Pointer to the corresponding mass matrix
    INTEGER                         :: LocElemType                            ! Type of element
    REAL                            :: Pf, b13, b33, b11                      ! fluid pressure and coeffcients for special initial loading in TPV26/TPV27
    INTEGER                         :: InterpolationScheme = 1                ! Select the interpolation scheme (linear=1, cubic=else)
    !--------------------------------------------------------------------------
    INTENT(IN)                      :: MESH,IC
    INTENT(OUT)                     :: MaterialVal
    INTENT(INOUT)                   :: EQN,OptionalFields
    ! -------------------------------------------------------------------------

    IF(EQN%nBackgroundVar.LE.0) THEN
       logInfo(*) 'No values specified. Exiting ini_MODEL.f90. '
       RETURN
    ENDIF

    ! Call the pre model hooks
    call call_hook_pre_model()

    ALLOCATE ( EQN%LocAnelastic(MESH%nElem) )                 ! Don't use anelasticity by default
    EQN%LocAnelastic(:) = 0
    ALLOCATE ( EQN%LocPoroelastic(MESH%nElem) )               ! Don't use poroelasticity by default
    EQN%LocPoroelastic(:) = 0

    IF (EQN%Plasticity .NE. 0) THEN
        ALLOCATE ( EQN%PlastCo(MESH%nElem) )
        EQN%PlastCo(:) = EQN%PlastCo_0 !assign constant value from parameter file
        !add element-dependent assignement for special lintypes in the following
        allocate(EQN%IniStress(6,MESH%nElem))
        EQN%IniStress(:,:) = 0.0D0
    ENDIF

      SELECT CASE(EQN%LinType)
      CASE(0)
        MaterialVal(:,1) = EQN%rho0
        MaterialVal(:,2) = EQN%mu
        MaterialVal(:,3) = EQN%lambda
        IF(EQN%Advection.EQ.1) THEN
            MaterialVal(:,4) = EQN%u0
            MaterialVal(:,5) = EQN%v0
            MaterialVal(:,6) = EQN%w0
        ENDIF
      CASE(1)
        IF(EQN%Anelasticity.EQ.0.AND.EQN%Poroelasticity.EQ.0)THEN
           DO iElem = 1, MESH%nElem
              iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
              MaterialVal(iElem,:) = EQN%MODEL(iLayer,:)   ! Set material properties for this zone.
           ENDDO

        ELSEIF(EQN%Anelasticity.EQ.1.AND.EQN%Poroelasticity.EQ.0)THEN
           DO iElem = 1, MESH%nElem
              iLayer = MESH%ELEM%Reference(0,iElem)                                  ! Zone number is given by reference 0
              MaterialTmp(:) = EQN%MODEL(iLayer,:)                                   ! Get material properties (temporary)
              MaterialVal(iElem,1:3) = EQN%MODEL(iLayer,1:3)                         ! Set elastic properties for this zone.
              !
              ! Check for anelastic material by evaluating Qp and Qs of the parameter file.
              IF(     (MaterialTmp(EQN%nAneMaterialVar-1).LT.9999.) &                 !
                  .OR.(MaterialTmp(EQN%nAneMaterialVar)  .LT.9999.) )THEN             !
                                                                                     !
                  EQN%LocAnelastic(iElem) = 1                                        ! Mark element with anelastic material
                                                                                     !
                  CALL ini_ATTENUATION(Theta,w_freq,Material_INF,MaterialTmp,EQN)    ! Initialize anelastic coefficients for this zone
                  !
                  MaterialVal(iElem,2:EQN%AneMatIni-1) = Material_INF(:)             ! Set unrelaxed material properties for this zone.                                                                      !
                  ! Fill MaterialVal vector for each element with anelastic coefficients w_freq and theta
                  DO iMech = 1, EQN%nMechanisms
                     MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1))             = w_freq(iMech)
                     MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1)+1:EQN%AneMatIni+4*(iMech-1)+3) = Theta(iMech,:)
                  ENDDO
              ELSE
                  MaterialVal(iElem,EQN%AneMatIni:EQN%nBackgroundVar) = 0.                       ! otherwise set all anelastic coefficients to zero
              ENDIF
           ENDDO
        ENDIF
        ! Poroelastic case
        IF(EQN%Anelasticity.EQ.0.AND.EQN%Poroelasticity.NE.0) THEN
           DO iElem = 1, MESH%nElem
              iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
              MaterialVal(iElem,:) = EQN%MODEL(iLayer,:)   ! Set material properties for this zone.
              IF(MaterialVal(iElem,27).GE.1.0e-8) THEN
                  EQN%LocPoroelastic(iElem) = EQN%Poroelasticity
              ELSE
                  MaterialVal(iElem,27) = 0.
              ENDIF
           ENDDO
        ENDIF
      CASE(2,11)     !special case for radially symmetric PREM data
        IF(EQN%Anelasticity.EQ.0)THEN
           nVar = 3
        ELSE
           nVar = 5
        ENDIF
        DO iElem=1,MESH%nElem
           x = MESH%ELEM%xyBary(1,iElem)
           y = MESH%ELEM%xyBary(2,iElem)
           z = MESH%ELEM%xyBary(3,iElem)
           rBary = SQRT(x*x + y*y + z*z)
           rBary = FLOAT(NINT(rBary * 1.0)) / 1.0
           mu_ct = 0
           la_ct = 0
           DO iVertex=1,MESH%nVertexMax
              x = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(iVertex,iElem))
              y = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(iVertex,iElem))
              z = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(iVertex,iElem))
              r = SQRT(x*x + y*y + z*z)
              r = FLOAT(NINT(r * 1.0)) / 1.0 !avoid roundoff error
              SELECT CASE(EQN%LinType)
                 CASE(2)
                    CALL prem_iso(r,rBary,nVar,Mat_Vert(iVertex,1:nVar),EQN%Anelasticity)
                 CASE(11)
                    DO iLayer = 1, EQN%nLayers
                       IF( (r.GE.EQN%MODEL(iLayer,1)).AND.(r.LE.EQN%MODEL(iLayer+1,1)) ) THEN
                          IF( r.EQ.EQN%MODEL(iLayer+1,1) ) THEN
                             IF(rBary.LT.r) THEN
                                Mat_Vert(iVertex,1:nVar) = EQN%MODEL(iLayer+1,2:nVar+1)
                             ELSE
                                Mat_Vert(iVertex,1:nVar) = EQN%MODEL(iLayer+2,2:nVar+1)
                             ENDIF
                          ELSE
                             Mat_Vert(iVertex,1:nVar) = EQN%MODEL(iLayer,2:nVar+1) + ( EQN%MODEL(iLayer+1,2:nVar+1) - &
                                                         EQN%MODEL(iLayer,2:nVar+1) ) / ( EQN%MODEL(iLayer+1,1) - &
                                                         EQN%MODEL(iLayer,1) ) * ( r - EQN%MODEL(iLayer,1) )
                          ENDIF
                          EXIT
                       ENDIF
                    ENDDO
              END SELECT
              IF( Mat_Vert(iVertex,2) .NE. 0.0 ) THEN
                 Mat_Vert(iVertex,2) = 1.0 / Mat_Vert(iVertex,2)
              ELSE
                 mu_ct = mu_ct+1
              ENDIF
              IF( Mat_Vert(iVertex,3) .NE. 0.0 ) THEN
                 Mat_Vert(iVertex,3) = 1.0 / Mat_Vert(iVertex,3)
              ELSE
                 la_ct = la_ct+1
              ENDIF
           ENDDO
           ! arithmetic mean of rho
           Mat_Sum(1) = SUM(Mat_Vert(:,1),Dim=1)
           MaterialVal(iElem,1) = (1.0 / FLOAT(MESH%nVertexMax)) * Mat_Sum(1)
           IF( mu_ct.EQ.4 ) THEN ! if mu=0 at all vertices
              MaterialVal(iElem,2) = 0.0
           ELSE ! harmonic mean
              Mat_Sum(2) = 1.0 / SUM(Mat_Vert(:,2),Dim=1)
              MaterialVal(iElem,2) = FLOAT(MESH%nVertexMax) * Mat_Sum(2)
           ENDIF
           IF( la_ct.EQ.4 ) THEN ! if lambda=0 at all vertices
              MaterialVal(iElem,3) = 0.0
           ELSE ! harmonic mean
              Mat_Sum(3) = 1.0 / SUM(Mat_Vert(:,3),Dim=1)
              MaterialVal(iElem,3) = FLOAT(MESH%nVertexMax) * Mat_Sum(3)
           ENDIF
           ! arithmetic mean of Qp, Qs
           IF(EQN%Anelasticity.EQ.1)THEN
              Mat_Sum(4:5) = SUM(Mat_Vert(:,4:5),Dim=1)
              MaterialVal(iElem,4:5) = (1.0 / FLOAT(MESH%nVertexMax)) * Mat_Sum(4:5)
              MaterialTmp(:) = MaterialVal(iElem,:)
              IF((MaterialTmp(EQN%nAneMaterialVar-1).LT.99999.) &
                 .OR.(MaterialTmp(EQN%nAneMaterialVar).LT.99999.) )THEN
                 EQN%LocAnelastic(iElem) = 1
                 CALL ini_ATTENUATION(Theta,w_freq,Material_INF,MaterialTmp,EQN)
                 MaterialVal(iElem,2:EQN%AneMatIni-1) = Material_INF(:)
                 DO iMech = 1, EQN%nMechanisms
                    MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1))             = w_freq(iMech)
                    MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1)+1:EQN%AneMatIni+4*(iMech-1)+3) = Theta(iMech,:)
                 ENDDO
              ENDIF
           ENDIF
        ENDDO

      CASE(3)  !special case for layered medium linear variation of material parameters
        DO ielem=1,MESH%nElem
           x = MESH%ELEM%xyBary(1,iElem)
           y = MESH%ELEM%xyBary(2,iElem)
           z = MESH%ELEM%xyBary(3,iElem)
           FoundLayer = .FALSE.
           DO iLayer = 1, EQN%nLayers - 1
             IF( (z.LE.EQN%MODEL(iLayer,1)).AND.(z.GE.EQN%MODEL(iLayer+1,1)) ) THEN
               FoundLayer = .TRUE.
               MaterialVal(iElem,1:3) = EQN%MODEL(iLayer,2:4) + (EQN%MODEL(iLayer+1,2:4) - EQN%MODEL(iLayer,2:4)) /  &
                                                                (EQN%MODEL(iLayer+1,1)   - EQN%MODEL(iLayer,1)  ) *  &
                                                                ( z - EQN%MODEL(iLayer,1) )
               EXIT
             ENDIF
           ENDDO
           IF(.NOT.FoundLayer) THEN
             logError(*) 'Layered Medium Error: No layer found for depth', z
             STOP
           ENDIF
        ENDDO

      CASE(4)     !special case for Sismovalp 2D benchmark test (model M2)
        DO iElem = 1, MESH%nElem
           iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           IF(iLayer.GE.6) THEN                         ! Bedrock properties
               MaterialVal(iElem,:) = (/ 2500. , 1.96e10 , 2.84e10 /)
           ELSE
               y = -1.*MESH%ELEM%xyBary(2,iElem)
               MaterialVal(iElem,1) = 1600. + 59.5*(y)**(1./3.)
               MaterialVal(iElem,2) = (260. + 30*SQRT(y))**2. * MaterialVal(iElem,1)
               MaterialVal(iElem,3) = (525. + 60*SQRT(y))**2. * MaterialVal(iElem,1) - 2*MaterialVal(iElem,2)
           ENDIF
        ENDDO
      CASE(5)     !special case for Sismovalp 2D benchmark test (model M2 SH-wave)
        DO iElem = 1, MESH%nElem
           iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           IF(iLayer.GE.6) THEN                         ! Bedrock properties
               MaterialVal(iElem,:) = (/ 2500. , 0.0 , 1.96e10 /)
           ELSE
               y = -1.*MESH%ELEM%xyBary(2,iElem)
               MaterialVal(iElem,1) = 1600. + 59.5*(y)**(1./3.)
               MaterialVal(iElem,2) = 0.0
               MaterialVal(iElem,3) = (260. + 30*SQRT(y))**2. * MaterialVal(iElem,1)
           ENDIF
        ENDDO

      CASE(6)     !special case for Sismovalp 3D benchmark test
                  !layered bedrock with linear variation of material parameters
                  !sediment rock with depth dependent material parameters
                  ! interpolated linearly

        allocate( BaryDist(MESH%nElem) )
        CALL calc_BaryDist(BaryDist,OptionalFields,EQN,MESH,IO)

         DO iElem = 1, MESH%nElem

            iLayer = MESH%ELEM%Reference(0,iElem)                                   ! Zone number is given by reference 0

            z = BaryDist(iElem)
            depths = (/ MESH%ELEM%xyBary(3,iElem)+z , -3000. , -27000. , -35000. /) ! Depths of top of layers that are

            IF(iLayer.EQ.1) THEN                                                    ! Sediment properties

               EQN%LocAnelastic(iElem) = 1                                          ! Set local anelasticity

               cs  = 300.0  + 19.0  * SQRT(z)                                       ! Set sediment properties
               cp  = 1450.0 + 1.2   * (z)
               rho = 2140.0 + 0.125 * (z)
               MaterialTmp(1) = rho
               MaterialTmp(2) = cs**2 * rho
               MaterialTmp(3) = cp**2 * rho - 2*MaterialTmp(2)
               MaterialTmp(4) = 50. * 3./4. / ((cs/cp)**2)
               MaterialTmp(5) = 50.

               CALL ini_ATTENUATION(Theta,w_freq,Material_INF,MaterialTmp,EQN)

               MaterialVal(iElem,1:3) = MaterialTmp(1:3)                          ! Set elastic properties
               DO iMech = 1, EQN%nMechanisms                                      ! Set anelastic coefficients w_freq and theta
                  MaterialVal(iElem,iMech*4)             = w_freq(iMech)
                  MaterialVal(iElem,iMech*4+1:iMech*4+3) = Theta(iMech,:)
               ENDDO
            ELSE                                                                  ! Bedrock properties
               FoundLayer = .FALSE.
               DO i = 1, 3
                 IF( (MESH%ELEM%xyBary(3,iElem).LE.depths(i)).AND.(MESH%ELEM%xyBary(3,iElem).GE.depths(i+1)) ) THEN
                    FoundLayer = .TRUE.
                    MaterialVal(iElem,1:3) = EQN%MODEL(i,1:3) + (EQN%MODEL(i+1,1:3) - EQN%MODEL(i,1:3)) /  &
                                                                (depths(i+1)        - depths(i)       ) *  &
                                                                ( MESH%ELEM%xyBary(3,iElem) - depths(i) )
                    MaterialVal(iElem,4:EQN%nBackgroundVar) = 0.
                    EXIT
                 ENDIF
               ENDDO
               IF(.NOT.FoundLayer) THEN
                 logError(*) 'Layered Medium Error: No layer found for depth', MESH%ELEM%xyBary(3,iElem)
                 STOP
               ENDIF
            ENDIF

         ENDDO
         deallocate( BaryDist )

      CASE(7)     !special case for Euroseistest benchmark I2b (volvi lake)
                  !layered bedrock with linear variation of material parameters
                  !sediment rock with zone-dependent material parameters
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/ 10000.0, 2600.0, 1.75760e10, 1.74980e10/)
         BedrockVelModel(2,:) = (/ -1000.0, 2600.0, 1.75760e10, 1.74980e10/)
         BedrockVelModel(3,:) = (/ -3000.0, 2736.0, 3.20014e10, 3.84725e10/)
         BedrockVelModel(4,:) = (/ -5000.0, 2744.0, 3.22829e10, 3.85452e10/)
         BedrockVelModel(5,:) = (/ -7000.0, 2752.0, 3.25661e10, 3.86172e10/)
         BedrockVelModel(6,:) = (/ -9000.0, 2760.0, 3.28509e10, 3.86883e10/)
         BedrockVelModel(7,:) = (/-10000.0, 2772.0, 3.32813e10, 3.89645e10/)
         !
         DO iElem = 1, MESH%nElem
            EQN%LocAnelastic(iElem) = 1                                             ! Anelastic material everywhere
            iLayer = MESH%ELEM%Reference(0,iElem)                                   ! Zone number is given by reference 0
            !
            IF(iLayer.GE.4)THEN

               !Bedrock material linearly interpolated with Qp = infty and Qs = 260 everywhere
               FoundLayer = .FALSE.
               DO i = 1, 6
                 IF( (MESH%ELEM%xyBary(3,iElem).LE.BedrockVelModel(i,1)).AND.(MESH%ELEM%xyBary(3,iElem).GE.BedrockVelModel(i+1,1)) ) THEN
                    FoundLayer = .TRUE.
                    MaterialTmp(1:3) = BedrockVelModel(i,2:4) + (BedrockVelModel(i+1,2:4)  - BedrockVelModel(i,2:4)) /  &
                                                                (BedrockVelModel(i+1,1)    - BedrockVelModel(i,1))   *  &
                                                                (MESH%ELEM%xyBary(3,iElem) - BedrockVelModel(i,1))
                    MaterialTmp(4) = 1.0e6
                    MaterialTmp(5) = 260.0
                    EXIT
                 ENDIF
               ENDDO
               !
               CALL ini_ATTENUATION(Theta,w_freq,Material_INF,MaterialTmp,EQN)
               !
               MaterialVal(iElem,1) = MaterialTmp(1)                              ! Set rho
               MaterialVal(iElem,2:EQN%AneMatIni-1) = Material_INF(:)             ! Set unrelaxed material properties for this zone.                                                                      !
               ! Fill MaterialVal vector for each element with anelastic coefficients w_freq and theta
               DO iMech = 1, EQN%nMechanisms                                      ! Set anelastic coefficients w_freq and theta
                  MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1))                               = w_freq(iMech)
                  MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1)+1:EQN%AneMatIni+4*(iMech-1)+3) = Theta(iMech,:)
               ENDDO
               !
               IF(.NOT.FoundLayer) THEN
                 logError(*) 'Layered Medium Error: No layer found for depth', MESH%ELEM%xyBary(3,iElem)
                 STOP
               ENDIF
               !
            ELSE

               !Basin structure consists of three constant zones
               MaterialTmp(:) = EQN%MODEL(iLayer,:)                               ! Get material properties (temporary)
               !
               CALL ini_ATTENUATION(Theta,w_freq,Material_INF,MaterialTmp,EQN)    ! Initialize anelastic coefficients for this zone
               !
               MaterialVal(iElem,1) = MaterialTmp(1)                              ! Set rho
               MaterialVal(iElem,2:EQN%AneMatIni-1) = Material_INF(:)             ! Set unrelaxed material properties for this zone.                                                                      !
               ! Fill MaterialVal vector for each element with anelastic coefficients w_freq and theta
               DO iMech = 1, EQN%nMechanisms
                  MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1))                               = w_freq(iMech)
                  MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1)+1:EQN%AneMatIni+4*(iMech-1)+3) = Theta(iMech,:)
               ENDDO
            ENDIF
            !
         ENDDO
         !
      CASE(8)     !special case for Sonic logging
                                                                                    ! interpolated linearly
         DO iElem = 1, MESH%nElem
               iLayer = MESH%ELEM%Reference(0,iElem)
               radius = SQRT(MESH%ELEM%xyBary(1,iElem)**2+MESH%ELEM%xyBary(2,iElem)**2) ! Radial distance from origin
               IF( (radius.GT.0.10795).AND.(radius.LT.2.0) ) THEN

                    MaterialVal(iElem,1:3) = EQN%MODEL(2,1:3) + (EQN%MODEL(iLayer,1:3) - EQN%MODEL(2,1:3)) / &
                                            (2.0-0.10795) * (radius-0.10795)
               ELSE
                    MaterialVal(iElem,1:3) = EQN%MODEL(iLayer,1:3)
               ENDIF
         ENDDO
         !
      CASE(9)    ! special case for a hemisphere with different material properties at the top of a box
         !
         ! center of hemisphere in 0/0/0 - radius 100.0 m
         radius_ref = 100.0D0
         !
         DO iElem = 1, MESH%nElem
            !
            radius = SQRT(MESH%ELEM%xyBary(1,iElem)**2+MESH%ELEM%xyBary(2,iElem)**2+MESH%ELEM%xyBary(3,iElem)**2)
            IF (radius .LE. radius_ref) THEN
               MaterialVal(iElem,1:3) = EQN%MODEL(1,:)
            ELSE
               MaterialVal(iElem,1:3) = EQN%MODEL(2,:)
            ENDIF
         ENDDO
         !
      CASE(10)

        logError(*) 'LinType parameter ',EQN%linType,' only implemented for 2D square.'
        STOP
      !
      CASE(12) ! special case for TPV13: depth dependent initial shear/normal stresses are defined for the whole domain for the plastic calculations
               ! but otherwise constant material properties

        MaterialVal(:,1) = EQN%rho0
        MaterialVal(:,2) = EQN%mu
        MaterialVal(:,3) = EQN%lambda

        DO iElem=1, MESH%nElem

                z = MESH%ELEM%xyBary(3,iElem) !average depth inside an element
                IF (z.GE. -11951.15D0) THEN !down dip distance less than 13 800m
                        EQN%IniStress(1,iElem)  = -11242.17 * abs(z) !xx, sigma_2 = 0.5*(sigma_3 + sigma_1)
                        EQN%IniStress(2,iElem)  = -5824.34 * abs(z) !yy, sigma_3 = 0.3496*sigma_1
                        EQN%IniStress(3,iElem)  = -16660 * abs(z) !zz,sigma_1 from the benchmark description minus the fluid pressure of 9800
                        EQN%IniStress(4,iElem)  = 0.0  !shear stress components
                        EQN%IniStress(5,iElem)  = 0.0  !shear stress components
                        EQN%IniStress(6,iElem)  = 0.0  !shear stress components
                        !
                        !
                ELSE !down dip distance more than 13 800m
                        EQN%IniStress(1,iElem)  = -16660 * abs(z) !xx, sigma_2
                        EQN%IniStress(2,iElem)  = -16660 * abs(z) !yy, sigma_3
                        EQN%IniStress(3,iElem)  = -16660 * abs(z) !zz,sigma_1
                        EQN%IniStress(4,iElem)  = 0.0  !shear stress components
                        EQN%IniStress(5,iElem)  = 0.0  !shear stress components
                        EQN%IniStress(6,iElem)  = 0.0  !shear stress components
                     !
              ENDIF
        ENDDO
        !

        CASE(26)
        ! S.Wollherr 2016
        ! 26 = special case for TPV27
        ! depth dependent initial shear/normal stresses are defined for the whole domain for the plastic calculations
        ! but otherwise constant material properties

        ! Initialisation of IniStress(6 stress components in 3D)
        !

        MaterialVal(:,1) = EQN%rho0
        MaterialVal(:,2) = EQN%mu
        MaterialVal(:,3) = EQN%lambda


        b11 = 0.926793
        b33 = 1.073206
        b13 = -0.169029

        DO iElem=1, MESH%nElem

              z = MESH%ELEM%xyBary(3,iElem) !average depth inside an element
              Pf = 9800.0D0* abs(z) !fluid pressure, hydrostatic with water table at the surface

              IF (z.GE. -15000.0D0) THEN !depth less than 15000m
                 omega = 1.0D0
              ELSEIF ((z.LT. -15000.0D0) .AND. (z .GE. -20000.0D0) ) THEN !depth between 15000 and 20000m
                 omega = (20000.0D0-abs(z))/5000.0D0
              ELSE ! depth more than 20000m
                 omega = 0.0D0
              ENDIF



              EQN%IniStress(3,iElem)  = -2670D0*9.8D0 * abs(z) !zz
              EQN%IniStress(1,iElem)  = omega*(b11*(EQN%IniStress(3,iElem)+Pf)-Pf)+(1-omega)*EQN%IniStress(3,iElem) !xx
              EQN%IniStress(2,iElem)  = omega*(b33*(EQN%IniStress(3,iElem)+Pf)-Pf)+(1-omega)*EQN%IniStress(3,iElem) !yy

              EQN%IniStress(4,iElem)  = omega*(b13*(EQN%IniStress(3,iElem)+Pf))  !shear stress xy
              EQN%IniStress(5,iElem)  = 0.0  !shear stress xz
              EQN%IniStress(6,iElem)  = 0.0  !shear stress yz

              !add fluid pressure
              EQN%IniStress(1,iElem)  = EQN%IniStress(1,iElem) + Pf
              EQN%IniStress(2,iElem)  = EQN%IniStress(2,iElem) + Pf
              EQN%IniStress(3,iElem)  = EQN%IniStress(3,iElem) + Pf

        ENDDO


      CASE(33)     ! T. Ulrich TPV33 14.01.16
        DO iElem = 1, MESH%nElem
           !iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           y = MESH%ELEM%xyBary(2,iElem) !average y inside an element
           IF(y.LT.-800d0) THEN                         ! zone -800
               MaterialVal(iElem,1) = 2670.
               MaterialVal(iElem,2) = 2.816717568E+10
               MaterialVal(iElem,3) = 2.817615756E+10
           ELSEIF ((y.GE.-800d0).AND.(y.LE.800d0)) THEN                     ! zone central
               MaterialVal(iElem,1) = 2670.
               MaterialVal(iElem,2) = 1.251489075E+10
               MaterialVal(iElem,3) = 1.251709350E+10
           ELSEIF(y.GT.800d0) THEN                                          ! zone + 800
               MaterialVal(iElem,1) = 2670.
               MaterialVal(iElem,2) = 3.203812032E+10
               MaterialVal(iElem,3) = 3.204375936E+10
           ELSE
              logError(*) iLayer, ":zone (region) unknown"
           ENDIF
        ENDDO

      CASE(60) ! special case of 1D layered medium, imposed without meshed layers for Landers 1992
               ! after Wald and Heaton 1994, Table 1
               ! Note that mesh coordinates are in km, but the scaling matrix is used in read_mesh

         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/ -1500.0, 2300.0, 0.9017e10, 1.5178e10/)
         BedrockVelModel(2,:) = (/ -4000.0, 2600.0, 2.5798e10, 2.7053e10/)
         BedrockVelModel(3,:) = (/ -26000.0, 2700.0, 3.3454e10,3.6880e10/)
         BedrockVelModel(4,:) = (/ -32000.0, 2870.0, 4.2100e10, 4.8509e10/)
         BedrockVelModel(5,:) = (/ -50000.0, 3500.0, 7.5354e10, 7.3293e10/)
         BedrockVelModel(6,:) = (/ 0.0, 0.0, 0.0, 0.0/)
         BedrockVelModel(7,:) = (/ 0.0, 0.0, 0.0, 0.0/)
         !
         DO iElem = 1, MESH%nElem
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(1,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
             ELSEIF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
             ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
             ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(4,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSE
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ENDIF
         ENDDO
      !
      CASE(61) ! special case of 1D layered medium, imposed without meshed layers for Landers 1992
               ! simplified model after Graves and Pitarka 2010, Table 4
               ! values are averaged respecting the thickness layer, except for layer 1: see info sheet
               ! Note that mesh coordinates are in km, but the scaling matrix is used in read_mesh

         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/ -100.0, 2300.0, 0.1766e10, 0.4999e10/) ! not correctly averaged value to respect the low velocities somehow
         BedrockVelModel(2,:) = (/ -300.0, 2300.0, 0.6936e10, 1.3872e10/)
         BedrockVelModel(3,:) = (/ -1000.0, 2600.0, 1.3717e10, 1.8962e10/)
         BedrockVelModel(4,:) = (/ -3000.0, 2700.0, 2.1168e10, 2.7891e10/)
         BedrockVelModel(5,:) = (/ -6000.0, 2870.0, 3.1041e10, 3.8591e10/)
         BedrockVelModel(6,:) = (/ -31000.0, 3500.0, 3.9847e10, 4.3525e10/)
         BedrockVelModel(7,:) = (/ -50000.0, 3200.0, 6.4800e10, 6.5088e10/)
         !
         DO iElem = 1, MESH%nElem
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(1,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
             ELSEIF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
             ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
             ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(4,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
             ELSE
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ENDIF
        ENDDO
      !
       IF (EQN%Plasticity.EQ.1) THEN
         DO iElem=1, MESH%nElem

                z = MESH%ELEM%xyBary(3,iElem) !average depth inside an element

                EQN%IniStress(1,iElem)  = -10.6723e6*(abs(z-2000.0D0))/1000.0D0
                EQN%IniStress(2,iElem)  = -29.3277e6*(abs(z-2000.0D0))/1000.0D0
                EQN%IniStress(3,iElem)  = -20.0000e6*(abs(z-2000.0D0))/1000.0D0
                EQN%IniStress(4,iElem)  = -3.7687e6*(abs(z-2000.0D0))/1000.0D0
                EQN%IniStress(5,iElem)  =  0.0D0
                EQN%IniStress(6,iElem)  =  0.0D0
                        !
        ENDDO
      ENDIF !Plasticity
      !
      CASE(62)! new velocity model for Landers after Graves/Pitarka 2010 with average over the first layers respecting
              ! the thickness of the layer, added more layers in depth

         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/ -300.0, 2349.3, 0.5868e10, 1.1728e10/) ! not correctly averaged value to respect the low velocities somehow
         BedrockVelModel(2,:) = (/ -1000.0, 2592.9, 1.3885e10, 1.8817e10/)
         BedrockVelModel(3,:) = (/ -3000.0, 2700.0, 2.1168e10, 2.7891e10/)
         BedrockVelModel(4,:) = (/ -5000.0, 2750.0, 2.9948e10, 3.9105e10/)
         BedrockVelModel(5,:) = (/ -6000.0, 2800.0, 3.3327e10, 3.7534e10/)
         BedrockVelModel(6,:) = (/ -11000.0, 2825.0, 3.6612e10, 3.3625e10/)
         BedrockVelModel(7,:) = (/ -16000.0, 2850.0, 3.7969e10, 3.7898e10/)
         BedrockVelModel(8,:) = (/ -21000.0, 2900.0, 3.9701e10, 4.5015e10/)
         BedrockVelModel(9,:) = (/ -31000.0, 2950.0, 4.2598e10, 5.1212e10/)
         BedrockVelModel(10,:) = (/ -50000.0, 3200.0, 6.4800e10, 6.5088e10/)
         !

         DO iElem = 1, MESH%nElem
             z = MESH%ELEM%xyBary(3,iElem)
             IF (z.GT.BedrockVelModel(1,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
             ELSEIF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
             ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
             ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(4,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
             ELSEIF ((z.LT.BedrockVelModel(6,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSEIF ((z.LT.BedrockVelModel(7,1)).AND.(z.GE.BedrockVelModel(8,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(8,2:4)
             ELSEIF ((z.LT.BedrockVelModel(8,1)).AND.(z.GE.BedrockVelModel(9,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(9,2:4)
             ELSE
                 MaterialVal(iElem,1:3) =   BedrockVelModel(10,2:4)
             ENDIF

             !Plasticity initializations
             IF (EQN%Plasticity.EQ.1) THEN

                !stress tensor for the whole domain
                IF (z.LT. 1500.0D0) THEN
                    EQN%IniStress(1,iElem)  = EQN%Bulk_xx_0*(abs(z-2000.0D0))/1000.0D0
                    EQN%IniStress(2,iElem)  = EQN%Bulk_yy_0*(abs(z-2000.0D0))/1000.0D0
                    EQN%IniStress(3,iElem)  = EQN%Bulk_zz_0*(abs(z-2000.0D0))/1000.0D0
                    EQN%IniStress(4,iElem)  = EQN%ShearXY_0*(abs(z-2000.0D0))/1000.0D0
                    EQN%IniStress(5,iElem)  =  0.0D0
                    EQN%IniStress(6,iElem)  =  0.0D0
                ELSE ! constant stress tensor for everything higher than 1500m
                    EQN%IniStress(1,iElem)  = EQN%Bulk_xx_0*(abs(-500.0D0))/1000.0D0
                    EQN%IniStress(2,iElem)  = EQN%Bulk_yy_0*(abs(-500.0D0))/1000.0D0
                    EQN%IniStress(3,iElem)  = EQN%Bulk_zz_0*(abs(-500.0D0))/1000.0D0
                    EQN%IniStress(4,iElem)  = EQN%ShearXY_0*(abs(-500.0D0))/1000.0D0
                    EQN%IniStress(5,iElem)  =  0.0D0
                    EQN%IniStress(6,iElem)  =  0.0D0
                ENDIF   !

                ! depth dependent plastic cohesion
                EQN%PlastCo(iElem) = MaterialVal(iElem,2)/10000.0D0
            ENDIF !Plasticity


       ENDDO


      CASE(99) ! special case of 1D layered medium, imposed without meshed layers
      ! Northridge regional 1D velocity structure for sediments sites after Wald et al. 1996
         !
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/ 10000.0, 1700.0, 1.53e08, 7.82e08/)
         !BedrockVelModel(2,:) = (/ -100.0, 1800.0, 4.5e08, 1.692e09/)
         BedrockVelModel(2,:) = (/ -250.0, 1950.0, 1.097e+09, 2.4911e+09/)! include for 250-mesh
         !BedrockVelModel(2,:) = (/ -300.0, 2100.0, 2.1e09, 3.381e09/)
         BedrockVelModel(3,:) = (/ -500.0, 2400.0, 9.60e09, 1.920e10/)
         BedrockVelModel(4,:) = (/ -1500.0, 2700.0, 2.76480e10, 2.63790e10/)
         BedrockVelModel(5,:) = (/ -4000.0, 2800.0, 3.62880e10, 3.85560e10/)
         BedrockVelModel(6,:) = (/ -270000.0, 2900.0, 4.41090e10, 4.58780e10/)
         BedrockVelModel(7,:) = (/ -400000.0, 3300.0, 6.68250e10, 6.71220e10/)
         !
         DO iElem = 1, MESH%nElem
         z = MESH%ELEM%xyBary(3,iElem)
         IF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
         ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
         ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(4,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
         ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
         ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
         ELSEIF ((z.LT.BedrockVelModel(6,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
         ELSE
             MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
         ENDIF

         ENDDO
      !
      CASE(100) ! special case of 1D layered medium, imposed without meshed layers
      ! Northridge regional 1D velocity structure for rock sites after Wald et al. 1996
         !
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/ 10000.0, 2100.0, 2.10e09, 3.3810e09/)
         BedrockVelModel(2,:) = (/ -500.0, 2400.0, 9.60e09, 1.920e10/)
         BedrockVelModel(3,:) = (/ -1500.0, 2700.0, 2.76480e10, 2.63790e10/)
         BedrockVelModel(4,:) = (/ -4000.0, 2800.0, 3.62880e10, 3.85560e10/)
         BedrockVelModel(5,:) = (/ -270000.0, 2900.0, 4.41090e10, 4.58780e10/)
         BedrockVelModel(6,:) = (/ -400000.0, 3300.0, 6.68250e10, 6.71220e10/)
         !
         DO iElem = 1, MESH%nElem
         z = MESH%ELEM%xyBary(3,iElem)
         IF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
         ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
         ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(4,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
         ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
         ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
             MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
         ELSE
             MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
         ENDIF

         ENDDO
      !
      CASE(101) ! special case of 3D complex medium, imposed without meshed layers
        ! media properties given as structured grid, which can be smaller than domain size
        ! interpolation to elements
        !
        ! e.g. SCEC 3D velocity model surrounding the Northridge fault
        !
        ! ATTENTION: zones in the mesh are ignored
        !
        call readVelocityField(eqn, mesh, materialVal(:,1:3))

     CASE(122)     ! T. Ulrich SUMATRA 2 x 1d 16.02.16
	 ! OCeanic Crust
	 ! Layer                   depth    rho     mu          lambda
	 BedrockVelModel(1,:) = (/  -6d3, 2550d0,18589500000d0,26571000000d0/)
	 BedrockVelModel(2,:) = (/  -8d3, 2850d0,39016500000d0,42379500000d0/)
	 BedrockVelModel(3,:) = (/ -12d3, 3050d0,50027625000d0,53695250000d0/)
	 ! Crustal Crust
	 ! Layer                   depth    rho     mu          lambda
	 BedrockVelModel(4,:) = (/-6d3,2720d0,33320000000d0,31280000000d0/)
	 BedrockVelModel(5,:) = (/-12d3,2860d0,41298400000d0,41984800000d0/)
	 BedrockVelModel(6,:) = (/-23d3,3050d0,46390500000d0,60969500000d0/)
	 !below 1d layers
	 BedrockVelModel(7,:) = (/ -5d10, 3330d0,65942325000d0,81235350000d0/)

        DO iElem = 1, MESH%nElem
           iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           SELECT CASE (iLayer)
            CASE(1)
            ! OCeanic Crust
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(1,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
             ELSEIF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
             ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
             ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(3)
            ! Crustal Crust

             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(4,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
             ELSEIF ((z.LT.BedrockVelModel(6,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(2)
            MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
           CASE DEFAULT
                 logError(*) "Material assignement: unkown region", iLayer
           END SELECT
        ENDDO

     CASE(1222)     ! T. Ulrich SUMATRA 2 x 1d 16.02.16 and one layer below fault
         ! OCeanic Crust
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/  -6d3, 2550d0,18589500000d0,26571000000d0/)
         BedrockVelModel(2,:) = (/  -8d3, 2850d0,39016500000d0,42379500000d0/)
         BedrockVelModel(3,:) = (/ -12d3, 3050d0,50027625000d0,53695250000d0/)
         ! Crustal Crust
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(4,:) = (/-6d3,2720d0,33320000000d0,31280000000d0/)
         BedrockVelModel(5,:) = (/-12d3,2860d0,41298400000d0,41984800000d0/)
         BedrockVelModel(6,:) = (/-23d3,3050d0,46390500000d0,60969500000d0/)
         !below 1d layers
         BedrockVelModel(7,:) = (/ -5d10, 3330d0,65942325000d0,81235350000d0/)

        DO iElem = 1, MESH%nElem
           iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           SELECT CASE (iLayer)
            CASE(1)
             MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
            CASE(2)
            ! OCeanic Crust
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(1,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
             ELSEIF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
             ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
             ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(3)
            ! Crustal Crust
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(4,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
             ELSEIF ((z.LT.BedrockVelModel(6,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(4)
            MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
           CASE DEFAULT
                 logError(*) "Material assignement: unkown region", iLayer
           END SELECT
        ENDDO

     CASE(1223)     ! T. Ulrich SUMATRA 2 x 1d 16.02.16 and 3 layers below fault, 1 layer above: 7 regions
         ! OCeanic Crust
         ! Layer                   depth    rho     mu          lambda
         !BedrockVelModel(1,:) = (/  -6d3, 2550d0,18589500000d0,26571000000d0/)
         BedrockVelModel(1,:) = (/  -6d3, 2310d0, 7401471000d0,13494558000d0/)
         BedrockVelModel(2,:) = (/  -8d3, 2850d0,39016500000d0,42379500000d0/)
         BedrockVelModel(3,:) = (/ -12d3, 3050d0,50027625000d0,53695250000d0/)
         ! Crustal Crust
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(4,:) = (/-6d3,2720d0,33320000000d0,31280000000d0/)
         BedrockVelModel(5,:) = (/-12d3,2860d0,41298400000d0,41984800000d0/)
         BedrockVelModel(6,:) = (/-23d3,3050d0,46390500000d0,60969500000d0/)
         !below 1d layers
         BedrockVelModel(7,:) = (/ -5d10, 3330d0,65942325000d0,81235350000d0/)

        DO iElem = 1, MESH%nElem
           iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           SELECT CASE (iLayer)
            CASE(5)
             MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
            CASE(6)
             MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
            CASE(4)
             MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
            CASE(7)
             MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
            CASE(3)
            ! OCeanic Crust
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(1,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
             ELSEIF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
             ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
             ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(1)
            ! Crustal Crust
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(4,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
             ELSEIF ((z.LT.BedrockVelModel(6,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(2) !below 80km
            MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
           CASE DEFAULT
                 logError(*) "Material assignement: unkown region", iLayer
           END SELECT
        ENDDO

     CASE(1224)     ! T. Ulrich SUMATRA 2 x 1d 13.07.16 and 2 layers below fault, fault in a LVZ (small model)
         ! OCeanic Crust
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(1,:) = (/  -6d3, 2550d0,18589500000d0,26571000000d0/)
         !BedrockVelModel(1,:) = (/  -6d3, 2310d0, 7401471000d0,13494558000d0/)
         BedrockVelModel(2,:) = (/  -8d3, 2850d0,39016500000d0,42379500000d0/)
         BedrockVelModel(3,:) = (/ -12d3, 3050d0,50027625000d0,53695250000d0/)
         ! Crustal Crust
         ! Layer                   depth    rho     mu          lambda
         BedrockVelModel(4,:) = (/-6d3,2720d0,33320000000d0,31280000000d0/)
         BedrockVelModel(5,:) = (/-12d3,2860d0,41298400000d0,41984800000d0/)
         BedrockVelModel(6,:) = (/-23d3,3050d0,46390500000d0,60969500000d0/)
         !below 1d layers
         BedrockVelModel(7,:) = (/ -5d10, 3330d0,65942325000d0,81235350000d0/)

        DO iElem = 1, MESH%nElem
           iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           SELECT CASE (iLayer)
            CASE(2) !LVZ
             MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
            CASE(1)
             MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
            CASE(3)
             MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
            CASE(4)
             MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
           CASE(6)
            ! Crustal Crust
             z = MESH%ELEM%xyBary(3,iElem) ! supported by Sebs new mesh reader
             IF (z.GT.BedrockVelModel(4,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
             ELSEIF ((z.LT.BedrockVelModel(6,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(5) !below 80km
            MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
           CASE DEFAULT
                 logError(*) "Material assignement: unkown region", iLayer
           END SELECT
        ENDDO

     CASE(1221)     ! T. Ulrich SUMATRA 2 x 1d 09.03.2016 GEO MESH
	 ! OCeanic Crust
	 ! Layer                   depth    rho     mu          lambda
	 BedrockVelModel(1,:) = (/  -6d3, 2550d0,18589500000d0,26571000000d0/)
	 BedrockVelModel(2,:) = (/  -8d3, 2850d0,39016500000d0,42379500000d0/)
	 BedrockVelModel(3,:) = (/ -12d3, 3050d0,50027625000d0,53695250000d0/)
	 ! Crustal Crust
	 ! Layer                   depth    rho     mu          lambda
	 BedrockVelModel(4,:) = (/-6d3,2720d0,33320000000d0,31280000000d0/)
	 BedrockVelModel(5,:) = (/-12d3,2860d0,41298400000d0,41984800000d0/)
	 BedrockVelModel(6,:) = (/-23d3,3050d0,46390500000d0,60969500000d0/)
	 !below 1d layers
	 BedrockVelModel(7,:) = (/ -5d10, 3330d0,65942325000d0,81235350000d0/)

        DO iElem = 1, MESH%nElem
           iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
           SELECT CASE (iLayer)
            CASE(1)
            ! OCeanic Crust
             ! R is taken at lon, lat  = (90,8)
             z = sqrt(MESH%ELEM%xyBary(1,iElem)**2+MESH%ELEM%xyBary(2,iElem)**2+MESH%ELEM%xyBary(3,iElem)**2)-6377726.19283
             IF (z.GT.BedrockVelModel(1,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(1,2:4)
             ELSEIF ((z.LT.BedrockVelModel(1,1)).AND.(z.GE.BedrockVelModel(2,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(2,2:4)
             ELSEIF ((z.LT.BedrockVelModel(2,1)).AND.(z.GE.BedrockVelModel(3,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(3,2:4)
             ELSEIF ((z.LT.BedrockVelModel(3,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(3)
            ! Crustal Crust

             ! R is taken at lon, lat  = (90,8)
             z = sqrt(MESH%ELEM%xyBary(1,iElem)**2+MESH%ELEM%xyBary(2,iElem)**2+MESH%ELEM%xyBary(3,iElem)**2)-6377726.19283
             IF (z.GT.BedrockVelModel(4,1)) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(4,2:4)
             ELSEIF ((z.LT.BedrockVelModel(4,1)).AND.(z.GE.BedrockVelModel(5,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(5,2:4)
             ELSEIF ((z.LT.BedrockVelModel(5,1)).AND.(z.GE.BedrockVelModel(6,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(6,2:4)
             ELSEIF ((z.LT.BedrockVelModel(6,1)).AND.(z.GE.BedrockVelModel(7,1))) THEN
                 MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
             ELSE
                 logError(*) "depth lower than",BedrockVelModel(7,1),iLayer,z
             ENDIF
           CASE(2)
            MaterialVal(iElem,1:3) =   BedrockVelModel(7,2:4)
           CASE DEFAULT
                 logError(*) "Material assignement: unkown region", iLayer
           END SELECT
        ENDDO

      CASE DEFAULT
        logError(*) 'Wrong linType for elastic wave equations.'
        STOP
      END SELECT

      !###################################################################################!
      !  Apply random field perturbation
      !###################################################################################!

      IF (EQN%RandomField_Flag.GT.0) THEN

          DO iRFFlag = 1,EQN%RandomField_Flag
             logInfo(*) 'Random Field perturbation is read from file : ',TRIM(IO%RF_Files(iRFFlag))
             CALL OpenFile(                                        &
                   UnitNr       = IO%UNIT%other01                , &
                   Name         = IO%RF_Files(iRFFlag)           , &
                   create       = .FALSE.                          )

             READ(IO%UNIT%other01,*) zoneNum, nPert               ! zoneNum specifies on which zone this perturbation is applied
             IF(zoneNum.GT.EQN%nLayers) THEN
                 logError(*) 'Zone ',zoneNum,' of random field ',TRIM(IO%RF_Files(iRFFlag)),   &
                             ' does not exist in the model!'
                 STOP
             ENDIF
             ALLOCATE( PertMaterial(nPert) )                      ! nPert   specifies on how many material parameters it is applied
                                                                  ! PertMaterial are the indices of the material parameters as
             READ(IO%UNIT%other01,*) NX, NY, NZ, PertMaterial(:)  ! ordered in MaterialVal

             ALLOCATE( xrf(NX), yrf(NY), zrf(NZ), pertrf(NX,NY,NZ,nPert) )
             pertrf = 0.0d0

             DO i=1,NX
               DO j=1,NY
                 DO k=1,NZ
                    READ(IO%UNIT%other01,*) xrf(i),yrf(j),zrf(k),pertrf(i,j,k,:)
                 ENDDO
               ENDDO
             ENDDO

             CLOSE(IO%UNIT%other01)

             xrf_max = MAXVAL(xrf(:))
             xrf_min = MINVAL(xrf(:))
             yrf_max = MAXVAL(yrf(:))
             yrf_min = MINVAL(yrf(:))
             zrf_max = MAXVAL(zrf(:))
             zrf_min = MINVAL(zrf(:))
             DO i = 1,nPert
                k = PertMaterial(i)
                pert_max = MAXVAL(pertrf(:,:,:,i))
                pert_min = MINVAL(pertrf(:,:,:,i))
                logInfo(*) 'Perturbation_max for MaterialValue ',k, ': ',pert_max
                logInfo(*) 'Perturbation_min for MaterialValue ',k, ': ',pert_min
             ENDDO
             logInfo('(" |   Perturbation domain x = [",E15.5," , ",E15.5,"]")')  xrf_min, xrf_max
             logInfo('(" |   Perturbation domain y = [",E15.5," , ",E15.5,"]")')  yrf_min, yrf_max
             logInfo('(" |   Perturbation domain z = [",E15.5," , ",E15.5,"]")')  zrf_min, zrf_max
             IF (InterpolationScheme.EQ.1) THEN
                logInfo('("Linear interpolation scheme is used.")')
             ELSE
                logInfo('("Cubic interpolation scheme is used.")')
             ENDIF

             ! If we use constant material apply perturbation to barycenter
             IF(DISC%Galerkin%nPolyMatOrig.EQ.0) THEN

                 ALLOCATE( posx(MESH%nElem), posy(MESH%nElem), posz(MESH%nElem), pert(MESH%nElem) )
                 counter = 0

                 DO iElem = 1,MESH%nElem
                     IF (MESH%ELEM%Reference(0,iElem).EQ.zoneNum) THEN
                        counter = counter + 1
                        posx(counter) = MESH%ELEM%xyBary(1,iElem)
                        posy(counter) = MESH%ELEM%xyBary(2,iElem)
                        posz(counter) = MESH%ELEM%xyBary(3,iElem)
                     ENDIF
                 ENDDO

                 posx_max = MAXVAL(posx(1:counter))
                 posx_min = MINVAL(posx(1:counter))
                 posy_max = MAXVAL(posy(1:counter))
                 posy_min = MINVAL(posy(1:counter))
                 posz_max = MAXVAL(posz(1:counter))
                 posz_min = MINVAL(posz(1:counter))

                 IF( posx_max.GT.xrf_max .OR. posy_max.GT.yrf_max .OR. posz_max.GT.zrf_max .OR.  &
                     posx_min.LT.xrf_min .OR. posy_min.LT.yrf_min .OR. posz_min.LT.zrf_min )THEN
                         logError(*) 'Random field ',TRIM(IO%RF_Files(iRFFlag)),   &
                                     ' does not fully include zone ',zoneNum
                     STOP
                 ENDIF

                 ALLOCATE( PerturbationVar(counter,nPert) )
                 logInfo(*) 'Interpolating random field ... for zone  ', zoneNum
                 !
                 ! Interpolation of the material perturbation of the
                 ! given random field onto the barycenter of the element
                 !
                 DO iMaterial = 1,nPert

                   IF (InterpolationScheme.EQ.1) THEN
                       ! Linear interpolation
                       CALL TrilinearFromRaster(pert(1:counter),posx(1:counter),posy(1:counter),posz(1:counter),       &
                                                pertrf(:,:,:,iMaterial),xrf(2)-xrf(1),yrf(2)-yrf(1),zrf(2)-zrf(1),xrf_min,yrf_min,zrf_min)
                   ELSE
                       ! Bicubic interpolation
                       CALL var3Dipol(pertrf(:,:,:,iMaterial),xrf,yrf,zrf,NX,NY,NZ,       &
                                      posx(1:counter),posy(1:counter),posz(1:counter),pert(1:counter),counter)
                   ENDIF

                   PerturbationVar(1:counter,iMaterial) = pert(1:counter)

                 ENDDO
                 !
                 counter = 0
                 DO iElem = 1,MESH%nElem
                    IF (MESH%ELEM%Reference(0,iElem).EQ.zoneNum) THEN
                       counter = counter + 1
                       DO iMaterial = 1,nPert
                           k = PertMaterial(iMaterial)
                           MaterialVal(iElem,k) = MaterialVal(iElem,k)+PerturbationVar(counter,iMaterial)
                            IF (MaterialVal(iElem,k).LE.0.0d0) THEN
                                logError('("  Error in interpolating variable ",I3)') iMaterial
                                logError('("       Element ",I6)') iElem
                                STOP
                            ENDIF
                       ENDDO
                    ENDIF
                 ENDDO

                 DEALLOCATE( posx, posy, posz, pert )
                 DEALLOCATE( PertMaterial, PerturbationVar )

             ENDIF ! DISC%Galerkin%nPolyMatOrig.EQ.0

             DEALLOCATE( xrf, yrf, zrf, pertrf )
          ENDDO ! iRFFlag = 1,EQN%RandomField_Flag

      ENDIF ! EQN%RandomField_Flag.GT.0

      !###################################################################################!
      !  Dynamic Rupture setup
      !###################################################################################!

      IF(EQN%DR.EQ.1) THEN

        CALL DR_setup(EQN,DISC,MESH,IO,BND)

      ENDIF ! EQN%DR.EQ.1

      ! Call the post model hooks
      call call_hook_post_model()

  END SUBROUTINE ini_MODEL


  SUBROUTINE ini_ATTENUATION(Theta,w_out,MatVal_INF,MaterialTmp,EQN)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tEquations)        :: EQN
    TYPE (tUnstructMesh)     :: MESH
    TYPE (tInputOutput)      :: IO
    TYPE (tInitialCondition) :: IC
    INTEGER                  :: i,j, counter
    INTEGER                  :: kmax                        ! Number of equations to solve (the system is overdetermined)
    REAL, POINTER            :: AP(:,:)   , AS(:,:)         ! System matrices for P- and S-waves
    REAL, POINTER            :: QPInv(:)  , QSInv(:)        ! Inverted values of QP and QS for each mechanism
    REAL, POINTER            :: Y_alpha(:), Y_beta(:)
    REAL                     :: Y_kappa, Y_mu
    REAL, POINTER            :: WM(:), VM(:,:), SolVec(:)
    REAL                     :: MaterialTmp(EQN%nAneMaterialVar), MatVal_INF(EQN%nAneMaterialVar-3)
    REAL                     :: cp, cs
    REAL                     :: w0, wmin, R, PSI1P, PSI2P, PSI1S, PSI2S
    REAL                     :: P_INF, mu_INF, Lambda_INF
    REAL                     :: QPLocVal,QSLocVal
    REAL, POINTER            :: w(:)
    REAL                     :: Theta(EQN%nMechanisms,3),w_out(EQN%nMechanisms)
    REAL                     :: AvP, AvMu, AvLambda     ! Averaged elastic properties
    ! Local variable declaration
    !--------------------------------------------------------------------------
    INTENT(IN)               :: MaterialTmp,EQN
    INTENT(OUT)              :: Theta, w_out, MatVal_INF
    ! -------------------------------------------------------------------------

    kmax    = 2.*EQN%nMechanisms - 1

    ALLOCATE( AP(kmax,EQN%nMechanisms),                  &
              AS(kmax,EQN%nMechanisms),                  &
              QPInv(kmax),                               &
              QSInv(kmax),                               &
              Y_alpha(EQN%nMechanisms),                  &
              Y_beta(EQN%nMechanisms),                   &
              w(kmax)                                    )
    ALLOCATE( WM(EQN%nMechanisms),                       &
              VM(EQN%nMechanisms,EQN%nMechanisms)        )

    ! Selection of the logarithmically equispaced frequencies
    ! i.e.: log(w) = log(wmin) + (i-1)*log(f_ratio)/(2*(n-1))

    w0     = 2.* EQN%Pi * EQN%FreqCentral
    wmin   = w0 / SQRT(EQN%FreqRatio)

    IF(EQN%nMechanisms.GT.1)THEN
       DO i = 1, kmax
           w(i) = EXP( LOG(wmin) + (i-1)/(2.*(EQN%nMechanisms-1)) * LOG(EQN%FreqRatio) )
       ENDDO
    ELSE
       w(:) = w0
    ENDIF

    QPLocVal = MaterialTmp(EQN%nAneMaterialVar-1)
    QSLocVal = MaterialTmp(EQN%nAneMaterialVar)

    counter = 1
    DO i = 1, kmax, 2
       w_out(counter) = w(i)
       counter = counter + 1
    ENDDO

    ! Build the system matrices
    DO i = 1, kmax
       DO j = 1, EQN%nMechanisms
           AP(i,j) = (w(2*j-1)*w(i)+w(2*j-1)**2 / QPLocVal) / (w(2*j-1)**2+w(i)**2)
       ENDDO
    ENDDO

    DO i = 1, kmax
       DO j = 1, EQN%nMechanisms
           AS(i,j) = (w(2*j-1)*w(i)+w(2*j-1)**2 / QSLocVal) / (w(2*j-1)**2+w(i)**2)
       ENDDO
    ENDDO

    ! Build the right hand sides
    DO i = 1, kmax
       QPInv(i) = 1./QPLocVal                 !Desired values of Q for each mechanism inverted (P wave)
    ENDDO

    DO i = 1, kmax
       QSInv(i) = 1./QSLocVal                 !Desired values of Q for each mechanism inverted (S wave)
    ENDDO

    ! Solving the overdetermined system for the anelastic coefficients Y_alpha and Y_beta
    VM = 0.
    WM = 0.

    CALL svdecomp(WM,VM,AP)                          ! Decompose matrix AP
    CALL svbacksb(WM,VM,AP,QPInv,Y_alpha)

    CALL svdecomp(WM,VM,AS)                          ! Decompose matrix AS
    CALL svbacksb(WM,VM,AS,QSInv,Y_beta)

    SELECT CASE(EQN%Anisotropy)
    CASE(0)
      !
      ! Computing unrelaxed moduli lambda and mu that give the
      ! specified wave velocities at the corresponding frequency
      !
      PSI1P = 1.
      PSI2P = 0.
      DO i = 1, EQN%nMechanisms
         PSI1P = PSI1P - Y_alpha(i) / (1+(w0/w(2*i-1))**2)
         PSI2P = PSI2P + Y_alpha(i) *     w0/w(2*i-1)/(1+(w0/w(2*i-1))**2)
      ENDDO
      R     = SQRT(PSI1P**2 + PSI2P**2)
      P_INF = (MaterialTmp(3)+2*MaterialTmp(2)) * (R + PSI1P)/(2*R**2)
      !
      PSI1S = 1.
      PSI2S = 0.
      DO i = 1, EQN%nMechanisms
         PSI1S = PSI1S - Y_beta(i)  / (1+(w0/w(2*i-1))**2)
         PSI2S = PSI2S + Y_beta(i)  *     w0/w(2*i-1)/(1+(w0/w(2*i-1))**2)
      ENDDO
      R    = SQRT(PSI1S**2 + PSI2S**2)
      ! unrelaxed moduli
      mu_INF     = MaterialTmp(2) * (R + PSI1S)/(2*R**2)
      Lambda_INF = P_INF - 2*mu_INF
      !
      MatVal_INF(1) = mu_INF
      MatVal_INF(2) = Lambda_INF
      !
      ! Compute Thetas for source term E
      DO i=1,EQN%nMechanisms
         ! Note to future self: Y_alpha != Y_lambda
         Theta(i,1) = - (Lambda_INF+2.d0*mu_INF) * Y_alpha(i)
         Theta(i,2) = - (Lambda_INF+2.d0*mu_INF) * Y_alpha(i) + 2.d0*mu_INF * Y_beta(i)
         Theta(i,3) = -2.d0 * mu_INF * Y_beta(i)
      ENDDO
      !
    CASE(1)
      !
      ! Compute Thetas for source term E
        ! Compute the average P (= lam + 2mu = 1/3*(c11+c22+c33)) and mu (= 1/3*(c44+c55+c66))
        AvP  = ( MaterialTmp(2) + MaterialTmp(8) + MaterialTmp(13)) / 3.d0
        AvMu     = ( MaterialTmp(17) + MaterialTmp(20) + MaterialTmp(22)) / 3.d0
        AvLambda = AvP - 2.d0 * AvMu
      !
      DO i=1,EQN%nMechanisms
           Theta(i,1) = - (AvP) * Y_alpha(i)
           Theta(i,2) = - (AvP) * Y_alpha(i) + 2.d0*AvMu * Y_beta(i)
           Theta(i,3) = - 2.d0*AvMu * Y_beta(i) !Remember that will have to be multiplied either by c44, c55 or c66 (in 3D)!!
      ENDDO
    END SELECT

  END SUBROUTINE ini_ATTENUATION


  SUBROUTINE generate_FacetList(iObject,iBNDType,OptionalFields,EQN,MESH,IO)

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tInputOutput)              :: IO
    TYPE(tUnstructOptionalFields)   :: OptionalFields
    INTEGER                         :: iObject
    INTEGER                         :: iBNDType
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i,j,k,iElem, iSide, iFacet, iIntGP
    INTEGER                         :: iSponge
    INTEGER                         :: iBNDElem, iBNDSide
    INTEGER                         :: counter
    INTEGER                         :: nearestFacet
    INTEGER                         :: SpongeList(MESH%nElem)
    INTEGER                         :: SpongeFacet(MESH%nElem,2)
    INTEGER                         :: VertexSide(4,EQN%DIMENSION)
    INTEGER                         :: NMAX
    INTEGER                         :: minl(1)
    REAL                            :: xP,yP,zP,dr,dphi,r1,r2, r,phi,state(5),time
    REAL, POINTER                   :: distance(:)
    REAL                            :: Factor
    REAL                            :: xGP,yGP,zGP
    REAL                            :: SpongeDistance(MESH%nElem)
    REAL                            :: Node1(EQN%DIMENSION), Node2(EQN%DIMENSION), Node3(EQN%DIMENSION)
    REAL                            :: Vec1(EQN%DIMENSION), Vec2(EQN%DIMENSION), Normal(EQN%DIMENSION)
    ! ------------------------------------------------------------------------!
    INTENT(IN)    :: iObject, iBNDType, EQN, MESH, IO
    INTENT(INOUT) :: OptionalFields
    ! ------------------------------------------------------------------------!
    !
        ! The unit tetrahedron has the following 4 local vertices:                !
        !                                                                         !
        ! 1 = (0,0,0)                                                             !
        ! 2 = (1,0,0)                                                             !
        ! 3 = (0,1,0)                                                             !
        ! 4 = (0,0,1)                                                             !
        !                                                                         !
        VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
        VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
        VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
        VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !
    !
    logInfo(*) 'Entering generate_FacetList...'
    !
    SpongeList(:)     = -1
    SpongeDistance(:) = 1e20
    SpongeFacet(:,:)  = -1
    !
    ! Identify boundary facets for this sponge
    !
    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
         IF(MESH%ELEM%Reference(iSide,iElem).EQ.iBNDType) THEN
           IF(MESH%ELEM%BoundaryToObject(iSide,iElem).EQ.iObject) THEN
             counter = counter + 1
           ENDIF
         ENDIF
      ENDDO
    ENDDO
    OptionalFields%nFacet = counter
    !
    logInfo(*) 'Free surface ', ' has ', OptionalFields%nFacet, ' boundary facets. '
    !
    ALLOCATE(OptionalFields%FacetList(OptionalFields%nFacet,2))
    ALLOCATE(OptionalFields%FacetBary(OptionalFields%nFacet,EQN%DIMENSION))
    counter = 0
    DO iElem = 1, MESH%nElem
      DO iSide = 1, MESH%GlobalElemType
         IF(MESH%ELEM%Reference(iSide,iElem).EQ.iBNDType) THEN
           IF(MESH%ELEM%BoundaryToObject(iSide,iElem).EQ.iObject) THEN
              counter = counter + 1
              OptionalFields%FacetList(counter,1) = iElem
              OptionalFields%FacetList(counter,2) = iSide
              Node1(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem))
              Node2(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,2),iElem))
              Node3(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,3),iElem))
              OptionalFields%FacetBary(counter,:) = 1./3.*( Node1(:)+Node2(:)+Node3(:) )
           ENDIF
         ENDIF
      ENDDO
    ENDDO
    !
    logInfo(*) 'Leaving generate_FacetList...'
    !
  END SUBROUTINE generate_FacetList

  SUBROUTINE calc_baryDist(BaryDist,OptionalFields,EQN,MESH,IO)

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tInputOutput)              :: IO
    TYPE (tUnstructOptionalFields)  :: OptionalFields
    REAL, allocatable, dimension(:) :: BaryDist
    INTEGER                         :: iObject
    INTEGER                         :: iBNDType
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i,j,k,iElem, iSide, iFacet, iIntGP
    INTEGER                         :: iSponge
    INTEGER                         :: iBNDElem, iBNDSide
    INTEGER                         :: counter
    INTEGER                         :: nearestFacet
    INTEGER                         :: SpongeList(MESH%nElem)
    INTEGER                         :: SpongeFacet(MESH%nElem,2)
    INTEGER                         :: VertexSide(4,EQN%DIMENSION)
    INTEGER                         :: NMAX
    INTEGER                         :: minl(1)
    REAL                            :: xP,yP,zP,dr,dphi,r1,r2, r,phi,state(5),time
    REAL, POINTER                   :: distance(:)
    REAL                            :: Factor
    REAL                            :: x(MESH%GlobalVrtxType) ! Element vertices in x-y space         !
    REAL                            :: y(MESH%GlobalVrtxType) ! Element vertices in x-y space         !
    REAL                            :: z(MESH%GlobalVrtxType) ! Element vertices in x-y space         !
    REAL                            :: xGP,yGP,zGP
    REAL                            :: SpongeDistance(MESH%nElem)
    REAL                            :: Node1(EQN%DIMENSION), Node2(EQN%DIMENSION), Node3(EQN%DIMENSION)
    REAL                            :: Vec1(EQN%DIMENSION), Vec2(EQN%DIMENSION), Normal(EQN%DIMENSION)
    ! ------------------------------------------------------------------------!
    INTENT(IN)    :: EQN, MESH, IO
    INTENT(INOUT) :: BaryDist
    INTENT(INOUT) :: OptionalFields
    ! ------------------------------------------------------------------------!
    !
    logInfo(*) 'Entering calc_baryDist...'
    !
    IF(OptionalFields%nFacet.LE.0) THEN
        logError(*) 'OptionalFields%nFacet <= 0. '
        STOP
    ENDIF
    !
    ALLOCATE( Distance(OptionalFields%nFacet) )
    !
    DO iElem = 1, MESH%nElem
      Node1(:) = MESH%ELEM%xyBary(:,iElem)
      nearestFacet = -1
      DO iFacet = 1, OptionalFields%nFacet
        Node2(:) = OptionalFields%FacetBary(iFacet,:)
        distance(iFacet) = (Node1(1)-Node2(1))**2 + (Node1(2)-Node2(2))**2
      ENDDO
      !
      minl = MINLOC(distance)
      nearestFacet = minl(1)
      !
      IF(nearestFacet.EQ.-1) THEN
         logError(*) 'Facet Error in calc_BaryDepth for element ', iElem
         STOP
      ENDIF
      !
      BaryDist(iElem) = ABS(OptionalFields%FacetBary(nearestFacet,3)-MESH%ELEM%xyBary(3,iElem))
      !
      IF(MESH%nElem.GT.20) THEN
         IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
             logInfo(*) iElem,' elements done...'
         END IF
      ENDIF
      !
    ENDDO
    !
    DEALLOCATE(OptionalFields%FacetList, OptionalFields%FacetBary, Distance)
    !
  END SUBROUTINE calc_baryDist

  SUBROUTINE var3Dipol(model,x,y,z,NX,NY,NZ,posx,posy,posz,perturbation,NRES)
    IMPLICIT NONE
    ! Argument list declaration
    integer :: i,j,NX,NY,NZ,NRES
    real :: model(NX,NY,NZ),x(NX),y(NY),z(NZ),posx(NRES),posy(NRES),posz(NRES)
    real :: perturbation(NRES),y2a(NX,NY,NZ),za(NZ),y2(NZ)
    INTENT(OUT):: perturbation
    DO i=1,NZ
       y2a(:,:,i) = spline_deriv2(x,y,model(:,:,i))
    ENDDO
    DO j=1,NRES
       DO i=1,NZ
          za(i) =spline_interp2(posx(j),posy(j),x,y,model(:,:,i),y2a(:,:,i))
       ENDDO
       perturbation(j) = spline_interp(posz(j),z,za,spline_deriv(1.e30,1.e30,z,za))
    ENDDO
  END SUBROUTINE var3Dipol

  pure function spline_deriv2(x1,x2,fx)result(dfx)
    implicit none
    real(rk),intent(in)::x1(:),x2(:),fx(size(x1),size(x2))
    real(rk) :: dfx(size(x1),size(x2))
    integer :: j
    forall(j=1:size(x1))dfx(j,:)=spline_deriv(fx_large,fx_large,x2,fx(j,:))
  end function spline_deriv2

  pure function spline_interp2(x_1,x_2,x1,x2,fx,fx2)result(sfx)
    implicit none
    real(rk),intent(in)::x_1,x_2,x1(:),x2(:),fx2(size(x1),size(x2)),fx(size(x1),size(x2))
    real(rk) :: sfx,tmp(size(x2))
    integer :: i
    forall(i=1:size(x1)) tmp(i)=spline_interp(x_2,x2,fx(i,:),fx2(i,:))
    sfx = spline_interp(x_1,x1,tmp,spline_deriv(fx_large,fx_large,x1,tmp))
  end function spline_interp2

  pure function spline_deriv(dfx1,dfxn,x,fx)result(dfx)
    implicit none
    real(rk), intent(in) :: dfx1,dfxn,x(:),fx(size(x))
    real(rk) :: dfx(size(x)),C,CN,coeff,W1,W(size(x))
    integer  :: i,n
    n=size(x)
    if(dfx1 < fx_large)then ;
       W(1)=(3._rk/(x(2)-x(1)))*((fx(2)-fx(1))/(x(2)-x(1))-dfx1); dfx(1)=-0.5_rk;
    else; W(1)=0._rk;dfx(1)=0._rk;
    end if
    do i=2,n-1
       coeff=(x(i)-x(i-1))/(x(i+1)-x(i-1)); C=dfx(i-1)*coeff+2._rk;dfx(i)=(coeff-1._rk)/C
       W(i)=(6._rk*((fx(i+1)-fx(i))/(x(i+1)-x(i)) &
            &      -(fx(i)-fx(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-coeff*W(i-1))/C
    enddo
    if(dfxn < fx_large)then;
       W1=(3._rk/(x(n)-x(n-1)))*(dfxn-(fx(n)-fx(n-1))/(x(n)-x(n-1)));CN=0.5_rk;
    else; W1=0._rk;CN=0._rk;
    end if
    dfx(n) = (W1 - CN*W(n-1))/(1._rk + CN*dfx(n-1))
    do i=n-1,1,-1
       dfx(i)=dfx(i)*dfx(i+1)+W(i)
    enddo
  end function spline_deriv

  pure function spline_interp(x1,x,fx,fx2)result(sfx)
    implicit none
    real(rk), intent(in)::x1,x(:),fx2(size(x)),fx(size(x))
    real(rk) :: sfx,der1,der2,dx
    integer :: myi,i2,i1
    real(rk),parameter :: one_o_6=real(1,rk)/real(6,rk)
    i1=1;i2=size(x)
    do while(i2-i1 > 1)
       myi=(i2+i1)/2
       if(x(myi)>x1)then;i2=myi
       else;i1=myi
       endif
    enddo
    dx = x(i2)-x(i1); !if(abs(dx)<epsilon(1._rk))stop 'error: spline_interp: wrong dx'
    der1 = (x(i2)-x1)/dx; der2 = (x1-x(i1))/dx
    sfx = fx(i1)*der1 + fx(i2)*der2 &
         &  + (dx*dx*one_o_6)*((der1*der1*der1-der1)*fx2(i1) &
         &  +                  (der2*der2*der2-der2)*fx2(i2))
  end function spline_interp

 SUBROUTINE ModelDefinition_new(A,B,C,E,xGP,yGP,zGP,MaterialVal,                   &
                                LocAnelastic,LocPoroelastic,EQN,MESH,DISC,SOURCE,IO)

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tDiscretization)           :: DISC
    TYPE(tSource)                   :: Source
    TYPE(tInputOutput)              :: IO
    REAL                            :: A(:,:)
    REAL                            :: B(:,:)
    REAL                            :: C(:,:)
    REAL                            :: E(:,:)
    REAL                            :: xGP,yGP,zGP
    REAL                            :: MaterialVal(:)
    INTEGER                         :: LocAnelastic
    INTEGER                         :: LocPoroelastic
    ! Local variable declaration
    REAL                            :: kx,ky,kz
    REAL                            :: rho,lambda,mu,cp
    REAL                            :: AvMu
    REAL                            :: ce(6,6)                                ! Elastic material constants
    REAL                            :: K_F, K_S, K_Mean, MM, Poro, Alpha(6)   !Porous parameters
    REAL                            :: rho_F, rho_S, nu, Kappa(3), Tor(3)     !Porous parameters
    REAL                            :: Rho1(3), Rho2(3), Beta1(3), Beta2(3)   !Porous parameters

    INTEGER                         :: i, j, k

    !--------------------------------------------------------------------------
    INTENT(IN)                      :: xGP,yGP,zGP,MaterialVal
    INTENT(IN)                      :: LocAnelastic,LocPoroelastic
    INTENT(IN)                      :: EQN,MESH,DISC,SOURCE,IO
    INTENT(INOUT)                   :: A,B,C,E
    ! -------------------------------------------------------------------------
    !
    IF (EQN%LinType.EQ.0 .AND. SOURCE%Type.EQ.1) THEN

        A = 0.
        B = 0.
        C = 0.
        E = 0.

        rho    = EQN%rho0  ! rho = 1, mu = 2, lambda = 3
        mu     = EQN%mu
        lambda = EQN%lambda
        cp     = SQRT((lambda+2*mu)/rho)
        !
        ! Jacobian in x-direction
        !
        A(1,7) = -lambda-2*mu
        A(2,7) = -lambda
        A(3,7) = -lambda
        A(4,8) = -mu
        A(6,9) = -mu
        A(7,1) = -1./rho
        A(8,4) = -1./rho
        A(9,6) = -1./rho
        !
        ! Jacobian in y-direction
        !
        B(1,8) = -lambda
        B(2,8) = -lambda-2*mu
        B(3,8) = -lambda
        B(4,7) = -mu
        B(5,9) = -mu
        B(7,4) = -1./rho
        B(8,2) = -1./rho
        B(9,5) = -1./rho
        !
        ! Jacobian in z-direction
        !
        C(1,9) = -lambda
        C(2,9) = -lambda
        C(3,9) = -lambda-2*mu
        C(5,8) = -mu
        C(6,7) = -mu
        C(7,6) = -1./rho
        C(8,5) = -1./rho
        C(9,3) = -1./rho
        !
        kx = SOURCE%CS%k1(1)
        ky = SOURCE%CS%k1(2)
        kz = SOURCE%CS%k1(3)
        !
        A(:,:) = A(:,:)*( 1. + SOURCE%CS%U0(1)*SIN(kx*xGP+ky*yGP+kz*zGP) )
        B(:,:) = B(:,:)*( 1. + SOURCE%CS%U0(2)*SIN(kx*xGP+ky*yGP+kz*zGP) )
        C(:,:) = C(:,:)*( 1. + SOURCE%CS%U0(3)*SIN(kx*xGP+ky*yGP+kz*zGP) )
        !
        RETURN
    ENDIF

    IF(EQN%Anisotropy.EQ.0.AND.EQN%Poroelasticity.EQ.0) THEN ! Isotropy assumption
        rho    = MaterialVal(1)  ! rho = 1, mu = 2, lambda = 3
        mu     = MaterialVal(2)
        lambda = MaterialVal(3)
        cp     = SQRT((lambda+2*mu)/rho)
        !
        ! Jacobian in x-direction
        !
        A(1,7) = -lambda-2*mu
        A(2,7) = -lambda
        A(3,7) = -lambda
        A(4,8) = -mu
        A(6,9) = -mu
        A(7,1) = -1./rho
        A(8,4) = -1./rho
        A(9,6) = -1./rho
        !
        ! Jacobian in y-direction
        !
        B(1,8) = -lambda
        B(2,8) = -lambda-2*mu
        B(3,8) = -lambda
        B(4,7) = -mu
        B(5,9) = -mu
        B(7,4) = -1./rho
        B(8,2) = -1./rho
        B(9,5) = -1./rho
        !
        ! Jacobian in z-direction
        !
        C(1,9) = -lambda
        C(2,9) = -lambda
        C(3,9) = -lambda-2*mu
        C(5,8) = -mu
        C(6,7) = -mu
        C(7,6) = -1./rho
        C(8,5) = -1./rho
        C(9,3) = -1./rho
        !
    ELSE !IF(EQN%Anisotropy.EQ.0.AND.EQN%Poroelasticity.EQ.0)
        SELECT CASE(LocPoroelastic)
        CASE(0)                  ! Triclinic system
            rho     = MaterialVal( 1)
            ce(1,1) = MaterialVal( 2)
            ce(1,2) = MaterialVal( 3)
            ce(1,3) = MaterialVal( 4)
            ce(1,4) = MaterialVal( 5)
            ce(1,5) = MaterialVal( 6)
            ce(1,6) = MaterialVal( 7)
            ce(2,2) = MaterialVal( 8)
            ce(2,3) = MaterialVal( 9)
            ce(2,4) = MaterialVal(10)
            ce(2,5) = MaterialVal(11)
            ce(2,6) = MaterialVal(12)
            ce(3,3) = MaterialVal(13)
            ce(3,4) = MaterialVal(14)
            ce(3,5) = MaterialVal(15)
            ce(3,6) = MaterialVal(16)
            ce(4,4) = MaterialVal(17)
            ce(4,5) = MaterialVal(18)
            ce(4,6) = MaterialVal(19)
            ce(5,5) = MaterialVal(20)
            ce(5,6) = MaterialVal(21)
            ce(6,6) = MaterialVal(22)
            !
            ! Jacobian in x-direction
            !
            A(1,7:9) = (/ -ce(1,1), -ce(1,6), -ce(1,5) /)
            A(2,7:9) = (/ -ce(1,2), -ce(2,6), -ce(2,5) /)
            A(3,7:9) = (/ -ce(1,3), -ce(3,6), -ce(3,5) /)
            A(4,7:9) = (/ -ce(1,6), -ce(6,6), -ce(5,6) /)
            A(5,7:9) = (/ -ce(1,4), -ce(4,6), -ce(4,5) /)
            A(6,7:9) = (/ -ce(1,5), -ce(5,6), -ce(5,5) /)
            A(7,1) = -1./rho
            A(8,4) = -1./rho
            A(9,6) = -1./rho
            !
            ! Jacobian in y-direction
            !
            B(1,7:9) = (/ -ce(1,6), -ce(1,2), -ce(1,4) /)
            B(2,7:9) = (/ -ce(2,6), -ce(2,2), -ce(2,4) /)
            B(3,7:9) = (/ -ce(3,6), -ce(2,3), -ce(3,4) /)
            B(4,7:9) = (/ -ce(6,6), -ce(2,6), -ce(4,6) /)
            B(5,7:9) = (/ -ce(4,6), -ce(2,4), -ce(4,4) /)
            B(6,7:9) = (/ -ce(5,6), -ce(2,5), -ce(4,5) /)
            B(7,4) = -1./rho
            B(8,2) = -1./rho
            B(9,5) = -1./rho
            !
            ! Jacobian in z-direction
            !
            C(1,7:9) = (/ -ce(1,5), -ce(1,4), -ce(1,3) /)
            C(2,7:9) = (/ -ce(2,5), -ce(2,4), -ce(2,3) /)
            C(3,7:9) = (/ -ce(3,5), -ce(3,4), -ce(3,3) /)
            C(4,7:9) = (/ -ce(5,6), -ce(4,6), -ce(3,6) /)
            C(5,7:9) = (/ -ce(4,5), -ce(4,4), -ce(3,4) /)
            C(6,7:9) = (/ -ce(5,5), -ce(4,5), -ce(3,5) /)
            C(7,6) = -1./rho
            C(8,5) = -1./rho
            C(9,3) = -1./rho
        CASE(1,2,3)                  ! Porous system
            rho_S   = MaterialVal( 1)
            ce(1,1) = MaterialVal( 2)
            ce(1,2) = MaterialVal( 3)
            ce(1,3) = MaterialVal( 4)
            ce(1,4) = MaterialVal( 5)
            ce(1,5) = MaterialVal( 6)
            ce(1,6) = MaterialVal( 7)
            ce(2,2) = MaterialVal( 8)
            ce(2,3) = MaterialVal( 9)
            ce(2,4) = MaterialVal(10)
            ce(2,5) = MaterialVal(11)
            ce(2,6) = MaterialVal(12)
            ce(3,3) = MaterialVal(13)
            ce(3,4) = MaterialVal(14)
            ce(3,5) = MaterialVal(15)
            ce(3,6) = MaterialVal(16)
            ce(4,4) = MaterialVal(17)
            ce(4,5) = MaterialVal(18)
            ce(4,6) = MaterialVal(19)
            ce(5,5) = MaterialVal(20)
            ce(5,6) = MaterialVal(21)
            ce(6,6) = MaterialVal(22)
            rho_F  = MaterialVal(23)
            K_F    = MaterialVal(24)
            nu     = MaterialVal(25)
            K_S    = MaterialVal(26)
            Poro   = MaterialVal(27)
            Kappa(1) = MaterialVal(28)
            Kappa(2) = MaterialVal(29)
            Kappa(3) = MaterialVal(30)
            Tor(1) = MaterialVal(31)
            Tor(2) = MaterialVal(32)
            Tor(3) = MaterialVal(33)
            rho    = rho_S * (1 - Poro) + Poro * rho_F
            K_Mean = 1./9.*(ce(1,1)+ce(2,2)+ce(3,3)+2*(ce(1,2)+ce(1,3)+ce(2,3)))
            MM     = K_S / ((1 - K_Mean/K_S) - Poro*(1 - K_S/K_F))
            Alpha(1) = 1 - (ce(1,1)+ce(1,2)+ce(1,3)) / (3.*K_S)
            Alpha(2) = 1 - (ce(1,2)+ce(2,2)+ce(2,3)) / (3.*K_S)
            Alpha(3) = 1 - (ce(1,3)+ce(2,3)+ce(3,3)) / (3.*K_S)
            Alpha(4) = - (ce(1,4)+ce(2,4)+ce(3,4)) / (3.*K_S)
            Alpha(5) = - (ce(1,5)+ce(2,5)+ce(3,5)) / (3.*K_S)
            Alpha(6) = - (ce(1,6)+ce(2,6)+ce(3,6)) / (3.*K_S)
            Rho1(:)  = rho - (rho_F**2 / (rho_F * Tor(:) / Poro))
            Rho2(:)  = rho_F - (rho_F * Tor(:) / Poro) * rho / rho_F
            Beta1(:) = rho_F / (rho_F * Tor(:) / Poro)
            Beta2(:)  = rho / rho_F
            !
            ! Computation of undrained ce(i,j) coeffs
            DO i=1,6
              DO j=1,6
                ce(i,j) = ce(i,j) + MM * Alpha(i)*Alpha(j)
              ENDDO
            ENDDO
            !
            ! Jacobian in x-direction
            !
            A(1,7:9) = (/ -ce(1,1), -ce(1,6), -ce(1,5) /)
            A(2,7:9) = (/ -ce(1,2), -ce(2,6), -ce(2,5) /)
            A(3,7:9) = (/ -ce(1,3), -ce(3,6), -ce(3,5) /)
            A(4,7:9) = (/ -ce(1,6), -ce(6,6), -ce(5,6) /)
            A(5,7:9) = (/ -ce(1,4), -ce(4,6), -ce(4,5) /)
            A(6,7:9) = (/ -ce(1,5), -ce(5,6), -ce(5,5) /)
            A(7,1) = -1./rho1(1)
            A(8,4) = -1./rho1(1)
            A(9,6) = -1./rho1(1)
            A(10,7:9) = (/ Alpha(1)*MM, Alpha(6)*MM, Alpha(5)*MM /)
            A(11,1) = -1./rho2(1)
            A(12,4) = -1./rho2(1)
            A(13,6) = -1./rho2(1)
            A(7,10) = -beta1(1)/rho1(1)
            A(11,10)= -beta2(1)/rho2(1)
            A(1,11) = -Alpha(1)*MM
            A(2,11) = -Alpha(2)*MM
            A(3,11) = -Alpha(3)*MM
            A(4,11) = -Alpha(6)*MM
            A(5,11) = -Alpha(4)*MM
            A(6,11) = -Alpha(5)*MM
            A(10,11)= MM
            !
            ! Jacobian in y-direction
            !
            B(1,7:9) = (/ -ce(1,6), -ce(1,2), -ce(1,4) /)
            B(2,7:9) = (/ -ce(2,6), -ce(2,2), -ce(2,4) /)
            B(3,7:9) = (/ -ce(3,6), -ce(2,3), -ce(3,4) /)
            B(4,7:9) = (/ -ce(6,6), -ce(2,6), -ce(4,6) /)
            B(5,7:9) = (/ -ce(4,6), -ce(2,4), -ce(4,4) /)
            B(6,7:9) = (/ -ce(5,6), -ce(2,5), -ce(4,5) /)
            B(7,4) = -1./rho1(2)
            B(8,2) = -1./rho1(2)
            B(9,5) = -1./rho1(2)
            B(10,7:9) = (/ Alpha(6)*MM, Alpha(2)*MM, Alpha(4)*MM /)
            B(11,4) = -1./rho2(2)
            B(12,2) = -1./rho2(2)
            B(13,5) = -1./rho2(2)
            B(8,10) = -beta1(2)/rho1(2)
            B(12,10)= -beta2(2)/rho2(2)
            B(1,12) = -Alpha(1)*MM
            B(2,12) = -Alpha(2)*MM
            B(3,12) = -Alpha(3)*MM
            B(4,12) = -Alpha(6)*MM
            B(5,12) = -Alpha(4)*MM
            B(6,12) = -Alpha(5)*MM
            B(10,12)= MM
            !
            ! Jacobian in z-direction
            !
            C(1,7:9) = (/ -ce(1,5), -ce(1,4), -ce(1,3) /)
            C(2,7:9) = (/ -ce(2,5), -ce(2,4), -ce(2,3) /)
            C(3,7:9) = (/ -ce(3,5), -ce(3,4), -ce(3,3) /)
            C(4,7:9) = (/ -ce(5,6), -ce(4,6), -ce(3,6) /)
            C(5,7:9) = (/ -ce(4,5), -ce(4,4), -ce(3,4) /)
            C(6,7:9) = (/ -ce(5,5), -ce(4,5), -ce(3,5) /)
            C(7,6) = -1./rho1(3)
            C(8,5) = -1./rho1(3)
            C(9,3) = -1./rho1(3)
            C(10,7:9) = (/ Alpha(5)*MM, Alpha(4)*MM, Alpha(3)*MM /)
            C(11,6) = -1./rho2(3)
            C(12,5) = -1./rho2(3)
            C(13,3) = -1./rho2(3)
            C(9,10) = -beta1(3)/rho1(3)
            C(13,10)= -beta2(3)/rho2(3)
            C(1,13) = -Alpha(1)*MM
            C(2,13) = -Alpha(2)*MM
            C(3,13) = -Alpha(3)*MM
            C(4,13) = -Alpha(6)*MM
            C(5,13) = -Alpha(4)*MM
            C(6,13) = -Alpha(5)*MM
            C(10,13)= MM
            !
            E(7,11)  = -beta1(1)/rho1(1)*nu/Kappa(1) !Minus sign because E is defined in the right-hand side of hyp. system
            E(8,12)  = -beta1(2)/rho1(2)*nu/Kappa(2)
            E(9,13)  = -beta1(3)/rho1(3)*nu/Kappa(3)
            E(11,11) = -beta2(1)/rho2(1)*nu/Kappa(1)
            E(12,12) = -beta2(2)/rho2(2)*nu/Kappa(2)
            E(13,13) = -beta2(3)/rho2(3)*nu/Kappa(3)
            !
        END SELECT
    ENDIF !IF(EQN%Anisotropy.EQ.0.AND.EQN%Poroelasticity.EQ.0)
    !
    IF(LocAnelastic.EQ.1) THEN
        !
        !Enlargening of the anelastic jacobians (ONLY used in the ADER time integration)
        !
        DO j = 1, EQN%nMechanisms
            A(10+EQN%nAneFuncperMech*(j-1),7)       = -MaterialVal(EQN%AneMatIni+4*(j-1))
            A(13+EQN%nAneFuncperMech*(j-1),8)       = -MaterialVal(EQN%AneMatIni+4*(j-1))/2.
            A(15+EQN%nAneFuncperMech*(j-1),9)       = -MaterialVal(EQN%AneMatIni+4*(j-1))/2.
            B(11+EQN%nAneFuncperMech*(j-1),8)       = -MaterialVal(EQN%AneMatIni+4*(j-1))
            B(13+EQN%nAneFuncperMech*(j-1),7)       = -MaterialVal(EQN%AneMatIni+4*(j-1))/2.
            B(14+EQN%nAneFuncperMech*(j-1),9)       = -MaterialVal(EQN%AneMatIni+4*(j-1))/2.
            C(12+EQN%nAneFuncperMech*(j-1),9)       = -MaterialVal(EQN%AneMatIni+4*(j-1))
            C(14+EQN%nAneFuncperMech*(j-1),8)       = -MaterialVal(EQN%AneMatIni+4*(j-1))/2.
            C(15+EQN%nAneFuncperMech*(j-1),7)       = -MaterialVal(EQN%AneMatIni+4*(j-1))/2.
            E(1,10+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+1)
            E(1,11+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+2)
            E(1,12+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+2)
            E(2,10+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+2)
            E(2,11+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+1)
            E(2,12+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+2)
            E(3,10+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+2)
            E(3,11+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+2)
            E(3,12+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+1)
            E(4,13+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+3)
            E(5,14+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+3)
            E(6,15+EQN%nAneFuncperMech*(j-1))       = -MaterialVal(EQN%AneMatIni+4*(j-1)+3)
            DO k=1,EQN%nAneFuncperMech
                E(9+EQN%nAneFuncperMech*(j-1)+k,9+EQN%nAneFuncperMech*(j-1)+k) = +MaterialVal(EQN%AneMatIni+4*(j-1))
            ENDDO
        ENDDO !j = 1, EQN%nMechanisms
        !
        IF(EQN%Anisotropy.EQ.1) THEN
            ! We have to correct the value for the (4,4), (5,5) and (6,6) entries of the complex modulus
            AvMu     = ( ce(4,4) + ce(5,5) + ce(6,6) ) / 3.d0
            DO j = 1, EQN%nMechanisms
                E(4,13+EQN%nAneFuncperMech*(j-1)) = E(4,13+EQN%nAneFuncperMech*(j-1)) / AvMu * ce(6,6)
                E(5,14+EQN%nAneFuncperMech*(j-1)) = E(5,14+EQN%nAneFuncperMech*(j-1)) / AvMu * ce(4,4)
                E(6,15+EQN%nAneFuncperMech*(j-1)) = E(6,15+EQN%nAneFuncperMech*(j-1)) / AvMu * ce(5,5)
            ENDDO
        ENDIF
        !
    ENDIF !IF(LocAnelastic.EQ.1)

 END SUBROUTINE ModelDefinition_new


  SUBROUTINE SolveRiemannProblem(FL,FR,A,B,C,T,iT,MaterialValLeft,MaterialValRight,Solver,Boundary)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    REAL                            :: FL(:,:)
    REAL                            :: FR(:,:)
    REAL                            :: A(:,:)
    REAL                            :: B(:,:)
    REAL                            :: C(:,:)
    REAL                            :: T(:,:)
    REAL                            :: iT(:,:)
    REAL                            :: MaterialValLeft(:)
    REAL                            :: MaterialValRight(:)
    INTEGER                         :: Solver
    INTEGER                         :: Boundary
    !--------------------------------------------------------------------------
    INTENT(IN)                      :: A,B,C,T,iT,MaterialValLeft,MaterialValRight
    INTENT(IN)                      :: Solver, Boundary
    INTENT(OUT)                     :: FL,FR
    !--------------------------------------------------------------------------
    ! Local variable declaration
    REAL, ALLOCATABLE               :: AbsA(:,:)
    REAL, ALLOCATABLE               :: S(:,:)                                  ! To build free surface boundary
    REAL                            :: rhoL, muL, lambdaL
    REAL                            :: cpL, csL
    REAL                            :: rhoR, muR, lambdaR
    REAL                            :: cpR, csR
    REAL                            :: ConstP, ConstS
    INTEGER                         :: nVar
    INTEGER                         :: iVar
    !--------------------------------------------------------------------------

    ! Left material
    rhoL    = MaterialValLeft(1)
    muL     = MaterialValLeft(2)
    lambdaL = MaterialValLeft(3)
    cpL     = sqrt((lambdaL+2.0d0*muL)/rhoL)
    csL     = sqrt(muL/rhoL)
    ! Right material
    rhoR    = MaterialValRight(1)
    muR     = MaterialValRight(2)
    lambdaR = MaterialValRight(3)
    cpR     = sqrt((lambdaR+2.0d0*muR)/rhoR)
    csR     = sqrt(muR/rhoR)
    !
    nVar = size(A, DIM = 1)
    !
    IF (Boundary.EQ.1) THEN ! Free surface
        ALLOCATE(S(nVar,nVar))
        S = 0.0d0
        S(1,1) = -1.      ! (xx) Normal stress vanishes
        S(2,2) =  1.      ! (yy)
        S(3,3) =  1.      ! (zz)
        S(4,4) = -1.      ! (xy) Shear in normal direction vanishes
        S(5,5) =  1.      ! (yz)
        S(6,6) = -1.      ! (xz) Shear in normal direction vanishes
        S(7,7) =  1.      ! (u)
        S(8,8) =  1.      ! (v)
        S(9,9) =  1.      ! (w)
    ENDIF
    !
    SELECT CASE(Solver)
    CASE(0) ! Godunov flux
        ConstP  = cpL*(lambdaR+2.0d0*muR)+cpR*(lambdaL+2.0d0*muL)
        ConstS  = csL*muR+csR*muL

        FR(:,:) = 0.0d0;
        FR(1,1) = cpR * (lambdaL+2.0d0*muL) / ConstP
        FR(2,1) = cpR *  lambdaL            / ConstP
        FR(3,1) = FR(2,1)
        FR(7,1) = cpR *  cpL                / ConstP
        !
        IF (ConstS.EQ.0.0d0) THEN
            CONTINUE
        ELSE
            FR(4,4) = csR *  muL            / ConstS
            FR(8,4) = csR *  csL            / ConstS
            !
            FR(6,6) = FR(4,4)
            FR(9,6) = FR(8,4)
        ENDIF
        !
        FR(1,7) = (lambdaR+2.0d0*muR) * (lambdaL+2.0d0*muL) / ConstP
        FR(2,7) = (lambdaR+2.0d0*muR) *  lambdaL            / ConstP
        FR(3,7) = FR(2,7)
        FR(7,7) = (lambdaR+2.0d0*muR) *  cpL                / ConstP
        !
        IF (ConstS.EQ.0.0d0) THEN
            CONTINUE
        ELSE
            FR(4,8) = muR *  muL            / ConstS
            FR(8,8) = muR *  csL            / ConstS
            !
            FR(6,9) = FR(4,8)
            FR(9,9) = FR(8,8)
        ENDIF

        FL(:,:) = -FR(:,:)
        DO iVar=1,9
            IF (muR.EQ.0.0d0 .AND. (iVar.GE.4 .AND. iVar.LE.6)) THEN
                CONTINUE
            ELSE
                FL(iVar,iVar) = 1.0d0+FL(iVar,iVar)
            ENDIF
        ENDDO ! iVar
        !
        FL = MATMUL(T,MATMUL(MATMUL(A,FL),iT))
        !
        IF (Boundary.EQ.1) FR = MATMUL(FR,S) ! Free surface
        !
        FR = MATMUL(T,MATMUL(MATMUL(A,FR),iT))
    CASE(1) ! Rusanov flux (|A|=diag(max(alpha_i)))
        ALLOCATE(AbsA(nVar,nVar))
        absA = 0.0d0
        DO iVar = 1,nVar
            absA(iVar,iVar) = cpL
        ENDDO
        FL = 0.5d0*MATMUL(T,MATMUL(A+AbsA,iT))
        FR = 0.5d0*MATMUL(T,MATMUL(A-AbsA,iT))
        DEALLOCATE(AbsA)
    CASE DEFAULT
        print*,' Error in Riemann Solver subroutine. Wrong solver'
        STOP
    END SELECT
  END SUBROUTINE SolveRiemannProblem

  SUBROUTINE prem_iso(r_vert,r,nVar,MaterialVal,Ani)
!   subroutine prem_iso(rho,drhodr,mu,lambda,Qkappa,Qmu)

  IMPLICIT NONE

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  REAL              :: MaterialVal(nVar)
  REAL              :: r_vert, r, rho, drhodr, vp, vs, Qkappa, Qmu, x,&
  RICB, RCMB, RTOPDDOUBLEPRIME, R600, R670, R220, R771, R400, &
  R80, RMOHO, RMIDDLE_CRUST,ROCEAN
  REAL              :: scaleval, mu, lambda
  REAL, PARAMETER   :: GRAV = 6.6723e-11
  REAL, PARAMETER   :: RHOAV = 5514.30
  REAL, PARAMETER   :: R_EARTH = 6371000.0
  INTEGER           :: nVar,Ani,io_error,i
  LOGICAL, PARAMETER :: NOCRUST = .false.

  ROCEAN = 6368000.0
  RMIDDLE_CRUST = 6356000.0
  RMOHO = 6346600.0
  R80  = 6291000.0
  R220 = 6151000.0
  R400 = 5971000.0
  R600 = 5771000.0
  R670 = 5701000.0
  R771 = 5600000.0
  RTOPDDOUBLEPRIME = 3630000.0
  RCMB = 3480000.0
  RICB = 1221000.0

  x = r_vert/R_EARTH

!  DO i=1,6371
!     x=(REAL(i)*1000.0)/R_EARTH
!     r=(REAL(i)*1000.0)

!
!--- inner core
!
  if(r >= 0.0 .and. r <= RICB) then
    drhodr=-2.00*8.83810*x
    rho=13.08850-8.83810*x*x
    vp=11.26220-6.36400*x*x
    vs=3.66780-4.44750*x*x
    Qmu=84.60
    Qkappa=1327.70
!
!--- outer core
!
  else if(r > RICB .and. r <= RCMB) then
    drhodr=-1.26380-2.00*3.64260*x-3.00*5.52810*x*x
    rho=12.58150-1.26380*x-3.64260*x*x-5.52810*x*x*x
    vp=11.04870-4.03620*x+4.80230*x*x-13.57320*x*x*x
    vs=1.0e-10
    Qmu=9999999.9
    Qkappa=57827.00
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    drhodr=-6.47610+2.00*5.52830*x-3.00*3.08070*x*x
    rho=7.95650-6.47610*x+5.52830*x*x-3.08070*x*x*x
    vp=15.38910-5.31810*x+5.52420*x*x-2.55140*x*x*x
    vs=6.92540+1.46720*x-2.08340*x*x+0.97830*x*x*x
    Qmu=312.00
    Qkappa=57827.00
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    drhodr=-6.47610+2.00*5.52830*x-3.00*3.08070*x*x
    rho=7.95650-6.47610*x+5.52830*x*x-3.08070*x*x*x
    vp=24.95200-40.46730*x+51.48320*x*x-26.64190*x*x*x
    vs=11.16710-13.78180*x+17.45750*x*x-9.27770*x*x*x
    Qmu=312.00
    Qkappa=57827.00
!
  else if(r > R771 .and. r <= R670) then
    drhodr=-6.47610+2.00*5.52830*x-3.00*3.08070*x*x
    rho=7.95650-6.47610*x+5.52830*x*x-3.08070*x*x*x
    vp=29.27660-23.60270*x+5.52420*x*x-2.55140*x*x*x
    vs=22.34590-17.24730*x-2.08340*x*x+0.97830*x*x*x
    Qmu=312.00
    Qkappa=57827.00
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R600) then
    drhodr=-1.48360
    rho=5.31970-1.48360*x
    vp=19.09570-9.86720*x
    vs=9.98390-4.93240*x
    Qmu=143.00
    Qkappa=57827.00
  else if(r > R600 .and. r <= R400) then
    drhodr=-8.02980
    rho=11.24940-8.02980*x
    vp=39.70270-32.61660*x
    vs=22.35120-18.58560*x
    Qmu=143.00
    Qkappa=57827.00
  else if(r > R400 .and. r <= R220) then
    drhodr=-3.80450
    rho=7.10890-3.80450*x
    vp=20.39260-12.25690*x
    vs=8.94960-4.45970*x
    Qmu=143.00
    Qkappa=57827.00
  else if(r > R220 .and. r <= R80) then
    drhodr=0.69240
    rho=2.69100+0.69240*x
    vp=4.18750+3.93820*x
    vs=2.15190+2.34810*x
    Qmu=80.00
    Qkappa=57827.00
  else
    if(NOCRUST) then
      drhodr=0.69240
      rho=2.69100+0.69240*(RMOHO / R_EARTH)
      vp=4.18750+3.93820*(RMOHO / R_EARTH)
      vs=2.15190+2.34810*(RMOHO / R_EARTH)
      Qmu=600.00
      Qkappa=57827.00
    else if (r > r80 .and. r <= RMOHO) then
      drhodr=0.69240
      rho=2.69100+0.69240*x
      vp=4.18750+3.93820*x
      vs=2.15190+2.34810*x
      Qmu=600.00
      Qkappa=57827.00
    else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
      drhodr=0.00
      rho=2.90
      vp=6.80
      vs=3.90
      Qmu=600.00
      Qkappa=57827.00
    else
      drhodr=0.00
      rho=2.60
      vp=5.80
      vs=3.20
      Qmu=600.00
      Qkappa=57827.00
    endif
  endif

  drhodr = drhodr*1000.00
  rho = rho*1000.00
  vp = vp*1000.00
  vs = vs*1000.00

  mu = (vs*vs)*rho
  lambda = ((vp*vp)*rho)-2.0*mu

!  open(unit=20,file='model.dat',status='old',action='write',position='append',iostat=io_error)
!  write(20,'(6(E15.8))') r,rho,mu,lambda,vp,vs
!  close(unit=20)

!  END DO

  MaterialVal(1) = rho
  MaterialVal(2) = mu
  MaterialVal(3) = lambda

  IF( Ani .EQ. 1 ) THEN
    MaterialVal(4) = Qkappa
    MaterialVal(5) = Qmu
  ENDIF

  END SUBROUTINE prem_iso


  subroutine svbacksb(wsvd,vsvd,usvd,b,x)
    !perform svd back substitution
    real :: usvd(:,:),vsvd(:,:),wsvd(:),b(:),x(:),t(size(x))
    t=0.0; where(abs(wsvd)>=1e-10)t=matmul(b,usvd)/wsvd; x=matmul(vsvd,t)
  end subroutine svbacksb

  subroutine svdecomp(wsvd,vsvd,usvd)
    !perform singular velue decomposition
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    real :: usvd(:,:),wsvd(:),vsvd(:,:)
    integer :: non,i,j,k,n1,n2,i1,i2,it,ic
    real    :: fac,f1,f2,f3,vc,vs,v1,v2,v3,dist
    real    :: n1pt(size(usvd,1)),svdrv(size(usvd,2)),n2pt(size(usvd,2))
    real,   parameter :: l1=1.0e-12,l2=1.0e-14,zero=0.0,one=1.0
    integer,parameter :: l3=41
    character(len=20), parameter :: message='svd: no convergence'
    f2=zero;fac=zero;
    n1=size(usvd,1);n2=size(usvd,2)
    do i=1,n2
       ic=i+1; svdrv(i)=fac*f2; f2=zero;fac=zero
       if (i <= n1) then
          fac=ar(usvd(i:n1,i)); if(abs(fac)>l1)i2=isub1(i1)
       end if
       wsvd(i) = fac*f2;f2=zero; fac=zero
       if ((i /= n2).and.(i <= n1)) then
          fac = ar(usvd(i,ic:n2)); if(abs(fac)>l1)i1=isub2(i2)
       end if
    end do
    dist=maxval(abs(wsvd)+abs(svdrv))
    do i=n2,1,-1
       if (i < n2) then
          if (abs(f2)>l1) i2=isub3(i1)
          vsvd(i,ic:n2)=zero;vsvd(ic:n2,i)=zero
       end if
       vsvd(i,i)=one; f2=svdrv(i); ic=i
    end do
    do i=min(n1,n2),1,-1
       i2=isub4(i)
    end do
    do k=n2,1,-1
       do it=1,l3
          do ic=k,1,-1
             non=ic-1
             if ((abs(svdrv(ic))+dist) == dist) exit
             if ((abs(wsvd(non))+dist) == dist) then
                i1=isub0(i2);exit
             end if
          end do
          v3=wsvd(k)
          if (ic == k) then
             if (v3 < zero) then
                wsvd(k)=-v3;vsvd(1:n2,k)=-vsvd(1:n2,k)
             end if
             exit
          end if
          if (it==l3) write(error_unit,*)message
          i2=isub5(i1);i1=isub6(i2);svdrv(ic)=zero;
          svdrv(k)=f1;wsvd(k)=v1
       end do
    end do
  contains
    integer function isub0(in)
      integer :: in
      vc=zero; vs=one;isub0=in;
      do i=ic,k
         f1=vs*svdrv(i); svdrv(i)=vc*svdrv(i)
         if ((abs(f1)+dist) == dist) exit
         f2=wsvd(i); f3=hypot(f1,f2)
         wsvd(i)=f3; f3=one/f3; vc= (f2*f3); vs=-(f1*f3)
         n1pt(1:n1)=usvd(1:n1,non)
         usvd(1:n1,non)=usvd(1:n1,non)*vc+usvd(1:n1,i)*vs
         usvd(1:n1,i)=-n1pt(1:n1)*vs+usvd(1:n1,i)*vc
      end do
    end function isub0
    integer function isub1(in)
      integer :: in
      isub1=in; usvd(i:n1,i)=usvd(i:n1,i)/fac
      vs=dp(usvd(i:n1,i),usvd(i:n1,i))
      f1=usvd(i,i);f2=-ss(vs,f1); f3=f1*f2-vs; usvd(i,i)=f1-f2
      n2pt(ic:n2)=matmul(usvd(i:n1,i),usvd(i:n1,ic:n2))/f3
      usvd(i:n1,ic:n2) = usvd(i:n1,ic:n2) + &
           &      spread(usvd(i:n1,i),dim=2,ncopies=size(n2pt(ic:n2))) * &
           &      spread(n2pt(ic:n2),dim=1,ncopies=size(usvd(i:n1,i)))
      usvd(i:n1,i)=fac*usvd(i:n1,i)
    end function isub1
    integer function isub2(in)
      integer :: in
      isub2=in; usvd(i,ic:n2)=usvd(i,ic:n2)/fac
      vs=dp(usvd(i,ic:n2),usvd(i,ic:n2))
      f1=usvd(i,ic);f2=-ss(vs,f1); f3=f1*f2-vs
      usvd(i,ic)=f1-f2; svdrv(ic:n2)=usvd(i,ic:n2)/f3
      n1pt(ic:n1)=matmul(usvd(ic:n1,ic:n2),usvd(i,ic:n2))
      usvd(ic:n1,ic:n2) = usvd(ic:n1,ic:n2) + &
           & spread(n1pt(ic:n1),dim=2,ncopies=size(svdrv(ic:n2))) * &
           & spread(svdrv(ic:n2),dim=1,ncopies=size(n1pt(ic:n1)))
      usvd(i,ic:n2)=fac*usvd(i,ic:n2)
    end function isub2
    integer function isub3(in)
      integer :: in
      isub3=in;vsvd(ic:n2,i)=(usvd(i,ic:n2)/usvd(i,ic))/f2
      n2pt(ic:n2)=matmul(usvd(i,ic:n2),vsvd(ic:n2,ic:n2))
      vsvd(ic:n2,ic:n2) = vsvd(ic:n2,ic:n2) +  &
           & spread(vsvd(ic:n2,i),dim=2,ncopies=size(n2pt(ic:n2))) * &
           & spread(n2pt(ic:n2),dim=1,ncopies=size(vsvd(ic:n2,i)))
    end function isub3
    integer function isub4(in)
      integer :: in
      ic=i+1; f2=wsvd(i);isub4=in
      usvd(i,ic:n2)=zero
      if (abs(f2)>l2) then
         f2=one/f2
         n2pt(ic:n2)=(matmul(usvd(ic:n1,i),usvd(ic:n1,ic:n2))/usvd(i,i))*f2
         usvd(i:n1,ic:n2) = usvd(i:n1,ic:n2)  &
              & + spread(usvd(i:n1,i),dim=2,ncopies=size(n2pt(ic:n2))) * &
              &   spread(n2pt(ic:n2),dim=1,ncopies=size(usvd(i:n1,i)))
         usvd(i:n1,i)=usvd(i:n1,i)*f2
      else
         usvd(i:n1,i)=zero
      end if
      usvd(i,i)=usvd(i,i)+one
    end function isub4
   real function dp(x,y)
     real :: x(:),y(:)
     dp = dot_product(x,y)
   end function dp
   integer function isub5(in)
     integer :: in
     isub5=in;v1=wsvd(ic); non=k-1;v2=wsvd(non)
     f2=svdrv(non); f3=svdrv(k)
     f1=((v2-v3)*(v2+v3)+(f2-f3)*(f2+f3))/(2.0*f3*v2)
     f2=hypot(f1,one)
     f1=((v1-v3)*(v1+v3)+f3*((v2/(f1+sign(f2,f1)))-f3))/v1
     vc=one; vs=one
   end function isub5
   integer function isub6(in)
     integer :: in
     do j=ic,non
        i=j+1; v2=wsvd(i)
        f2=svdrv(i); f3=vs*f2; f2=vc*f2
        v3=hypot(f1,f3)
        svdrv(j)=v3; vc=f1/v3; vs=f3/v3
        f1= (v1*vc)+(f2*vs);f2=-(v1*vs)+(f2*vc);f3=v2*vs
        v2=v2*vc;n2pt(1:n2)=vsvd(1:n2,j)
        vsvd(1:n2,j)=vsvd(1:n2,j)*vc+vsvd(1:n2,i)*vs
        vsvd(1:n2,i)=-n2pt(1:n2)*vs+vsvd(1:n2,i)*vc
        v3=hypot(f1,f3);wsvd(j)=v3;isub6=in
        if (abs(v3)>l1) then
           v3=one/v3;vc=f1*v3;vs=f3*v3
        end if
        f1= (vc*f2)+(vs*v2);v1=-(vs*f2)+(vc*v2)
        n1pt(1:n1)=usvd(1:n1,j)
        usvd(1:n1,j)=usvd(1:n1,j)*vc+usvd(1:n1,i)*vs
        usvd(1:n1,i)=-n1pt(1:n1)*vs+usvd(1:n1,i)*vc
     end do
   end function isub6
   elemental real function ss(x,y)
     real , intent(in) :: x,y
     ss = sign(sqrt(x),y)
   end function ss
   real function ar(x)
     real :: x(:)
     ar=sum(abs(x))
   end function ar
 end subroutine svdecomp


END MODULE ini_MODEL_mod
