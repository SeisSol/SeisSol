!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Josep de la Puente Alvarez (josep.delapuente AT bsc.es, http://www.geophysik.uni-muenchen.de/Members/jdelapuente)
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

MODULE DGSponge_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  INTERFACE ini_DGSponge
     MODULE PROCEDURE ini_DGSponge
  END INTERFACE
  INTERFACE DGSponge
     MODULE PROCEDURE DGSponge
  END INTERFACE
  ! 
  !---------------------------------------------------------------------------!
  PUBLIC  :: ini_DGSponge, DGSponge
  !---------------------------------------------------------------------------!

CONTAINS


  SUBROUTINE ini_DGSponge(EQN,DISC,MESH,IO)
    USE COMMON_operators_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tDiscretization)           :: DISC
    TYPE(tInputOutput)              :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i,j,k,iElem, iSide, iFacet, iIntGP
    INTEGER                         :: iSponge, iObject, iBNDType
    INTEGER                         :: iBNDElem, iBNDSide
    INTEGER                         :: counter, nFacet
    INTEGER                         :: nearestFacet
    INTEGER                         :: SpongeList(MESH%nElem)
    INTEGER                         :: SpongeFacet(MESH%nElem,2)
    INTEGER                         :: VertexSide_Tri(3,EQN%DIMENSION)
    INTEGER                         :: VertexSide_Quad(4,EQN%DIMENSION)
    INTEGER                         :: VertexSide(4,EQN%DIMENSION)
    INTEGER                         :: NMAX
    INTEGER, POINTER                :: FacetList(:,:)
    REAL                            :: xP,yP,zP,dr,dphi,r1,r2, r,phi,state(5),time
    REAL                            :: distance, distance_new
    REAL                            :: Factor
    REAL                            :: x(MESH%nVertexMax) ! Element vertices in x-y space         !
    REAL                            :: y(MESH%nVertexMax) ! Element vertices in x-y space         !
    REAL                            :: z(MESH%nVertexMax) ! Element vertices in x-y space         !
    REAL                            :: xGP,yGP,zGP
    REAL                            :: SpongeDistance(MESH%nElem)
    REAL                            :: Node1(EQN%DIMENSION), Node2(EQN%DIMENSION), Node3(EQN%DIMENSION)
    REAL                            :: Vec1(EQN%DIMENSION), Vec2(EQN%DIMENSION), Normal(EQN%DIMENSION)
    REAL, POINTER                   :: FacetBary(:,:)
    ! ------------------------------------------------------------------------!
    INTENT(IN)    :: EQN, MESH, IO
    INTENT(INOUT) :: DISC
    ! ------------------------------------------------------------------------!
    !
    !                                                                         !

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
    ! Per default, no element is in a sponge yet
    SpongeList(:)     = -1
    SpongeDistance(:) = 1e20
    SpongeFacet(:,:)  = -1
    !
    DO iSponge = 1, DISC%Galerkin%nDGSponge
       !
       iObject  = DISC%Galerkin%DGSponge(iSponge)%SpongeObject
       iBNDType = DISC%Galerkin%DGSponge(iSponge)%SpongeBND
       !
       ! Identify boundary facets for this sponge
       !
       counter = 0
       DO iElem = 1, MESH%nElem
          DO iSide = 1, MESH%LocalElemType(iElem)
             IF(MESH%ELEM%Reference(iSide,iElem).EQ.iBNDType) THEN
               IF(MESH%ELEM%BoundaryToObject(iSide,iElem).EQ.iObject) THEN
                 counter = counter + 1
               ENDIF
             ENDIF
          ENDDO
       ENDDO
       nFacet = counter
       !
       logInfo(*) 'Sponge ', iSponge, ' has ', nFacet, ' boundary facets. '
       !
       ALLOCATE(FacetList(nFacet,2))
       ALLOCATE(FacetBary(nFacet,EQN%DIMENSION))
       counter = 0
       DO iElem = 1, MESH%nElem
          DO iSide = 1, MESH%LocalElemType(iElem)
             IF(MESH%ELEM%Reference(iSide,iElem).EQ.iBNDType) THEN
               IF(MESH%ELEM%BoundaryToObject(iSide,iElem).EQ.iObject) THEN
                  counter = counter + 1
                  FacetList(counter,1) = iElem
                  FacetList(counter,2) = iSide
                  SELECT CASE(MESH%LocalElemType(iElem))
                  CASE(3)
                    Node1(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tri(iSide,1),iElem))
                    Node2(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tri(iSide,2),iElem))
                  CASE(4)
                    Node1(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Quad(iSide,1),iElem))
                    Node2(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Quad(iSide,2),iElem))
                  END SELECT
                    Node3(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iSide,3),iElem))
                    FacetBary(counter,:) = 1./3.*( Node1(:)+Node2(:)+Node3(:) )
               ENDIF
             ENDIF
          ENDDO
       ENDDO
       !
       ! Identify the elements affected by this sponge
       !
       Factor = 1.+DISC%Galerkin%DGSpongeTol
       !
       DO iElem = 1, MESH%nElem
          ! If barycenter lies within a certain tolerance w.r.t. the sponge
          ! thickness, consider this element for sponge action
          Node1(:) = MESH%ELEM%xyBary(:,iElem)
          distance = 1e20
          nearestFacet = -1
          DO iFacet = 1, nFacet
            Node2(:) = FacetBary(iFacet,:)
            distance_new = SUM( (Node1(:)-Node2(:))**2 ) 
            IF(distance_new.LT.distance) THEN
               distance = distance_new
               nearestFacet = iFacet
            ENDIF
          ENDDO
          IF(SQRT(distance).LE.(Factor*DISC%Galerkin%DGSponge(iSponge)%SpongeDelta)) THEN
             IF(SQRT(distance).LE.SpongeDistance(iElem)) THEN
               SpongeList(iElem)     = iSponge
               SpongeDistance(iElem) = SQRT(distance)
               IF(nearestFacet.EQ.-1) THEN
                  logError(*) 'Facet Error in DGSponge.'
                  STOP
               ENDIF
               SpongeFacet(iElem,:)  = FacetList(nearestFacet,:)
             ENDIF
          ENDIF
          !
          IF(MESH%nElem.GT.20) THEN
             IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
                 logInfo(*) iElem,' elements done...'
             END IF
          ENDIF
       ENDDO
       !
       DEALLOCATE(FacetList, FacetBary)
       !   
    ENDDO
    !
    DO iSponge = 1, DISC%Galerkin%nDGSponge
       !
       counter = 0
       DO iElem = 1, MESH%nElem
         IF(SpongeList(iElem).EQ.iSponge) THEN
           counter = counter + 1
         ENDIF
       ENDDO
       DISC%Galerkin%DGSponge(iSponge)%nSpongeElements = counter
       ALLOCATE(DISC%Galerkin%DGSponge(iSponge)%SpongeElements(counter) )
       ALLOCATE(DISC%Galerkin%DGSponge(iSponge)%SpongeDistances(counter,DISC%Galerkin%nIntGP) )
       logInfo(*) 'Sponge ', iSponge, ' affects ', counter, ' Elements.'
       !
       counter = 0
       DO iElem = 1, MESH%nElem
         IF(SpongeList(iElem).EQ.iSponge) THEN
           counter = counter + 1
           DISC%Galerkin%DGSponge(iSponge)%SpongeElements(counter) = iElem
         ENDIF
       ENDDO
       !
       ! Alle Gausspunkte in Ebenengleichung des n�chsten Facets einsetzen und damit
       ! (approximativen) Abstand ausrechnen
       ! 
       DO i = 1, DISC%Galerkin%DGSponge(iSponge)%nSpongeElements

          iElem    = DISC%Galerkin%DGSponge(iSponge)%SpongeElements(i)
          DO j=1,MESH%LocalElemType(iElem)
             x(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(j,iElem))
             y(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(j,iElem))
             z(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(j,iElem))
          ENDDO
          ! Nearest facet
          iBNDElem = SpongeFacet(iElem,1)
          iBNDSide = SpongeFacet(iElem,2)
          ! Node coordinates of the nearest facet
          SELECT CASE(MESH%LocalElemType(iElem))
          CASE(3)
                Node1(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tri(iBNDSide,1),iBNDElem))
                Node2(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tri(iBNDSide,2),iBNDElem))
          CASE(4)
                Node1(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Quad(iBNDSide,1),iBNDElem))
                Node2(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Quad(iBNDSide,2),iBNDElem))
          END SELECT
            Node3(:) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide(iBNDSide,3),iBNDElem))
          ! Compute normal vector of facet
          Vec1(:) = Node2(:)-Node1(:)
          Vec2(:) = Node3(:)-Node1(:)
          Normal(:) = Vec1(:).x.Vec2(:)
          Normal(:) = Normal(:)/SQRT( SUM(Normal(:)**2) )    
          ! Insert Gausspoints into Hessian form to compute distance           
          ! distance = (x-x0)*n (dotproduct)
          DO iIntGP = 1, DISC%Galerkin%nIntGP
             !
                 xGP = x(1) + (x(2)-x(1))*DISC%Galerkin%intGaussP(1,iIntGP) +  &
                              (x(3)-x(1))*DISC%Galerkin%intGaussP(2,iIntGP) +  &
                              (x(4)-x(1))*DISC%Galerkin%intGaussP(3,iIntGP) 
                 !
                 yGP = y(1) + (y(2)-y(1))*DISC%Galerkin%intGaussP(1,iIntGP) +  &
                              (y(3)-y(1))*DISC%Galerkin%intGaussP(2,iIntGP) +  &
                              (y(4)-y(1))*DISC%Galerkin%intGaussP(3,iIntGP) 
                 !
                 zGP = z(1) + (z(2)-z(1))*DISC%Galerkin%intGaussP(1,iIntGP) +  &
                              (z(3)-z(1))*DISC%Galerkin%intGaussP(2,iIntGP) +  &
                              (z(4)-z(1))*DISC%Galerkin%intGaussP(3,iIntGP) 
                 ! Vector of the Gaussian point
                 vec1(:) = (/ xGP, yGP, zGP /)
             ! Difference vector from the facet plane (we take node 1)
             vec2(:) = vec1(:)-Node1(:)
             distance = DOT_PRODUCT(vec2(:),Normal(:))
             DISC%Galerkin%DGSponge(iSponge)%SpongeDistances(i,iIntGP) = ABS(distance)
          ENDDO
          !
       ENDDO
       !
     ENDDO
     !
     !
  END SUBROUTINE ini_DGSponge


  SUBROUTINE DGSponge(time, EQN, DISC, MESH, BND, IC, SOURCE, IO)  
    !-------------------------------------------------------------------------!
    USE DGBasis_mod
    USE COMMON_InitialField_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    TYPE(tBoundary)          :: BND
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tInputOutput)       :: IO
    TYPE(tInitialCondition)  :: IC
    TYPE(tSource)            :: SOURCE
    REAL                     :: time
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER :: iElem                  ! Element number                        !
    INTEGER :: iSponge                ! Sponge number                         !
    INTEGER :: iType                  ! Sponge type (4)=inflow, (5)=outflow   ! 
    INTEGER :: iObject                ! Sponge object number                  ! 
    INTEGER :: iIntGP                 ! Index of internal Gausspoint          !
    INTEGER :: iDegFr                 ! Degree of freedom                     !
    INTEGER :: iVar                   ! Variable number                       ! 
    INTEGER :: i,j                    ! Loop variables                        !
    REAL    :: xGP, yGP, zGP          ! Coordinates of GP in x-y-z            !
    REAL    :: x(MESH%nVertexMax)     ! Element vertices in x-y space         !
    REAL    :: y(MESH%nVertexMax)     ! Element vertices in x-y space         !
    REAL    :: z(MESH%nVertexMax)     ! Element vertices in x-y space         !
    REAL    :: state_ref(EQN%nVar)    ! Reference state     (ref)
    REAL    :: state_ane(3*EQN%nMechanisms)  ! Initial anelastic state vector at Gausspoint    
    REAL    :: state_raw(EQN%nVar)    ! State before sponge (raw)
    REAL    :: state_spo(EQN%nVar)    ! State after sponge  (spo)
    REAL    :: SpongeDOF(DISC%Galerkin%nDegFr,EQN%nVar)
    REAL    :: sigma
    REAL    :: distance
    REAL    :: phi
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: EQN, MESH, IO
    INTENT(INOUT) :: DISC
    !-------------------------------------------------------------------------!
    !
    DO iSponge = 1, DISC%Galerkin%nDGSponge
       !
       iType   = DISC%Galerkin%DGSponge(iSponge)%SpongeBND
       iObject = DISC%Galerkin%DGSponge(iSponge)%SpongeObject
       !
       DO i = 1, DISC%Galerkin%DGSponge(iSponge)%nSpongeElements
          !
          iElem    = DISC%Galerkin%DGSponge(iSponge)%SpongeElements(i)
          DO j=1,MESH%LocalElemType(iElem)
             x(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(j,iElem))
             y(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(j,iElem))
             z(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(j,iElem))
          ENDDO

          SpongeDOF(:,:) = 0.
          DO iIntGP = 1,DISC%Galerkin%nIntGP
                !
                !
                    xGP = x(1) + (x(2)-x(1))*DISC%Galerkin%intGaussP(1,iIntGP) +  &
                                 (x(3)-x(1))*DISC%Galerkin%intGaussP(2,iIntGP) +  &
                                 (x(4)-x(1))*DISC%Galerkin%intGaussP(3,iIntGP) 
                    !
                    yGP = y(1) + (y(2)-y(1))*DISC%Galerkin%intGaussP(1,iIntGP) +  &
                                 (y(3)-y(1))*DISC%Galerkin%intGaussP(2,iIntGP) +  &
                                 (y(4)-y(1))*DISC%Galerkin%intGaussP(3,iIntGP) 
                    !
                    zGP = z(1) + (z(2)-z(1))*DISC%Galerkin%intGaussP(1,iIntGP) +  &
                                 (z(3)-z(1))*DISC%Galerkin%intGaussP(2,iIntGP) +  &
                                 (z(4)-z(1))*DISC%Galerkin%intGaussP(3,iIntGP) 
                !
                !
                distance = DISC%Galerkin%DGSponge(iSponge)%SpongeDistances(i,iIntGP)
                IF(distance.LE.DISC%Galerkin%DGSponge(iSponge)%SpongeDelta) THEN
                    sigma =   DISC%Galerkin%DGSponge(iSponge)%SigmaMax *                &
                           ( (DISC%Galerkin%DGSponge(iSponge)%SpongeDelta-distance) /   &
                             (DISC%Galerkin%DGSponge(iSponge)%SpongeDelta         ) )** &
                              DISC%Galerkin%DGSponge(iSponge)%SpongePower
                ELSE
                    sigma = 0. ! No sponge
                ENDIF
                !
                !CALL GetStateGP(state_raw,iElem,iIntGP,EQN,DISC)
                !
                 CALL GetStateXiEta(                                &
                       state  = state_raw,                          &
                       xi     = MESH%ELEM%GP_Tri(1,iIntGP),         &
                       eta    = MESH%ELEM%GP_Tri(2,iIntGP),         &
                       nDegFr = DISC%Galerkin%nDegFr,               &
                       nVar   = EQN%nVar,                           &
                       DOF    = DISC%Galerkin%DGvar(:,:,iElem,1),   &
                       EQN    = EQN,                                &
                       DISC   = DISC,                               & 
                       eType  = MESH%LocalElemType(iElem)           )   


                IF(iType.EQ.4) THEN  ! Inflow
                   SELECT CASE(BND%ObjInflow(iObject)%InflowType)
                   CASE(0)    ! Constant data
                     state_ref(:) = BND%ObjInflow(i)%u0_in(:)
                   CASE(1,2)
                     CALL InitialField(state_ref(:), state_ane(:), time, xGP, yGP, zGP, iElem, EQN, IC, SOURCE, IO)                     
                   CASE DEFAULT
                        PRINT *, 'Error in DGSponge. Inflow boundary type not known.'
                        STOP
                   END SELECT
                ELSE 
                   state_ref(:) = 0.
                ENDIF
                !
                state_spo(:) = state_raw(:)*(1.-sigma) + sigma*state_ref(:)
                !
                !
                DO iDegFr = 1, DISC%Galerkin%nDegFr
                   !
                   ! Projection onto the basis functions
                   phi = MESH%ELEM%BF_GP_Tri(iDegFr,iIntGP)
                   !
                   SpongeDOF(iDegFr,:) = SpongeDOF(iDegFr,:) + &
                                         MESH%ELEM%GW_Tri(iIntGP)*state_spo(:)*phi
                   !
                ENDDO
                !
          ENDDO
          ! 
          DO iDegFr = 1,DISC%Galerkin%nDegFr
             SpongeDOF(iDegFr,:) = SpongeDOF(iDegFr,:)*DISC%Galerkin%iMassMatrix_Tri(iDegFr,iDegFr,DISC%Galerkin%nPoly)
          ENDDO
          !
          ! Update DOF 
          !
          DISC%Galerkin%dgvar( :,1:EQN%nVar,iElem,1) = SpongeDOF(:,1:EQN%nVar)
          IF(DISC%Galerkin%DGMethod.NE.3) THEN
            DISC%Galerkin%dgwork(:,1:EQN%nVar,iElem  ) = SpongeDOF(:,1:EQN%nVar)
          ENDIF
          !
        ENDDO
        !
    ENDDO
    !
  END SUBROUTINE DGSponge


END MODULE DGSponge_mod

