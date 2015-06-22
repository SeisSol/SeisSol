!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2011, SeisSol Group
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

MODULE Galerkin_source_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  INTERFACE GalerkinSource3D
     MODULE PROCEDURE GalerkinSource3D
  END INTERFACE
  !
  !---------------------------------------------------------------------------!
  PUBLIC  :: GalerkinSource3D
  !---------------------------------------------------------------------------!

CONTAINS


  SUBROUTINE GalerkinSource3D(src,iElem,time,dt,LocDegFr,OptionalFields,EQN,DISC,MESH,SOURCE,IO,IC,src_deriv)
    !-------------------------------------------------------------------------!
    USE DGBasis_mod
    USE COMMON_operators_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tSource)            :: SOURCE
    TYPE(tInitialCondition),OPTIONAL  :: IC
    TYPE(tInputOutput)       :: IO
    TYPE (tUnstructOptionalFields) :: OptionalFields
    INTEGER                        :: iElem, LocDegFr
    INTEGER                        :: i_left, i_right
    REAL                           :: time
    REAL                           :: dt
    REAL                           :: src(DISC%Galerkin%nDegFr,EQN%nVar)
    REAL,OPTIONAL                  :: src_deriv(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)
    ! Local variable declaration
    REAL                           :: x(MESH%GlobalVrtxType)
    REAL                           :: y(MESH%GlobalVrtxType)
    REAL                           :: z(MESH%GlobalVrtxType)
    REAL                           :: phi
    REAL                           :: xGP, yGP, zGP
    REAL                           :: xi, eta, zeta
    REAL                           :: strk, dip, rake
    REAL                           :: f, integral, M0_integral, tau, t, value
    REAL                           :: a1, a2, m, hw2,fac
    REAL                           :: LocalMomentTensor(3,3)
    REAL                           :: t_proj, t_left, d_t
    REAL                           :: A(EQN%nVarTotal,EQN%nVarTotal) 
    REAL                           :: B(EQN%nVarTotal,EQN%nVarTotal) 
    REAL                           :: C(EQN%nVarTotal,EQN%nVarTotal)
    REAL                           :: EV(EQN%nVarTotal) 
    REAL                           :: LambdaA(EQN%nVarTotal) 
    REAL                           :: dA(EQN%nVarTotal,EQN%nVarTotal) 
    REAL                           :: dB(EQN%nVarTotal,EQN%nVarTotal) 
    REAL                           :: dC(EQN%nVarTotal,EQN%nVarTotal)
    REAL                           :: U0S(EQN%nVarTotal) 
    REAL                           :: U0C(EQN%nVarTotal) 
    REAL                           :: Ux(EQN%nVarTotal) 
    REAL                           :: Uy(EQN%nVarTotal) 
    REAL                           :: Uz(EQN%nVarTotal)
    REAL                           :: kx,ky,kz,omega
    REAL                           :: rho,mu,lambda,cp,cs 
    REAL                           :: n(3)
    REAL                           :: src_derivatives(LocDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)
    REAL                           :: src0S(LocDegFr,EQN%nVarTotal)
    REAL                           :: src0C(LocDegFr,EQN%nVarTotal)
    REAL                           :: dSS(EQN%nVarTotal),dSC(EQN%nVarTotal)
    REAL                           :: sig(0:10)
    REAL                           :: JacobiDetVolum 
    REAL, POINTER                  :: cPoly3D(:,:,:,:,:)         => NULL()
    INTEGER, POINTER               :: NonZeroCPoly(:,:)          => NULL()
    INTEGER, POINTER               :: NonZeroCPolyIndex(:,:,:,:) => NULL()
    INTEGER                        :: nVert
    INTEGER                        :: iVar, iDegFr, iDirac, iRicker
    INTEGER                        :: iTimeGP, iIntGP
    INTEGER                        :: k
    INTEGER                        :: LocElemType
    INTEGER                        :: iWindow
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: iElem, time, dt
    INTENT(OUT)   :: src,src_deriv 
    INTENT(INOUT) :: OptionalFields
    !-------------------------------------------------------------------------!
    src = 0.0
    !
    LocElemType  = MESH%LocalElemType(iElem)
    SELECT CASE(LocElemType)
    CASE(4)
        cPoly3D           => DISC%Galerkin%cPoly3D_Tet
        NonZeroCPoly      => DISC%Galerkin%NonZeroCPoly_Tet
        NonZeroCPolyIndex => DISC%Galerkin%NonZeroCPolyIndex_Tet
        nVert             =  MESH%nVertices_Tet
        JacobiDetVolum    =  6.0d0 * MESH%ELEM%Volume(iElem)
    CASE(6)
        cPoly3D           => DISC%Galerkin%cPoly3D_Hex
        NonZeroCPoly      => DISC%Galerkin%NonZeroCPoly_Hex
        NonZeroCPolyIndex => DISC%Galerkin%NonZeroCPolyIndex_Hex
        nVert             =  MESH%nVertices_Tet
        JacobiDetVolum    =  MESH%ELEM%Volume(iElem)
    ENDSELECT
    !
    x(1:nVert) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:nVert,iElem))
    y(1:nVert) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:nVert,iElem))
    z(1:nVert) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:nVert,iElem))
    !
    SELECT CASE(SOURCE%Type)
    CASE(1)
        ! Source term needed for the convergence studies
        src_derivatives = 0.
        src             = 0.
        IF(PRESENT(src_deriv)) THEN
            src_deriv = 0.
        ENDIF
        !
        sig(:) = (/ 1., -1., -1., 1., 1., -1., -1., 1., 1., -1., -1. /) 
        !
        A = 0.
        B = 0.
        C = 0.
        !
        rho    = EQN%rho0
        mu     = EQN%mu
        lambda = EQN%Lambda
        cp     = SQRT((lambda+2*mu)/rho)
        cs     = SQRT(mu/rho)
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
        n(:)   = IC%PW%n(:)                                              
        LambdaA(:) = (/ -cp, -cs, -cs, 0., 0., 0., cs, cs, cp /)
        !
        EV(1) = -EQN%rho0*(-2*n(2)**2*EQN%mu-2*n(3)**2*EQN%mu+EQN%lambda+2*EQN%mu)
        EV(2) = -EQN%rho0*(2*n(2)**2*EQN%mu+EQN%lambda)
        EV(3) = -EQN%rho0*(2*n(3)**2*EQN%mu+EQN%lambda)
        EV(4) = -2*n(2)*EQN%mu*n(1)*EQN%rho0
        EV(5) = -2*EQN%mu*n(2)*EQN%rho0*n(3)
        EV(6) = -2*EQN%mu*n(1)*EQN%rho0*n(3)
        EV(7) = n(1)*SQRT(EQN%rho0*(EQN%lambda+2*EQN%mu))
        EV(8) = n(2)*sqrt(EQN%rho0*(EQN%lambda+2*EQN%mu))
        EV(9) = SQRT(EQN%rho0*(EQN%lambda+2*EQN%mu))*n(3)
        !
        kx = SOURCE%CS%k1(1) 
        ky = SOURCE%CS%k1(2) 
        kz = SOURCE%CS%k1(3) 
        !
        omega = LambdaA(9)*SQRT(SUM(IC%PW%k_vec(:)**2))
        !
        src0C = 0.
        src0S = 0.
        DO iIntGP = 1, DISC%Galerkin%nIntGP
          xi   = DISC%Galerkin%intGaussP(1,iIntGP)
          eta  = DISC%Galerkin%intGaussP(2,iIntGP)
          zeta = DISC%Galerkin%intGaussP(3,iIntGP)
          CALL HexaTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,x,y,z)                  
          dA(:,:) = A(:,:)*SOURCE%CS%U0(1)*SIN(kx*xGP+ky*yGP+kz*zGP) 
          dB(:,:) = B(:,:)*SOURCE%CS%U0(2)*SIN(kx*xGP+ky*yGP+kz*zGP) 
          dC(:,:) = C(:,:)*SOURCE%CS%U0(3)*SIN(kx*xGP+ky*yGP+kz*zGP) 
          U0S = EV*IC%PW%ampfield(1)*SIN( IC%PW%k_vec(1)*xGP + IC%PW%k_vec(2)*yGP + IC%PW%k_vec(3)*zGP - omega*time ) 
          U0C = EV*IC%PW%ampfield(1)*COS( IC%PW%k_vec(1)*xGP + IC%PW%k_vec(2)*yGP + IC%PW%k_vec(3)*zGP - omega*time ) 
          Ux   = U0S * IC%PW%k_vec(1)
          Uy   = U0S * IC%PW%k_vec(2)
          Uz   = U0S * IC%PW%k_vec(3)
          dSS  = MATMUL(dA,Ux) + MATMUL(dB,Uy) + MATMUL(dC,Uz)
          Ux   = U0C * IC%PW%k_vec(1)
          Uy   = U0C * IC%PW%k_vec(2)
          Uz   = U0C * IC%PW%k_vec(3)
          dSC  = MATMUL(dA,Ux) + MATMUL(dB,Uy) + MATMUL(dC,Uz)
          DO iDegFr = 1, DISC%Galerkin%nDegFr
            phi = DISC%Galerkin%IntGPBaseFunc(iDegFr,iIntGP,DISC%Galerkin%nPoly)
            src0S(iDegFr,:) = src0S(iDegFr,:) + DISC%Galerkin%IntGaussW(iIntGP)*phi*dSS
            src0C(iDegFr,:) = src0C(iDegFr,:) + DISC%Galerkin%IntGaussW(iIntGP)*phi*dSC
          ENDDO
        ENDDO
        !
        src(:,:) = 0. 
        DO k = 0, DISC%Galerkin%nPoly
            src_derivatives(:,:,k) = sig(k)*(-omega)**k*( src0C(:,:)*MOD(k+1,2) + src0S(:,:)*MOD(k,2) ) 
            src(:,:)               = src(:,:) + dt**(k+1)/DISC%Galerkin%Faculty(k+1)*src_derivatives(:,1:EQN%nVar,k) 
        ENDDO
        IF(PRESENT(src_deriv)) THEN
            DO iDegFr = 1, LocDegFr
               src_deriv(iDegFr,:,:) = src_derivatives(iDegFr,:,:) * DISC%Galerkin%iMassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly)
            ENDDO
        ENDIF
        !
        src = src * JacobiDetVolum
        !
    CASE(16,18)
          ! By default no source. Only when Dirac acts.
          !
          src(:,:) = 0.
          !
          ! Ricker source
          !
          DO iRicker = 1, SOURCE%Ricker%nRicker
             IF( SOURCE%Ricker%Element(iRicker).EQ.iElem) THEN
                !
                CALL TrafoXYZ2XiEtaZeta(                                              &
                                    xi    = xi,                                       &
                                    eta   = eta,                                      &
                                    zeta  = zeta,                                     &
                                    xP    = SOURCE%Ricker%SpacePosition(1,iRicker),   &
                                    yP    = SOURCE%Ricker%SpacePosition(2,iRicker),   &
                                    zP    = SOURCE%Ricker%SpacePosition(3,iRicker),   &
                                    x     = x,                                        &
                                    y     = y,                                        &
                                    z     = z,                                        &
                                    vType = MESH%LocalVrtxType(iElem)                 )

                !
                iVar = SOURCE%Ricker%EqnNr(iRicker)
                !
                a1  = SOURCE%Ricker%a1(iRicker)
                f   = SOURCE%Ricker%f(iRicker)
                a2  = -(EQN%Pi*f)*(EQN%Pi*f)
                hw2 = 1.0/(f*f)
                !
                ! Do Gaussian integration of the source term in time
                !
                integral = 0.
                DO iTimeGP=1,DISC%Galerkin%nTimeGP
                   tau      = time + DISC%Galerkin%TimeGaussP(iTimeGP) -          &
                                  SOURCE%Ricker%Delay(iRicker)
                   !
                   ! "Value" contains the evaluation of the Ricker wavelet at the Gaussian point
                   !
                   IF(SOURCE%Type.EQ.16) THEN
                     value = a1*(0.5+a2*tau*tau)*EXP(a2*tau*tau)/OptionalFields%BackgroundValue(iElem,1)
                   ELSEIF(SOURCE%Type.EQ.18) THEN
                     value = a1*EXP(-tau*tau*hw2)
                   ENDIF
                   !
                   ! Integral contains the time integral of the source term
                   !
                   integral = integral + DISC%Galerkin%TimeGaussW(iTimeGP)*value
                ENDDO
                !
                ! Value is already integrated in time. This must be taken
                ! into account by the calling subroutine !
                !
                DO iDegFr = 1, LocDegFr
                   CALL BaseFunc3D(phi,iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly,cPoly3D,NonZeroCPoly,NonZeroCPolyIndex)
                   src(iDegFr,iVar) = src(iDegFr,iVar) + phi*integral
                ENDDO
             ENDIF
          ENDDO
          !
    CASE(20)
             ! By default no source. Only when Dirac acts.
             !
             src(:,:) = 0.
             !
             ! Point source with individual source time function
             ! specified in the FSRM-file read from readpar.f90
             !
             IF(time.LT.SOURCE%RP%T_max) THEN

               DO iRicker = 1, SOURCE%Ricker%nRicker
                  IF( SOURCE%RP%Element(iRicker).EQ.iElem) THEN
                     !
                     CALL TrafoXYZ2XiEtaZeta(                                        &
                                       xi    = xi,                                   &
                                       eta   = eta,                                  &
                                       zeta  = zeta,                                 &
                                       xP    = SOURCE%RP%SpacePosition(1,iRicker),   &
                                       yP    = SOURCE%RP%SpacePosition(2,iRicker),   &
                                       zP    = SOURCE%RP%SpacePosition(3,iRicker),   &
                                       x     = x,                                    &
                                       y     = y,                                    &
                                       z     = z,                                    &
                                       vType = MESH%LocalVrtxType(iElem)             )

                     iVar = SOURCE%Ricker%EqnNr(iRicker)
                     !
                     ! Gaussian integration of the source term in time
                     integral = 0.
                     DO iTimeGP=1,DISC%Galerkin%nTimeGP
                        !
                        t      = time + DISC%Galerkin%TimeGaussP(iTimeGP)
                        !
                        ! "Value" contains the evaluation of the Source time function at the Gaussian point
                        !         and is interpolated linearly between discrete values given through the
                        !         time history
                        !
                        t_proj = t
                        !
                        value  = 0.
                        !
                        IF(t_proj.LT.SOURCE%RP%T_max) THEN

                           ! We assume equispaced sampling of the Time History!
                           i_left  = FLOOR(t_proj / SOURCE%RP%t_samp)+1
                           i_right = i_left + 1

                           t_left  = (i_left-1) * SOURCE%RP%t_samp
                           d_t     = t_proj - t_left
                           !
                           ! Linear interpolation between samples of the time history
                           m     = ( SOURCE%RP%TimeHist(i_right,iRicker) - SOURCE%RP%TimeHist(i_left,iRicker) ) / SOURCE%RP%t_samp
                           value = SOURCE%RP%TimeHist(i_left,iRicker) + m * d_t
                           !
                        ENDIF
                        !
                        value = value/OptionalFields%BackgroundValue(iElem,1)
                        integral = integral + DISC%Galerkin%TimeGaussW(iTimeGP)*value
                        !
                     ENDDO
                     !
                     ! Value is already integrated in time. This must be taken
                     ! into account by the calling subroutine !
                     !
                     DO iDegFr = 1, LocDegFr
                       CALL BaseFunc3D(phi,iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly,cPoly3D,NonZeroCPoly,NonZeroCPolyIndex)
                       src(iDegFr,iVar) = src(iDegFr,iVar) + phi*integral
                     ENDDO
                  ENDIF
               ENDDO
             ENDIF
             !
    CASE(30)
             !
             ! By default no source. Only when Dirac acts.
             !
             src(:,:) = 0.
             !
             SELECT CASE(SOURCE%RP%Type)
                !
             CASE(3)       !Only one Rupture Plane Segment with StrikeSlip, DipSlip, Onset times and Rise times (triangular sliprate!)
                !
                IF(MESH%IncludesFSRP(iElem)) THEN
                    DO iDirac = 1, SOURCE%RP%nxRP*SOURCE%RP%nzRP
                       IF( SOURCE%RP%Element(iDirac).EQ.iElem) THEN
                           !
                           CALL TrafoXYZ2XiEtaZeta(                                              &
                                               xi    = xi,                                       &
                                               eta   = eta,                                      &
                                               zeta  = zeta,                                     &
                                               xP    = SOURCE%RP%SpacePosition(1,iDirac),        &
                                               yP    = SOURCE%RP%SpacePosition(2,iDirac),        &
                                               zP    = SOURCE%RP%SpacePosition(3,iDirac),        &
                                               x     = x,                                        &
                                               y     = y,                                        &
                                               z     = z,                                        &
                                               vType = MESH%LocalVrtxType(iElem)                 )

                           !
                           ! Do Gaussian integration of the source term in time
                           !
                           integral = 0.

                           DO iTimeGP=1,DISC%Galerkin%nTimeGP
                              !
                              tau      = time + DISC%Galerkin%TimeGaussP(iTimeGP)
                              !                           !
                              ! "Value" contains the evaluation of the Source time function at the Gaussian point
                              !
                              value = 0.
                              DO iWindow = 1,SOURCE%RP%nTWindow
                                  ! left (up-going) ramp of triangle
                                  a1    = heaviside( tau -  SOURCE%RP%Tonset(iDirac) )
                                  a2    = heaviside( tau - (SOURCE%RP%Tonset(iDirac)+SOURCE%RP%TRise(iDirac)*0.5) )
                                  fac   = a1 - a2
                                  m     = 2.*SOURCE%RP%Sliprate(iDirac,iWindow) / (SOURCE%RP%TRise(iDirac)*0.5)
                                  value = value + fac * m * (tau - SOURCE%RP%Tonset(iDirac))
                                  ! right (down-going) ramp of triangle
                                  a1    = heaviside( tau - (SOURCE%RP%Tonset(iDirac)+SOURCE%RP%TRise(iDirac)*0.5) )
                                  a2    = heaviside( tau - (SOURCE%RP%Tonset(iDirac)+SOURCE%RP%TRise(iDirac)    ) )
                                  fac   = a1 - a2
                                  m     = -m
                                  value = value + fac * ( 2.*SOURCE%RP%Sliprate(iDirac,iWindow) +   &
                                                    m * (tau - (SOURCE%RP%Tonset(iDirac)+SOURCE%RP%TRise(iDirac)*0.5)) )
                              ENDDO
                              !
                              ! Then multiply by rigidity \mu and subfault area A. (M_0 = D_t * \mu * A)
                              value = value * OptionalFields%BackgroundValue(iElem,2) * SOURCE%RP%dxRP*SOURCE%RP%dzRP
                              !
                              ! Integral contains the time integral of the source term
                              !
                              integral = integral + DISC%Galerkin%TimeGaussW(iTimeGP)*value
                              !
                           ENDDO
                           !
                           ! Get orientation angles
                           !
                           strk = SOURCE%RP%strk
                           dip  = SOURCE%RP%dip
                           rake = SOURCE%RP%rake(iDirac,1)
                           !
                           ! Set moment tensor in epicentral coordinate system (with double couple assumption)
                           ! Ref: Jose Pujol; Elastic Wave Propagation and Generation in Seismology, Cambridge, 2003.
                           !      Chapter 10.7
                           !
                           SOURCE%RP%MomentTensor(1,1) = sin(dip)*cos(rake)*sin(2.*strk)+sin(2.*dip)*sin(rake)*sin(strk)**2.
                           SOURCE%RP%MomentTensor(1,2) = sin(dip)*cos(rake)*cos(2.*strk)+sin(2.*dip)*sin(rake)*sin(2.*strk)/2.
                           SOURCE%RP%MomentTensor(1,3) = cos(dip)*cos(rake)*cos(strk)   +cos(2.*dip)*sin(rake)*sin(strk)
                           SOURCE%RP%MomentTensor(2,1) = SOURCE%RP%MomentTensor(1,2)
                           SOURCE%RP%MomentTensor(2,2) = -sin(dip)*cos(rake)*sin(2.*strk)+sin(2.*dip)*sin(rake)*cos(strk)**2.
                           SOURCE%RP%MomentTensor(2,3) = -cos(dip)*cos(rake)*sin(strk)   +cos(2.*dip)*sin(rake)*cos(strk)
                           SOURCE%RP%MomentTensor(3,1) = SOURCE%RP%MomentTensor(1,3)
                           SOURCE%RP%MomentTensor(3,2) = SOURCE%RP%MomentTensor(2,3)
                           SOURCE%RP%MomentTensor(3,3) = -sin(2.*dip)*sin(rake)
                           !
                           ! Value is already integrated in time. This must be taken
                           ! into account by the calling subroutine !
                           !
                           ! Distribution of the Seismic Moment Tensor as source term for the stresses
                           !
                           DO iDegFr = 1, LocDegFr
                              CALL BaseFunc3D(phi,iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly,cPoly3D,NonZeroCPoly,NonZeroCPolyIndex)
                              src(iDegFr,1) = src(iDegFr,1) + phi*integral*SOURCE%RP%MomentTensor(1,1)
                              src(iDegFr,2) = src(iDegFr,2) + phi*integral*SOURCE%RP%MomentTensor(2,2)
                              src(iDegFr,3) = src(iDegFr,3) + phi*integral*SOURCE%RP%MomentTensor(3,3)
                              src(iDegFr,4) = src(iDegFr,4) + phi*integral*SOURCE%RP%MomentTensor(1,2)
                              src(iDegFr,5) = src(iDegFr,5) + phi*integral*SOURCE%RP%MomentTensor(2,3)
                              src(iDegFr,6) = src(iDegFr,6) + phi*integral*SOURCE%RP%MomentTensor(1,3)
                           ENDDO

                       ENDIF
                    ENDDO
                ENDIF
                !
             CASE(4)       !Only one Rupture Plane Segment with StrikeSlip, DipSlip, Onset times and Rise times (boxcar sliprate!)
                !
                IF(MESH%IncludesFSRP(iElem)) THEN
                    DO iDirac = 1, SOURCE%RP%nxRP*SOURCE%RP%nzRP
                       IF( SOURCE%RP%Element(iDirac).EQ.iElem) THEN
                           !
                           CALL TrafoXYZ2XiEtaZeta(                                              &
                                               xi    = xi,                                       &
                                               eta   = eta,                                      &
                                               zeta  = zeta,                                     &
                                               xP    = SOURCE%RP%SpacePosition(1,iDirac),        &
                                               yP    = SOURCE%RP%SpacePosition(2,iDirac),        &
                                               zP    = SOURCE%RP%SpacePosition(3,iDirac),        &
                                               x     = x,                                        &
                                               y     = y,                                        &
                                               z     = z,                                        &
                                               vType = MESH%LocalVrtxType(iElem)                 )
                             !
                           ! Do Gaussian integration of the source term in time
                           !
                           integral = 0.

                           DO iTimeGP=1,DISC%Galerkin%nTimeGP
                              !
                              tau      = time + DISC%Galerkin%TimeGaussP(iTimeGP)
                              !                           !
                              ! "Value" contains the evaluation of the Source time function at the Gaussian point
                              !
                              value = 0.
                              a1    = heaviside( tau -  SOURCE%RP%Tonset(iDirac) )
                              a2    = heaviside( tau - (SOURCE%RP%Tonset(iDirac)+SOURCE%RP%TRise(iDirac)) )
                              fac   = a1 - a2
                              value = value + fac * SOURCE%RP%Sliprate(iDirac,1)
                              !
                              ! Then multiply by rigidity \mu and subfault area A. (M_0 = D_t * \mu * A)
                              value = value * OptionalFields%BackgroundValue(iElem,2) * SOURCE%RP%dxRP*SOURCE%RP%dzRP
                              !
                              ! Integral contains the time integral of the source term
                              !
                              integral = integral + DISC%Galerkin%TimeGaussW(iTimeGP)*value
                              !
                           ENDDO
                           !
                           ! Get orientation angles
                           !
                           strk = SOURCE%RP%strk
                           dip  = SOURCE%RP%dip
                           rake = SOURCE%RP%rake(iDirac,1)
                           !
                           ! Set moment tensor in epicentral coordinate system (with double couple assumption)
                           ! Ref: Jose Pujol; Elastic Wave Propagation and Generation in Seismology, Cambridge, 2003.
                           !      Chapter 10.7
                           !
                           SOURCE%RP%MomentTensor(1,1) = sin(dip)*cos(rake)*sin(2.*strk)+sin(2.*dip)*sin(rake)*sin(strk)**2.
                           SOURCE%RP%MomentTensor(1,2) = sin(dip)*cos(rake)*cos(2.*strk)+sin(2.*dip)*sin(rake)*sin(2.*strk)/2.
                           SOURCE%RP%MomentTensor(1,3) = cos(dip)*cos(rake)*cos(strk)   +cos(2.*dip)*sin(rake)*sin(strk)
                           SOURCE%RP%MomentTensor(2,1) = SOURCE%RP%MomentTensor(1,2)
                           SOURCE%RP%MomentTensor(2,2) = -sin(dip)*cos(rake)*sin(2.*strk)+sin(2.*dip)*sin(rake)*cos(strk)**2.
                           SOURCE%RP%MomentTensor(2,3) = -cos(dip)*cos(rake)*sin(strk)   +cos(2.*dip)*sin(rake)*cos(strk)
                           SOURCE%RP%MomentTensor(3,1) = SOURCE%RP%MomentTensor(1,3)
                           SOURCE%RP%MomentTensor(3,2) = SOURCE%RP%MomentTensor(2,3)
                           SOURCE%RP%MomentTensor(3,3) = -sin(2.*dip)*sin(rake)
                           !
                           ! Value is already integrated in time. This must be taken
                           ! into account by the calling subroutine !
                           !
                           ! Distribution of the Seismic Moment Tensor as source term for the stresses
                           !
                           DO iDegFr = 1, LocDegFr
                              CALL BaseFunc3D(phi,iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly,cPoly3D,NonZeroCPoly,NonZeroCPolyIndex)
                              src(iDegFr,1) = src(iDegFr,1) + phi*integral*SOURCE%RP%MomentTensor(1,1)
                              src(iDegFr,2) = src(iDegFr,2) + phi*integral*SOURCE%RP%MomentTensor(2,2)
                              src(iDegFr,3) = src(iDegFr,3) + phi*integral*SOURCE%RP%MomentTensor(3,3)
                              src(iDegFr,4) = src(iDegFr,4) + phi*integral*SOURCE%RP%MomentTensor(1,2)
                              src(iDegFr,5) = src(iDegFr,5) + phi*integral*SOURCE%RP%MomentTensor(2,3)
                              src(iDegFr,6) = src(iDegFr,6) + phi*integral*SOURCE%RP%MomentTensor(1,3)
                           ENDDO

                       ENDIF
                    ENDDO
                ENDIF
                !
             CASE DEFAULT
                !
                logError(*)  'The format type of the Finite Source Rupture Model is unknown! '
                STOP                                     ! STOP
                !
             END SELECT
             !
    CASE(50) ! Finite rupture planes with individual source time function 
             ! specified in the FSRM-file read from readpar.f90
             !
             src(:,:) = 0.
             !
             IF(MESH%IncludesFSRP(iElem).AND.time.LT.SOURCE%RP%T_max) THEN
                 DO iDirac = 1, SOURCE%RP%nSbfs(1)
                    IF( SOURCE%RP%Element(iDirac).EQ.iElem) THEN
                        !
                        CALL TrafoXYZ2XiEtaZeta(                                              &
                                            xi    = xi,                                       &
                                            eta   = eta,                                      &
                                            zeta  = zeta,                                     &
                                            xP    = SOURCE%RP%SpacePosition(1,iDirac),        &
                                            yP    = SOURCE%RP%SpacePosition(2,iDirac),        &
                                            zP    = SOURCE%RP%SpacePosition(3,iDirac),        &
                                            x     = x,                                        &
                                            y     = y,                                        &
                                            z     = z,                                        &
                                            vType = MESH%LocalVrtxType(iElem)                 )
                        !
                        ! Gaussian integration of the source term in time
                        M0_integral = 0.

                        DO iTimeGP=1,DISC%Galerkin%nTimeGP
                           !
                           t      = time + DISC%Galerkin%TimeGaussP(iTimeGP)
                           !
                           ! "Value" contains the evaluation of the Source time function at the Gaussian point
                           !         and is interpolated linearly between discrete values given through the 
                           !         time history
                           !
                           t_proj = t
                           !
                           value  = 0.
                           !
                           IF(t_proj.LT.SOURCE%RP%T_max) THEN

                              ! We assume equispaced sampling of the Time History!
                              i_left  = FLOOR(t_proj / SOURCE%RP%t_samp)+1
                              i_right = i_left + 1

                              t_left  = (i_left-1) * SOURCE%RP%t_samp
                              d_t     = t_proj - t_left
                              !
                              ! Linear interpolation between samples of the time history
                              m     = ( SOURCE%RP%TimeHist(i_right,iDirac) - SOURCE%RP%TimeHist(i_left,iDirac) ) / SOURCE%RP%t_samp           
                              value = SOURCE%RP%TimeHist(i_left,iDirac) + m * d_t
                              !
                           ENDIF  
                           !
                           ! Then multiply by a factor including the subfault area and already the rigidity value
                           value = value * SOURCE%RP%Area(iDirac)
                           !
                           M0_integral = M0_integral + DISC%Galerkin%TimeGaussW(iTimeGP)*value
                           !
                        ENDDO
                        !
                        ! Get orientation angles of the subfault
                        !
                        strk = SOURCE%RP%strks(1,iDirac)
                        dip  = SOURCE%RP%dips(1,iDirac)
                        rake = SOURCE%RP%rake(1,iDirac)
                        !
                        ! Set moment tensor in epicentral coordinate system (with double couple assumption)
                        !
                        SOURCE%RP%TensorRotation(1,1) =  cos(rake)*cos(strk) + sin(rake)*cos(dip)*sin(strk)
                        SOURCE%RP%TensorRotation(1,2) = -cos(rake)*sin(strk) + sin(rake)*cos(dip)*cos(strk)
                        SOURCE%RP%TensorRotation(1,3) = -sin(rake)*sin(dip)
                        SOURCE%RP%TensorRotation(2,1) = -sin(rake)*cos(strk) + cos(rake)*cos(dip)*sin(strk)
                        SOURCE%RP%TensorRotation(2,2) =  sin(rake)*sin(strk) + cos(rake)*cos(dip)*cos(strk)
                        SOURCE%RP%TensorRotation(2,3) = -cos(rake)*sin(dip)
                        SOURCE%RP%TensorRotation(3,1) =  sin(dip) *sin(strk)
                        SOURCE%RP%TensorRotation(3,2) =  sin(dip) *cos(strk)
                        SOURCE%RP%TensorRotation(3,3) =  cos(dip)

                        SOURCE%RP%TensorRotationT     = TRANSPOSE(SOURCE%RP%TensorRotation)

                        LocalMomentTensor             = MATMUL(SOURCE%RP%TensorRotationT, &
                                                        MATMUL(SOURCE%RP%MomentTensor,SOURCE%RP%TensorRotation))
                        !
                        ! Value is already integrated in time. This must be taken
                        ! into account by the calling subroutine
                        !
                        ! Distribution of the Seismic Moment Tensor as source term for the stresses
                        !
                        DO iDegFr = 1, LocDegFr
                           CALL BaseFunc3D(phi,iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly,cPoly3D,NonZeroCPoly,NonZeroCPolyIndex)
                           src(iDegFr,1) = src(iDegFr,1) + phi*M0_integral*LocalMomentTensor(1,1)
                           src(iDegFr,2) = src(iDegFr,2) + phi*M0_integral*LocalMomentTensor(2,2)
                           src(iDegFr,3) = src(iDegFr,3) + phi*M0_integral*LocalMomentTensor(3,3)
                           src(iDegFr,4) = src(iDegFr,4) + phi*M0_integral*LocalMomentTensor(1,2)
                           src(iDegFr,5) = src(iDegFr,5) + phi*M0_integral*LocalMomentTensor(2,3)
                           src(iDegFr,6) = src(iDegFr,6) + phi*M0_integral*LocalMomentTensor(1,3)
                        ENDDO

                    ENDIF
                 ENDDO
             ENDIF
             !
    END SELECT
    !
!    NULLIFY(cPoly3D)
!    NULLIFY(NonZeroCPoly)
!    NULLIFY(NonZeroCPolyIndex)
    !
  END SUBROUTINE GalerkinSource3D

END MODULE Galerkin_source_mod

