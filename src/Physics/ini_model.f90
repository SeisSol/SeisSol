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
  interface fitAttenuation
    module procedure fitAttenuation
  end interface
  INTERFACE ini_MODEL
     MODULE PROCEDURE ini_MODEL
  END INTERFACE
  INTERFACE ini_ATTENUATION
     MODULE PROCEDURE ini_ATTENUATION
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: fitAttenuation,           &
             ini_MODEL,                &
             ini_ATTENUATION
CONTAINS
  subroutine fitAttenuation(material, materialFitted, EQN)
    implicit none
    !--------------------------------------------------------------------------
    TYPE (tEquations)   :: EQN
    real                :: material(EQN%nAneMaterialVar)
    real                :: materialFitted(EQN%nBackgroundVar)
    !--------------------------------------------------------------------------
    real                :: Material_INF(2)
    real                :: Theta(EQN%nMechanisms,3)
    real                :: w_freq(EQN%nMechanisms)
    integer             :: iMech
    !--------------------------------------------------------------------------
    intent(in)          :: material, EQN
    intent(out)         :: materialFitted

    call ini_ATTENUATION(Theta, w_freq, Material_INF, material, EQN)
    materialFitted(1) = material(1)
    materialFitted(2:3) = Material_INF(:)
    do iMech=1, EQN%nMechanisms
      materialFitted(4+4*(iMech-1)) = w_freq(iMech)
      materialFitted(4+4*(iMech-1)+1:4+4*(iMech-1)+3) = Theta(iMech,:)
    end do
  end subroutine

  SUBROUTINE ini_MODEL(MaterialVal,EQN,MESH,IO,DISC,BND)
    !--------------------------------------------------------------------------

    USE COMMON_operators_mod, ONLY: OpenFile, XYinTriangle
    USE TrilinearInterpolation_mod
    USE ini_model_DR_mod
    use modules
    use f_ftoc_bind_interoperability

    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tEquations)               :: EQN
    TYPE (tUnstructMesh)            :: MESH
    TYPE (tInputOutput)             :: IO
    TYPE (tDiscretization)          :: DISC
    TYPE (tBoundary)                :: BND
    REAL                            :: MaterialVal(MESH%nElem,EQN%nBackgroundVar)
    !--------------------------------------------------------------------------
    integer                         :: iElem
    real                            :: material(EQN%nAneMaterialVar)
    real                            :: materialFitted(EQN%nBackgroundVar)
    !--------------------------------------------------------------------------
    INTENT(IN)                      :: MESH
    INTENT(OUT)                     :: MaterialVal
    INTENT(INOUT)                   :: EQN
    ! -------------------------------------------------------------------------

    IF(EQN%nBackgroundVar.LE.0) THEN
       logInfo(*) 'No values specified. Exiting ini_MODEL.f90. '
       RETURN
    ENDIF

    ! Call the pre model hooks
    call call_hook_pre_model()
    
    IF (EQN%Plasticity .NE. 0) THEN
      allocate (EQN%BulkFriction(MESH%nElem), EQN%PlastCo(MESH%nElem), EQN%IniStress(6,MESH%nElem))
    else
      allocate (EQN%BulkFriction(0), EQN%PlastCo(0), EQN%IniStress(0,0))
    ENDIF

    call c_interoperability_initializeModel(trim(EQN%MaterialFileName) // c_null_char, EQN%Anelasticity, EQN%Plasticity, MaterialVal, EQN%BulkFriction, EQN%PlastCo, EQN%IniStress)
    
    if (EQN%Anelasticity == 1) then
      do iElem=1, MESH%nElem
        material(:) = MaterialVal(iElem,1:EQN%nAneMaterialVar)
        call fitAttenuation(material, materialFitted, EQN)
        MaterialVal(iElem,:) = materialFitted(:)
      end do
    end if

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
