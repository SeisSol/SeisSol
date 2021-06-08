!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Josep de la Puente Alvarez (josep.delapuente AT bsc.es, http://www.geophysik.uni-muenchen.de/Members/jdelapuente)
!!
!! @section LICENSE
!! Copyright (c) 2008-2013, SeisSol Group
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

module DGBasis_mod
  !---------------------------------------------------------------------------!
  use TypesDef
  !---------------------------------------------------------------------------!
  implicit none
  private
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  interface XYInElement
     module procedure XYInElement
  end interface
  interface TrafoXY2XiEta
     module procedure TrafoXY2XiEta
  end interface
  interface TrafoXiEta2XY
     module procedure TrafoXiEta2XY
  end interface
  interface TriTrafoXY2XiEta
     module procedure TriTrafoXY2XiEta
  end interface
  interface QuadTrafoXY2XiEta
     module procedure QuadTrafoXY2XiEta
  end interface
  interface TriTrafoXiEta2XY
     module procedure TriTrafoXiEta2XY
  end interface
  interface QuadTrafoXiEta2XY
     module procedure QuadTrafoXiEta2XY
  end interface
  interface TrafoXYZ2XiEtaZeta
     module procedure TrafoXYZ2XiEtaZeta
  end interface
  interface TetraTrafoXiEtaZeta2XYZ
     module procedure TetraTrafoXiEtaZeta2XYZ
  end interface
  interface HexaTrafoXiEtaZeta2XYZ
     module procedure HexaTrafoXiEtaZeta2XYZ
  end interface
  interface TriTrafoChi2XiEta
     module procedure TriTrafoChi2XiEta
  end interface
  interface QuadTrafoChi2XiEta
     module procedure QuadTrafoChi2XiEta
  end interface
  interface TrafoChiTau2XiEtaZeta
     module procedure TrafoChiTau2XiEtaZeta
  end interface
  interface HexaTrafoChiTau2XiEtaZeta
     module procedure HexaTrafoChiTau2XiEtaZeta
  end interface
  interface HexaTrafoXiEtaZetaGrad
     module procedure HexaTrafoXiEtaZetaGrad
  end interface
  interface TetraTrafoXiEtaZetaGrad
     module procedure TetraTrafoXiEtaZetaGrad
  end interface
  interface QuadTrafoXiEtaGrad
     module procedure QuadTrafoXiEtaGrad
  end interface
  interface TriTrafoXiEtaGrad
     module procedure TriTrafoXiEtaGrad
  end interface
  interface SpaceTimeBasis
     module procedure SpaceTimeBasis
  end interface
  interface HexaTrafoXiEtaZetaGradGLL
     module procedure HexaTrafoXiEtaZetaGradGLL
  end interface
  interface Angle_Vector_Vector_2D
     module procedure Angle_Vector_Vector_2D
  end interface
  interface Mirror_P_atLine
     module procedure Mirror_P_atLine
  end interface

  !---------------------------------------------------------------------------!
  ! 2D Tri
  public  :: XYInElement
  public  :: TriTrafoChi2XiEta
  public  :: TrafoXiEta2XY
  public  :: TrafoXY2XiEta
  public  :: TriTrafoXiEta2XY
  public  :: TriTrafoXY2XiEta
  public  :: TriTrafoXiEtaGrad
  ! 2D Quad
  public  :: QuadTrafoXiEta2XY
  public  :: QuadTrafoXY2XiEta
  public  :: QuadTrafoChi2XiEta
  public  :: QuadTrafoXiEtaGrad
  ! 3D Tetra
  public  :: XYZInElement
  public  :: TetraTrafoXiEtaZeta2XYZ
  public  :: TetraTrafoXiEtaZetaGrad
  public  :: TrafoXYZ2XiEtaZeta
  public  :: TrafoChiTau2XiEtaZeta
  ! 3D Hexa
  public  :: HexaTrafoXiEtaZeta2XYZ
!  public  :: HexaTrafoXYZ2XiEtaZeta
  public  :: HexaTrafoChiTau2XiEtaZeta
  public  :: HexaTrafoXiEtaZetaGrad
  ! 3D Hexa GLL
  public  :: HexaTrafoXiEtaZetaGradGLL
  ! 2D/3D
  ! Local space-time DG
  public  :: SpaceTimeBasis
  ! Geometry
  public  :: Angle_Vector_Vector_2D
  public  :: Mirror_P_atLine
  !---------------------------------------------------------------------------!
  integer, parameter, private :: rk=kind(1.)
contains
 
  ! ***************************************************************** !
  ! *                                                               * !
  ! * XYInElement returns .TRUE. if the point (x/y) is inside the   * !
  ! * element iElem, else it returns .FALSE.                        * !
  ! *                                                               * !
  ! ***************************************************************** !

  function XYInElement(xS,yS,iElem,epsilon,MESH)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Function result type
    logical                  :: XYInElement
    ! Argument list declaration
    type(tUnstructMesh)      :: MESH
    type(tDiscretization)    :: DISC
    type(tInputOutput)       :: IO
    integer     :: iElem
    real        :: xS, yS
    real        :: epsilon
    real        :: refFactor
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    real        :: x(MESH%LocalElemType(iElem)), y(MESH%LocalElemType(iElem))
    real        :: xi, eta
    real        :: prod
    integer     :: i,k
    real        :: dist1, dist2
    !-------------------------------------------------------------------------!
    intent(in)  :: xS, yS, iElem, epsilon
    !-------------------------------------------------------------------------!
  
    do i=1,MESH%LocalElemType(iElem)
       x(i) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(i,iElem))
       y(i) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(i,iElem))
    enddo

    select case(MESH%LocalElemType(iElem))
    case(3) !Triangular element

      refFactor = 0.5

      xi  = refFactor/MESH%ELEM%Volume(iElem)*( ( x(MESH%LocalElemType(iElem))*y(1) - x(1)*y(MESH%LocalElemType(iElem)) ) + &
            xS*(y(MESH%LocalElemType(iElem))-y(1))                  + &
            yS*(x(1)-x(MESH%LocalElemType(iElem)))              )
      eta = refFactor/MESH%ELEM%Volume(iElem)*( ( x(1)*y(2) - x(2)*y(1) )         + &
            xS*(y(1)-y(2))                      + &
            yS*(x(2)-x(1))              )
      
      !-------------------------------------------------------------------------!
      ! Because of numerical errors, it is possible that a point which lies     !
      ! exactly on the boundary is considered outside the element               !
      ! So we set a tolerance value of epsilon which allows points lo lie       !
      ! slightly out of the element. This is necessary because fluxes are       !
      ! calculated on boundary points!                                          !
      ! For epsilon we choose 1e-5                                              !
      !-------------------------------------------------------------------------!
      
      if((xi.lt.(0.-epsilon)).or.(eta.lt.(0.-epsilon)).or. &
         (eta.gt.(1.-xi+epsilon))) then
          XYInElement = .false.
      else
          XYInElement = .true.
      endif
      !
    case(4) !Quadrilateral element
      !
      ! dist 1 gives a messure how far is the point (xS,yS) from the element
      dist1 = sqrt((xS-x(1))**2+(yS-y(1))**2)
      ! dist 2 gives a messure how big is the element
      dist2 = sqrt((x(3)-x(1))**2+(y(3)-y(1))**2)

      if(dist1<3*dist2) then

          call QuadTrafoXY2XiEta(xi,eta,xS,yS,x,y)

          if ((xi .lt. (0.-epsilon)).or.(eta .lt.(0.-epsilon)) .or.        &
              (xi .gt. (1.+epsilon)).or.(eta .gt.(1.+epsilon)))     then
              XYInElement = .false.
          else
              XYInElement = .true.
          endif

      else
          XYInElement = .false.
      endif

!      XYInElement = .TRUE.
!      DO i = 1,MESH%GlobalElemType
!         k    = MOD(i,4) + 1
!         prod = (x(k)-x(i))*(yS-y(i)) - (y(k)-y(i))*(xS-x(i))
!         
!         IF (prod.LT.epsilon) THEN
!           XYInElement = .FALSE.
!           EXIT
!         ENDIF
!      ENDDO
      !

      !
      call TrafoXY2XiEta(xi,eta,xS,yS,x,y,4)

      if ((xi .lt. (0.-epsilon)).or.(eta .lt.(0.-epsilon)) .or.        &
          (xi .gt. (1.+epsilon)).or.(eta .gt.(1.+epsilon)))     then
          XYInElement = .false.
      else
          XYInElement = .true.
      endif
      !
    end select
    ! 
  end function XYInElement

  ! ***************************************************************** !
  ! *                                                               * !
  ! * XYZInElement returns .TRUE. if the point (x/y/z) is inside    * !
  ! * the element iElem, else it returns .FALSE.                    * !
  ! *                                                               * !
  ! ***************************************************************** !

  function XYZInElement(xP,yP,zP,iElem,epsilon,MESH,DISC)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Function result type
    logical                  :: XYZInElement
    ! Argument list declaration
    type(tUnstructMesh)     :: MESH
    type(tDiscretization)   :: DISC
    integer     :: iElem
    real        :: xP, yP, zP
    real        :: epsilon
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    integer     :: iSide
    integer     :: VertexSide(MESH%LocalElemType(iElem),MESH%LocalElemType(iElem))           
    real        :: x(MESH%LocalElemType(iElem))
    real        :: y(MESH%LocalElemType(iElem))
    real        :: z(MESH%LocalElemType(iElem))
    real        :: xi, eta, zeta
    real        :: J
    real        :: prod(MESH%LocalElemType(iElem))
    real        :: vecP(3), x0(3)
    logical     :: inside
    !-------------------------------------------------------------------------!
    intent(IN)  :: xP, yP, zP, iElem, epsilon, MESH
    !-------------------------------------------------------------------------!
  
    select case(MESH%LocalElemType(iElem))
    case(4)!Tet
        x(:) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(:,iElem))
        y(:) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(:,iElem))
        z(:) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(:,iElem))
        !
        J = -z(1)*x(2)*y(3)+x(1)*z(3)*y(4)+x(1)*z(2)*y(3)+z(4)*x(2)*y(3)+             &
             z(3)*x(4)*y(2)+z(2)*y(4)*x(3)+z(4)*y(1)*x(3)-z(3)*y(1)*x(4)+z(3)*        &
             y(1)*x(2)-x(1)*z(2)*y(4)+x(1)*z(4)*y(2)+z(1)*x(2)*y(4)+z(2)*y(1)*        &
             x(4)-z(4)*y(1)*x(2)-x(1)*z(4)*y(3)-z(3)*y(4)*x(2)-z(1)*x(3)*y(4)-        &
             z(4)*x(3)*y(2)-z(2)*x(4)*y(3)+z(1)*x(3)*y(2)-z(2)*y(1)*x(3)-z(1)*x(4)*   &
             y(2)-x(1)*z(3)*y(2)+z(1)*x(4)*y(3)
        !
        xi = ( z(4)*y(3)+z(1)*y(4)-z(3)*y(4)-z(1)*y(3)-z(4)*y(1)+z(3)*y(1))/J*xP +    &
             (-x(1)*z(3)+z(4)*x(1)-z(4)*x(3)+z(1)*x(3)-z(1)*x(4)+x(4)*z(3))/J*yP +    &
             (-y(1)*x(3)+y(3)*x(1)+y(1)*x(4)+y(4)*x(3)-y(3)*x(4)-y(4)*x(1))/J*zP +    &
             ( x(1)*z(3)*y(4)-z(3)*y(1)*x(4)+z(1)*x(4)*y(3)-x(1)*z(4)*y(3) -          &
               z(1)*x(3)*y(4)+z(4)*y(1)*x(3))/J
        !
        eta = -( z(4)*y(2)-z(1)*y(2)+z(1)*y(4)+z(2)*y(1)-z(2)*y(4)-z(4)*y(1))/J*xP -  &
               ( z(2)*x(4)-z(2)*x(1)-z(1)*x(4)-z(4)*x(2)+z(4)*x(1)+z(1)*x(2))/J*yP -  &
               (-y(2)*x(4)+y(1)*x(4)+y(4)*x(2)-y(1)*x(2)-y(4)*x(1)+y(2)*x(1))/J*zP -  &
               (-z(2)*y(1)*x(4)+x(1)*z(2)*y(4)-x(1)*z(4)*y(2)-z(1)*x(2)*y(4) +        &
                 z(1)*x(4)*y(2)+z(4)*y(1)*x(2))/J
        !
        zeta = (z(1)*y(3)-z(3)*y(1)-z(1)*y(2)+z(3)*y(2)+z(2)*y(1)-z(2)*y(3))/J*xP +   &
               (x(1)*z(3)-z(2)*x(1)-z(3)*x(2)+z(1)*x(2)+z(2)*x(3)-z(1)*x(3))/J*yP +   &
               (y(2)*x(1)-y(3)*x(1)+x(2)*y(3)-y(1)*x(2)+y(1)*x(3)-x(3)*y(2))/J*zP +   &
               (-x(1)*z(3)*y(4)-z(4)*x(2)*y(3)+z(1)*x(3)*y(4)-z(3)*x(4)*y(2) -        &
                 z(2)*y(4)*x(3)-z(4)*y(1)*x(3)+z(3)*y(1)*x(4)-z(2)*y(1)*x(4) +        &
                 x(1)*z(2)*y(4)-x(1)*z(4)*y(2)-z(1)*x(2)*y(4)+z(1)*x(4)*y(2) +        &
                 z(4)*y(1)*x(2)+x(1)*z(4)*y(3)+z(3)*y(4)*x(2)+z(4)*x(3)*y(2) -        &
                 z(1)*x(4)*y(3)+z(2)*x(4)*y(3)+J)/J

        !-------------------------------------------------------------------------!
        ! Because of numerical errors, it is possible that a point which lies     !
        ! exactly on the boundary is considered outside the element               !
        ! So we set a tolerance value of epsilon which allows points lo lie       !
        ! slightly out of the element.                                            ! 
        !-------------------------------------------------------------------------!

        if( (xi  .lt.(0.-epsilon)).or. &
            (eta .lt.(0.-epsilon)).or. &
            (zeta.lt.(0.-epsilon)).or. &
            (zeta.gt.(1.-xi-eta+epsilon)) ) then
                !
                XYZInElement = .false.
                !
        else
                !
                XYZInElement = .true.
                !
        endif
    
    case(6)!Hex

        VertexSide(1,:) =  (/ 1,2,6,5 /)   ! Local hex. vertices of side I        !
        VertexSide(2,:) =  (/ 2,4,8,6 /)   ! Local hex. vertices of side II       !
        VertexSide(3,:) =  (/ 4,3,7,8 /)   ! Local hex. vertices of side III      !
        VertexSide(4,:) =  (/ 3,1,5,7 /)   ! Local hex. vertices of side IV       !
        VertexSide(5,:) =  (/ 2,1,3,4 /)   ! Local hex. vertices of side V        !
        VertexSide(6,:) =  (/ 5,6,8,7 /)   ! Local hex. vertices of side VI       !

        inside = .true. 
        do iSide = 1, MESH%LocalElemType(iElem)
            vecP(1) = xP - MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem)) 
            vecP(2) = yP - MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem)) 
            vecP(3) = zP - MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,1),iElem)) 
            prod(iSide) = dot_product( vecP, DISC%Galerkin%geoNormals(:,iSide,iElem) )
            if(prod(iSide).gt.epsilon) then
                inside = .false. 
            else
                continue
            endif
        enddo

        XYZInElement = inside
        
    end select 

  end function XYZInElement


  !
  ! 3D Functions
  !



  ! ***************************************************************** !
  ! *                                                               * !
  ! * TrafoXiEta2XY maps xi-eta to xy coords. in 2D                 * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine TrafoXiEta2XY(xP,yP,xi,eta,x,y,eType)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    integer     :: eType                               ! Element type         !
    real        :: xP, yP                              ! Input:  xP, yP       !
    real        :: xi, eta, zeta                       ! Output: xi eta       !
    real        :: x(eType), y(eType)                  ! Vertex coordinates   !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xi,eta,eType
    intent(OUT) :: xP,yP
    !-------------------------------------------------------------------------!
    !
    if(eType.eq.3) then
      call TriTrafoXiEta2XY(xP,yP,xi,eta,x,y)
    elseif(eType.eq.4) then
      call QuadTrafoXiEta2XY(xP,yP,xi,eta,x,y)
    endif
    !
  end subroutine TrafoXiEta2XY


  ! ***************************************************************** !
  ! *                                                               * !
  ! * TrafoXY2XiEta maps xy to xi-eta coords. in 2D                 * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine TrafoXY2XiEta(xi,eta,xP,yP,x,y,eType)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    integer     :: eType                               ! Element type         !
    real        :: xP, yP                              ! Input:  xP, yP       !
    real        :: xi, eta, zeta                       ! Output: xi eta       !
    real        :: x(eType), y(eType)                  ! Vertex coordinates   !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xP,yP,eType
    intent(OUT) :: xi,eta
    !-------------------------------------------------------------------------!
    !
    if(eType.eq.3) then
      call TriTrafoXY2XiEta(xi,eta,xP,yP,x,y)
    elseif(eType.eq.4) then
      call QuadTrafoXY2XiEta(xi,eta,xP,yP,x,y)
    endif
    !
  end subroutine TrafoXY2XiEta

  ! ***************************************************************** !
  ! *                                                               * !
  ! * TriTrafoXiEta2XY maps xi-eta to xy coords. in triangles       * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine TriTrafoXiEta2XY(xP,yP,xi,eta,x,y)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    integer     :: eType                               ! Element type         !
    real        :: xP, yP                              ! Input:  xP, yP       !
    real        :: xi, eta, zeta                       ! Output: xi eta       !
    real        :: x(3), y(3)                          ! Vertex coordinates   !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xi,eta
    intent(OUT) :: xP,yP
    !-------------------------------------------------------------------------!
    !
    xP = x(1) + (x(2)-x(1))*xi + (x(3)-x(1))*eta
    yP = y(1) + (y(2)-y(1))*xi + (y(3)-y(1))*eta
    !
  end subroutine TriTrafoXiEta2XY


  ! ***************************************************************** !
  ! *                                                               * !
  ! * QuadTrafoXiEta2XY maps xi-eta to xy coords. in quadrilaterals * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine QuadTrafoXiEta2XY(xP,yP,xi,eta,x,y)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xP, yP                              ! Output: xP, yP       !
    real        :: xi, eta                             ! Input:  xi eta       !
    real        :: x(4), y(4)                          ! Vertex coordinates   !
    ! Local Variable declaration
    integer     :: i                                   ! Loop index           ! 
    real        :: xiU, etaU                           ! modified xi eta zeta !
    real        :: L01xi, L01eta                       ! Lagrangian Polynomial!
    real        :: L11xi, L11eta                       ! Lagrangian Polynomial!
    real        :: SF(4)                               ! Shape Function       !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xi,eta
    intent(OUT) :: xP,yP
    !-------------------------------------------------------------------------!
    !
    ! Compute Lagrangian Polynomials of degree 1 for 4-node quadrilateral (straight edges)
    !
    L01xi   = 1  -xi
    L11xi   =     xi
    L01eta  = 1 -eta
    L11eta  =    eta
    !
    ! Compute shape functions of 4-node quadrilateral
    !
    SF(1) = L01xi * L01eta 
    SF(2) = L11xi * L01eta 
    SF(3) = L11xi * L11eta 
    SF(4) = L01xi * L11eta 
    !
    ! Coordinate Transformation from xiU, etaU to X, Y
    !
    xP = 0.
    yP = 0.
    !
    do i = 1,4
        xP = xP + SF(i) * x(i)
        yP = yP + SF(i) * y(i)
    enddo
    !
  end subroutine QuadTrafoXiEta2XY


  ! ***************************************************************** !
  ! *                                                               * !
  ! * TriTrafoXY2XiEta maps xy to xi-eta coords. on triangles       * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine TriTrafoXY2XiEta(xi,eta,xP,yP,x,y)
    
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    type(tUnstructMesh) :: MESH
    real        :: xP, yP                              ! Input:  xP, yP       !
    real        :: xi, eta, zeta                       ! Output: xi eta       !
    real        :: x(3)                                !
    real        :: y(3)                                !
    ! Local variable declaration
    real        :: J                                   ! Jacobi determinant   !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xP,yP
    intent(OUT) :: xi,eta
    !-------------------------------------------------------------------------!
    !
    J   = x(2)*y(3)-x(2)*y(1)-x(1)*y(3)-x(3)*y(2)+x(3)*y(1)+x(1)*y(2)
    xi  = ( (x(3)*y(1)-x(1)*y(3)) + xP*(y(3)-y(1)) + yP*(x(1)-x(3)) )/J
    eta = ( (x(1)    *y(2)-x(2)*y(1)    ) + xP*(y(1)    -y(2)) + yP*(x(2)-x(1))     )/J
    !
  end subroutine TriTrafoXY2XiEta


  ! ***************************************************************** !
  ! *                                                               * !
  ! * QuadTrafoXY2XiEta maps xy to xi-eta coords. on triangles      * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine QuadTrafoXY2XiEta(xi,eta,xP,yP,x,y)

    use COMMON_operators_mod
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xP, yP                              ! Input:  xP, yP       !
    real        :: xi, eta                             ! Output: xi eta       !
    real        :: x(4)              !
    real        :: y(4)              !
    real        :: x2(3)              !
    real        :: y2(3)              !

    ! Local variable declaration
    ! Local Variable declaration
    integer     :: i                                   ! Loop index           ! 
    real        :: solution(2)                         ! modified xi eta      !
    real        :: xT,yT 
    logical     :: check

    real        :: xi_old, eta_old                     ! Iteration values for newton
    real        :: xp_old, yp_old                      ! Iteration values for newton
    real        :: grad(2,2), JacobiT(2,2)             ! dx/dxi   dxi/dx
    real        :: F(2), dX(2)
    real        :: TOL                                 ! Tolerance 
    real        :: ERROR, ERROR_xi, ERROR_eta
    integer     :: count                               ! counter
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xP,yP
    intent(OUT) :: xi,eta
    !-------------------------------------------------------------------------!
    !
    ! Compute the initial guess using the tetrahedral mapping
    !x2(1:2)=x(1:2)
    !x2(3)=x(4)

    !y2(1:2)=y(1:2)
    !y2(3)=y(4)

    !CALL TriTrafoXY2XiEta(xi,eta,xP,yP,x2,y2)
    !
    !solution(:) = (/ xi,eta /)
    ! Initialize goal function of the Newton algorithm
    !CALL iniFuncV2D(xP,yP,x,y)
    ! 
    !CALL newt(solution,2,check)
    !xi   = solution(1)
    !eta  = solution(2)
    !
    !CONTINUE

    ! Guess value baricenter
    ! Baricenter of the element
    ![xrk , yrk] = TriTrafoXY2XiEta(xp,yp,X(1:3),Y(1:3));
    xi_old = 0.5;   eta_old = 0.5;
    dX  = (/ 0, 0 /)
    
    Tol   = 1.0e-5;
    Error = 1.0e5;
    count = 0;
    do while (Error > Tol)
        count = count +1;
        call QuadTrafoXiEta2XY(xp_old,yp_old,xi_old,eta_old,x,y);

        !         dX/dxi
        ! | dx/dxi     dx/deta |
        ! | dy/dxi     dy/deta |
        call QuadTrafoXiEtaGrad(grad,xi_old,eta_old,x,y)
        grad = -grad

        !         dxi/dX     This cannot be computed directly for general quadrangular elements
        ! | dxi/dx     dxi/dy  |
        ! | deta/dx    deta/dy |
        call MatrixInverse2x2(JacobiT,grad)

        F = (/ xP - xp_old , yP - yp_old /);
        
        dX = - matmul(JacobiT,F);
        
        if (abs(xi_old)<Tol) then
            Error_xi = abs(dX(1));
        else
            Error_xi = abs(dX(1))/(0.5*(2*xi_old + dX(1)));
        endif
        if (abs(eta_old)<Tol) then
            Error_eta = abs(dX(2));
        else
            Error_eta = abs(dX(2))/(0.5*(2*eta_old + dX(2)));
        endif
        
        Error = abs(max(Error_xi,Error_eta));
        
        xi_old  = xi_old + dX(1);
        eta_old = eta_old + dX(2);
        if(count>20) then
            write(*,*) ' | Error in QuadTrafoXY2XiEta. '
            write(*,*) ' |  Inverse of Jacobian not founded.'
            call MPI_ABORT(MPI%commWorld, 134)
        endif
    enddo
    
    xi  = xi_old;
    eta = eta_old;

    ! Rounding the value with tolerance 1.0e10
    !aux = xr *  1.e10;
    !aux = round(aux);
    !xr  = aux * 1.e-10;
    
    !aux = yr *  1.e10;
    !aux = round(aux);
    !yr  = aux * 1.e-10;

    !
  end subroutine QuadTrafoXY2XiEta


  ! *************************************************************************************** !
  ! *                                                                                     * !
  ! * TetraTrafoXiEtaZeta2XYZ maps xi-eta-zeta coords. to xyz coords. in a tetrahedron    * !
  ! *                                                                                     * !
  ! *************************************************************************************** !

  subroutine TetraTrafoXiEtaZeta2XYZ(xP,yP,zP,xi,eta,zeta,x,y,z)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xP, yP, zP                          ! Input:  xP, yP, zP   !
    real        :: xi, eta, zeta                       ! Output: xi eta zeta  !
    real        :: x(4), y(4), z(4)                    ! Vertex coordinates   !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,z,xi,eta,zeta
    intent(OUT) :: xP,yP,zP
    !-------------------------------------------------------------------------!
    !
    xP = x(1) + (x(2)-x(1))*xi + (x(3)-x(1))*eta + (x(4)-x(1))*zeta
    yP = y(1) + (y(2)-y(1))*xi + (y(3)-y(1))*eta + (y(4)-y(1))*zeta
    zP = z(1) + (z(2)-z(1))*xi + (z(3)-z(1))*eta + (z(4)-z(1))*zeta 
    !
  end subroutine TetraTrafoXiEtaZeta2XYZ

  ! *************************************************************************************** !
  ! *                                                                                     * !
  ! * HexaTrafoXiEtaZeta2XYZ maps xi-eta-zeta coords. to xyz coords. in a hexahedron      * !
  ! *                                                                                     * !
  ! *************************************************************************************** !

  subroutine HexaTrafoXiEtaZeta2XYZ(xP,yP,zP,xi,eta,zeta,x,y,z)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xP, yP, zP                          ! Output: xP, yP, zP   !
    real        :: xi, eta, zeta                       ! Input:  xi eta zeta  !
    real        :: x(8), y(8), z(8)                    ! Vertex coordinates   !
    ! Local Variable declaration
    integer     :: i                                   ! Loop index           ! 
    real        :: xiU, etaU, zetaU                    ! modified xi eta zeta !
    real        :: L01xi, L01eta, L01zeta              ! Lagrangian Polynomial!
    real        :: L11xi, L11eta, L11zeta              ! Lagrangian Polynomial!
    real        :: SF(8)                               ! Shape Function       !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,z,xi,eta,zeta
    intent(OUT) :: xP,yP,zP
    !-------------------------------------------------------------------------!
    !
    ! Compute Lagrangian Polynomials of degree 1 for 8-node hexahedron (straight edges)
    !
    L01xi   = 1  -xi
    L11xi   =     xi
    L01eta  = 1 -eta
    L11eta  =    eta
    L01zeta = 1-zeta
    L11zeta =   zeta
    !
    ! Compute shape functions of 8-node hexahedron
    !
    SF(1) = L01xi * L01eta * L01zeta
    SF(2) = L11xi * L01eta * L01zeta
    SF(3) = L01xi * L11eta * L01zeta
    SF(4) = L11xi * L11eta * L01zeta
    SF(5) = L01xi * L01eta * L11zeta
    SF(6) = L11xi * L01eta * L11zeta
    SF(7) = L01xi * L11eta * L11zeta
    SF(8) = L11xi * L11eta * L11zeta
    !
    ! Coordinate Transformation from xiU, etaU, zetaU to X, Y, Z
    !
    xP = 0.
    yP = 0.
    zP = 0.
    !
    do i = 1,8
        xP = xP + SF(i) * x(i)
        yP = yP + SF(i) * y(i)
        zP = zP + SF(i) * z(i)
    enddo
    !
  end subroutine HexaTrafoXiEtaZeta2XYZ


  ! ********************************************************************************** !
  ! *                                                                                * !
  ! * TrafoXYZ2XiEtaZeta maps xyz coords. to xi, eta, zeta  coordinates              * !
  ! *                                                                                * !
  ! ********************************************************************************** !

  subroutine TrafoXYZ2XiEtaZeta(xi,eta,zeta,xP,yP,zP,x,y,z,vType)
    
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    integer     :: vType
    real        :: xi, eta, zeta                       ! Output: xi eta zeta  !
    real        :: xP, yP, zP                          ! Input:  xP, yP, zP   !
    real        :: x(vType)                            !
    real        :: y(vType)                            !
    real        :: z(vType)                            !     
    ! Local variable declaration
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,z,xP,yP,zP,vType
    intent(OUT) :: xi,eta,zeta
    !-------------------------------------------------------------------------!
    !
!    select case(vType)
!    case(4)
        call TetraTrafoXYZ2XiEtaZeta(xi,eta,zeta,xP,yP,zP,x,y,z,vType)
!    case(8)
      !  call HexaTrafoXYZ2XiEtaZeta(xi,eta,zeta,xP,yP,zP,x,y,z,vType)
      ! 'Hexas not supported at the moment!'
!    end select
    !
  end subroutine TrafoXYZ2XiEtaZeta
  
  ! ********************************************************************************** !
  ! *                                                                                * !
  ! * TetraTrafoXYZ2XiEtaZeta maps xyz coords. to xi...  coords. in a tetrahedron    * !
  ! *                                                                                * !
  ! ********************************************************************************** !

  subroutine TetraTrafoXYZ2XiEtaZeta(xi,eta,zeta,xP,yP,zP,x,y,z,vType)
    
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    integer     :: vType
    real        :: xi, eta, zeta                       ! Output: xi eta zeta  !
    real        :: xP, yP, zP                          ! Input:  xP, yP, zP   !
    real        :: x(vType)                            !
    real        :: y(vType)                            !
    real        :: z(vType)                            !     
    ! Local variable declaration
    real        :: J                                   ! Jacobi determinant   !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,z,xP,yP,zP,vType
    intent(OUT) :: xi,eta,zeta
    !-------------------------------------------------------------------------!
    !
    select case(vType)
    case(4)
        J = -z(1)*x(2)*y(3)+x(1)*z(3)*y(4)+x(1)*z(2)*y(3)+z(4)*x(2)*y(3)+             &
             z(3)*x(4)*y(2)+z(2)*y(4)*x(3)+z(4)*y(1)*x(3)-z(3)*y(1)*x(4)+z(3)*        &
             y(1)*x(2)-x(1)*z(2)*y(4)+x(1)*z(4)*y(2)+z(1)*x(2)*y(4)+z(2)*y(1)*        &
             x(4)-z(4)*y(1)*x(2)-x(1)*z(4)*y(3)-z(3)*y(4)*x(2)-z(1)*x(3)*y(4)-        &
             z(4)*x(3)*y(2)-z(2)*x(4)*y(3)+z(1)*x(3)*y(2)-z(2)*y(1)*x(3)-z(1)*x(4)*   &
             y(2)-x(1)*z(3)*y(2)+z(1)*x(4)*y(3)
        !
        xi = ( z(4)*y(3)+z(1)*y(4)-z(3)*y(4)-z(1)*y(3)-z(4)*y(1)+z(3)*y(1))/J*xP +    &
             (-x(1)*z(3)+z(4)*x(1)-z(4)*x(3)+z(1)*x(3)-z(1)*x(4)+x(4)*z(3))/J*yP +    &
             (-y(1)*x(3)+y(3)*x(1)+y(1)*x(4)+y(4)*x(3)-y(3)*x(4)-y(4)*x(1))/J*zP +    &
             ( x(1)*z(3)*y(4)-z(3)*y(1)*x(4)+z(1)*x(4)*y(3)-x(1)*z(4)*y(3) -          &
               z(1)*x(3)*y(4)+z(4)*y(1)*x(3))/J
        !
        eta = -( z(4)*y(2)-z(1)*y(2)+z(1)*y(4)+z(2)*y(1)-z(2)*y(4)-z(4)*y(1))/J*xP -  &
               ( z(2)*x(4)-z(2)*x(1)-z(1)*x(4)-z(4)*x(2)+z(4)*x(1)+z(1)*x(2))/J*yP -  &
               (-y(2)*x(4)+y(1)*x(4)+y(4)*x(2)-y(1)*x(2)-y(4)*x(1)+y(2)*x(1))/J*zP -  &
               (-z(2)*y(1)*x(4)+x(1)*z(2)*y(4)-x(1)*z(4)*y(2)-z(1)*x(2)*y(4) +        &
                 z(1)*x(4)*y(2)+z(4)*y(1)*x(2))/J
        !
        zeta = (z(1)*y(3)-z(3)*y(1)-z(1)*y(2)+z(3)*y(2)+z(2)*y(1)-z(2)*y(3))/J*xP +   &
               (x(1)*z(3)-z(2)*x(1)-z(3)*x(2)+z(1)*x(2)+z(2)*x(3)-z(1)*x(3))/J*yP +   &
               (y(2)*x(1)-y(3)*x(1)+x(2)*y(3)-y(1)*x(2)+y(1)*x(3)-x(3)*y(2))/J*zP +   &
               (-x(1)*z(3)*y(4)-z(4)*x(2)*y(3)+z(1)*x(3)*y(4)-z(3)*x(4)*y(2) -        &
                 z(2)*y(4)*x(3)-z(4)*y(1)*x(3)+z(3)*y(1)*x(4)-z(2)*y(1)*x(4) +        &
                 x(1)*z(2)*y(4)-x(1)*z(4)*y(2)-z(1)*x(2)*y(4)+z(1)*x(4)*y(2) +        &
                 z(4)*y(1)*x(2)+x(1)*z(4)*y(3)+z(3)*y(4)*x(2)+z(4)*x(3)*y(2) -        &
                 z(1)*x(4)*y(3)+z(2)*x(4)*y(3)+J)/J
    case(8) ! Used as initial guess for mapping in hexas
        J = -z(1)*x(2)*y(3)+x(1)*z(3)*y(5)+x(1)*z(2)*y(3)+z(5)*x(2)*y(3)+             &
             z(3)*x(5)*y(2)+z(2)*y(5)*x(3)+z(5)*y(1)*x(3)-z(3)*y(1)*x(5)+z(3)*        &
             y(1)*x(2)-x(1)*z(2)*y(5)+x(1)*z(5)*y(2)+z(1)*x(2)*y(5)+z(2)*y(1)*        &
             x(5)-z(5)*y(1)*x(2)-x(1)*z(5)*y(3)-z(3)*y(5)*x(2)-z(1)*x(3)*y(5)-        &
             z(5)*x(3)*y(2)-z(2)*x(5)*y(3)+z(1)*x(3)*y(2)-z(2)*y(1)*x(3)-z(1)*x(5)*   &
             y(2)-x(1)*z(3)*y(2)+z(1)*x(5)*y(3)
        !
        xi = ( z(5)*y(3)+z(1)*y(5)-z(3)*y(5)-z(1)*y(3)-z(5)*y(1)+z(3)*y(1))/J*xP +    &
             (-x(1)*z(3)+z(5)*x(1)-z(5)*x(3)+z(1)*x(3)-z(1)*x(5)+x(5)*z(3))/J*yP +    &
             (-y(1)*x(3)+y(3)*x(1)+y(1)*x(5)+y(5)*x(3)-y(3)*x(5)-y(5)*x(1))/J*zP +    &
             ( x(1)*z(3)*y(5)-z(3)*y(1)*x(5)+z(1)*x(5)*y(3)-x(1)*z(5)*y(3) -          &
               z(1)*x(3)*y(5)+z(5)*y(1)*x(3))/J
        !
        eta = -( z(5)*y(2)-z(1)*y(2)+z(1)*y(5)+z(2)*y(1)-z(2)*y(5)-z(5)*y(1))/J*xP -  &
               ( z(2)*x(5)-z(2)*x(1)-z(1)*x(5)-z(5)*x(2)+z(5)*x(1)+z(1)*x(2))/J*yP -  &
               (-y(2)*x(5)+y(1)*x(5)+y(5)*x(2)-y(1)*x(2)-y(5)*x(1)+y(2)*x(1))/J*zP -  &
               (-z(2)*y(1)*x(5)+x(1)*z(2)*y(5)-x(1)*z(5)*y(2)-z(1)*x(2)*y(5) +        &
                 z(1)*x(5)*y(2)+z(5)*y(1)*x(2))/J
        !
        zeta = (z(1)*y(3)-z(3)*y(1)-z(1)*y(2)+z(3)*y(2)+z(2)*y(1)-z(2)*y(3))/J*xP +   &
               (x(1)*z(3)-z(2)*x(1)-z(3)*x(2)+z(1)*x(2)+z(2)*x(3)-z(1)*x(3))/J*yP +   &
               (y(2)*x(1)-y(3)*x(1)+x(2)*y(3)-y(1)*x(2)+y(1)*x(3)-x(3)*y(2))/J*zP +   &
               (-x(1)*z(3)*y(5)-z(5)*x(2)*y(3)+z(1)*x(3)*y(5)-z(3)*x(5)*y(2) -        &
                 z(2)*y(5)*x(3)-z(5)*y(1)*x(3)+z(3)*y(1)*x(5)-z(2)*y(1)*x(5) +        &
                 x(1)*z(2)*y(5)-x(1)*z(5)*y(2)-z(1)*x(2)*y(5)+z(1)*x(5)*y(2) +        &
                 z(5)*y(1)*x(2)+x(1)*z(5)*y(3)+z(3)*y(5)*x(2)+z(5)*x(3)*y(2) -        &
                 z(1)*x(5)*y(3)+z(2)*x(5)*y(3)+J)/J
    end select
    !
  end subroutine TetraTrafoXYZ2XiEtaZeta

  ! ********************************************************************************* !
  ! *                                                                               * !
  ! * HexaTrafoXYZ2XiEtaZeta maps xyz to xi-eta-zeta  coords. in a hexahedron       * !
  ! *                                                                               * !
  ! ********************************************************************************* !

!  subroutine HexaTrafoXYZ2XiEtaZeta(xi,eta,zeta,xP,yP,zP,x,y,z,vType)
!
!    use minoopack, only: hybrd_type
!
!    !-------------------------------------------------------------------------!
!    implicit none
!    !-------------------------------------------------------------------------!
!    ! Argument list declaration
!    integer     :: vType
!    real        :: xP, yP, zP                          ! Output: xP, yP, zP   !
!    real        :: xi, eta, zeta                       ! Input:  xi eta zeta  !
!    real        :: x(8), y(8), z(8)                    ! Vertex coordinates   !
!    ! Local Variable declaration
!    integer     :: i                                   ! Loop index           !
!    real        :: solution(3)                         ! modified xi eta zeta !
!    real        :: xT,yT,zT
!    logical     :: check
!    !variables for minoopack solution:
!    real :: tol,fvec(3)
!    integer :: info
!    type(hybrd_type):: hybrd
!    !-------------------------------------------------------------------------!
!    intent(IN)  :: x,y,z,xP,yP,zP,vType
!    intent(OUT) :: xi,eta,zeta
!    !-------------------------------------------------------------------------!
!
!    ! Compute the initial guess using the tetrahedral mapping
!    call TetraTrafoXYZ2XiEtaZeta(xi,eta,zeta,xP,yP,zP,x,y,z,vType)
!
!
!    solution = [xi,eta,zeta]; tol = 0.00001D+00
!    call hybrd%init(myf,pars=[x,y,z,xP,yP,zP])
!    call hybrd%solve(solution,fvec,tol,info)
!
!    !-----------------------------------------------------------------------
!    xi  =solution(1)
!    eta =solution(2)
!    zeta=solution(3)
!
!  contains
!    subroutine myf(slf,x,fvec,iflag,fjac)
!      use minoopack
!      implicit none
!      class(hybrd_type)::slf
!      real :: x(:),fvec(size(x))
!      integer :: iflag
!      real, optional :: fjac(size(x),size(x))
!      real :: xi,eta,zeta
!      real :: xP,yP,zP
!      xi = x(1); eta = x(2); zeta = x(3)
!      call HexaTrafoXiEtaZeta2XYZ(xP,yP,zP,xi,eta,zeta,slf%pars(1:8),slf%pars(9:16),slf%pars(17:24))
!      fvec(1:3) = (slf%pars(25:27) - [xP,yP,zP]) !**2
!    end subroutine myf
!
!  end subroutine HexaTrafoXYZ2XiEtaZeta

  ! ********************************************************************************* !
  ! *                                                                               * !
  ! * HexaTrafoXiEtaZetaGrad computes the gradient of the hexahedral mapping        * !
  ! *                                                                               * !
  ! ********************************************************************************* !

  subroutine HexaTrafoXiEtaZetaGrad(grad,xi,eta,zeta,x,y,z)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: grad(3,3)                           ! Output:  grad        !
    real        :: xi, eta, zeta                       ! Input: xi eta zeta   !
    real        :: x(8), y(8), z(8)                    ! Vertex coordinates   !
    ! Local Variable declaration
    integer     :: i,j                                 ! Loop index           ! 
    real        :: L01xi, L01eta, L01zeta              ! Lagrangian Polynomial!
    real        :: L11xi, L11eta, L11zeta              ! Lagrangian Polynomial!
    real        :: SF(8,0:3)                           ! Shape Function       !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,z,xi,eta,zeta
    intent(OUT) :: grad
    !-------------------------------------------------------------------------!
    !
    ! Compute Lagrangian Polynomials of degree 1 for 8-node hexahedron (straight edges)
    !
    L01xi   = 1  -xi
    L11xi   =     xi
    L01eta  = 1 -eta
    L11eta  =    eta
    L01zeta = 1-zeta
    L11zeta =   zeta
    !
    ! Compute the shape functions of an 8-node hexahedron
    !
    SF(1,0) = L01xi * L01eta * L01zeta
    SF(2,0) = L11xi * L01eta * L01zeta
    SF(3,0) = L01xi * L11eta * L01zeta
    SF(4,0) = L11xi * L11eta * L01zeta
    SF(5,0) = L01xi * L01eta * L11zeta
    SF(6,0) = L11xi * L01eta * L11zeta
    SF(7,0) = L01xi * L11eta * L11zeta
    SF(8,0) = L11xi * L11eta * L11zeta
    !
    ! Compute the xi-derivative of the shape functions of an 8-node hexahedron
    !
    SF(1,1) = (-1) * L01eta * L01zeta
    SF(2,1) =        L01eta * L01zeta
    SF(3,1) = (-1) * L11eta * L01zeta
    SF(4,1) =        L11eta * L01zeta
    SF(5,1) = (-1) * L01eta * L11zeta
    SF(6,1) =        L01eta * L11zeta
    SF(7,1) = (-1) * L11eta * L11zeta
    SF(8,1) =        L11eta * L11zeta
    !
    ! Compute the eta-derivative of the shape functions of an 8-node hexahedron
    !
    SF(1,2) = L01xi * (-1) * L01zeta
    SF(2,2) = L11xi * (-1) * L01zeta
    SF(3,2) = L01xi        * L01zeta
    SF(4,2) = L11xi        * L01zeta
    SF(5,2) = L01xi * (-1) * L11zeta
    SF(6,2) = L11xi * (-1) * L11zeta
    SF(7,2) = L01xi        * L11zeta
    SF(8,2) = L11xi        * L11zeta
    !
    ! Compute the zeta-derivative of the shape functions of an 8-node hexahedron
    !
    SF(1,3) = L01xi * L01eta * (-1)
    SF(2,3) = L11xi * L01eta * (-1)
    SF(3,3) = L01xi * L11eta * (-1)
    SF(4,3) = L11xi * L11eta * (-1)
    SF(5,3) = L01xi * L01eta  
    SF(6,3) = L11xi * L01eta  
    SF(7,3) = L01xi * L11eta  
    SF(8,3) = L11xi * L11eta  
    !
    ! Gradient of the coordinate Transformation from xi, eta, zeta to X, Y, Z
    !
    grad(:,:) = 0.
    !
    do j = 1, 3
     do i = 1,8
        grad(1,j) = grad(1,j) + SF(i,j) * x(i)
        grad(2,j) = grad(2,j) + SF(i,j) * y(i)
        grad(3,j) = grad(3,j) + SF(i,j) * z(i)
     enddo
    enddo
    !
  end subroutine HexaTrafoXiEtaZetaGrad

  ! ********************************************************************************* !
  ! *                                                                               * !
  ! * TetraTrafoXiEtaZetaGrad computes the gradient of the tetrahedral mapping      * !
  ! *  which is constant inside element                                             * !
  ! *                                                                               * !
  ! ********************************************************************************* !

  subroutine TetraTrafoXiEtaZetaGrad(grad,x,y,z)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: grad(3,3)                           ! Output:  grad        !
    real        :: x(4), y(4), z(4)                    ! Vertex coordinates   !
    ! Local Variable declaration
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,z
    intent(OUT) :: grad
    !-------------------------------------------------------------------------!
    !
    grad(1,1) = x(2)-x(1); grad(1,2) = x(3)-x(1); grad(1,3) = x(4)-x(1)
    grad(2,1) = y(2)-y(1); grad(2,2) = y(3)-y(1); grad(2,3) = y(4)-y(1)
    grad(3,1) = z(2)-z(1); grad(3,2) = z(3)-z(1); grad(3,3) = z(4)-z(1)
    !
  end subroutine TetraTrafoXiEtaZetaGrad

  ! ********************************************************************************* !
  ! *                                                                               * !
  ! * QuadTrafoXiEtaGrad computes the gradient of the quadrangular mapping          * !
  ! *                                                                               * !
  ! ********************************************************************************* !

  subroutine QuadTrafoXiEtaGrad(grad,xi,eta,x,y)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: grad(2,2)                           ! Output:  grad        !
    real        :: xi, eta                             ! Input: xi eta        !
    real        :: x(4), y(4)                          ! Vertex coordinates   !
    ! Local Variable declaration
    integer     :: i,j                                 ! Loop index           ! 
    real        :: L01xi, L01eta                       ! Lagrangian Polynomial!
    real        :: L11xi, L11eta                       ! Lagrangian Polynomial!
    real        :: SF(4,0:2)                           ! Shape Function       !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xi,eta
    intent(OUT) :: grad
    !-------------------------------------------------------------------------!
    !
    ! Compute Lagrangian Polynomials of degree 1 for 4-node square (straight edges)
    !
    L01xi   = 1  -xi
    L11xi   =     xi
    L01eta  = 1 -eta
    L11eta  =    eta
    !
    ! Compute the shape functions of an 4-node square
    !
    SF(1,0) = L01xi * L01eta
    SF(2,0) = L11xi * L01eta
    SF(3,0) = L11xi * L11eta
    SF(4,0) = L01xi * L11eta
    !
    ! Compute the xi-derivative of the shape functions of an 4-node square
    !
    SF(1,1) = (-1) * L01eta
    SF(2,1) =        L01eta
    SF(3,1) =        L11eta
    SF(4,1) = (-1) * L11eta
    !
    ! Compute the eta-derivative of the shape functions of an 4-node square
    !
    SF(1,2) = L01xi * (-1)
    SF(2,2) = L11xi * (-1)
    SF(3,2) = L11xi
    SF(4,2) = L01xi
    !
    ! Gradient of the coordinate Transformation from xi, eta  to X, Y 
    !
    !         dX/dxi
    ! | dx/dxi     dx/deta |
    ! | dy/dxi     dy/deta |
    !
    grad(:,:) = 0.
    !
    do j = 1, 2
     do i = 1,4
        grad(1,j) = grad(1,j) + SF(i,j) * x(i)
        grad(2,j) = grad(2,j) + SF(i,j) * y(i)
     enddo
    enddo
    !
  end subroutine QuadTrafoXiEtaGrad


  ! ********************************************************************************* !
  ! *                                                                               * !
  ! * TriTrafoXiEtaGrad computes the gradient of the triangular mapping             * !
  ! *                                                                               * !
  ! ********************************************************************************* !

  subroutine TriTrafoXiEtaGrad(grad,xi,eta,x,y)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: grad(2,2)                           ! Output:  grad        !
    real        :: xi, eta                             ! Input: xi eta        !
    real        :: x(3), y(3)                          ! Vertex coordinates   !
    ! Local Variable declaration
     !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,xi,eta
    intent(OUT) :: grad
    !-------------------------------------------------------------------------!
    !
    !         dX/dxi
    ! | dx/dxi     dx/deta |
    ! | dy/dxi     dy/deta |
    grad(1,1) = x(2)     - x(1)
    grad(1,2) = x(3)     - x(1)
    grad(2,1) = y(2)     - y(1)
    grad(2,2) = y(3)     - y(1)

  end subroutine TriTrafoXiEtaGrad
  
  
  ! ***************************************************************** !
  ! *                                                               * !
  ! * TriTrafoChi2XiEta maps side-local chi coordinates             * ! 
  ! * to reference xi-eta coordinates on triangles                  * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine TriTrafoChi2XiEta(xi,eta,chi,iSide,iNeighborSide)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xi, eta                             ! Output: xi eta       !
    real        :: chi                                 ! Input:  chi          !
    integer     :: iSide, iNeighborSide                ! Element orientation  !
    ! Local variable declaration
    real        :: chi1                                ! Intermediate results !
    !-------------------------------------------------------------------------!
    intent(IN)  :: chi, iSide, iNeighborSide
    intent(OUT) :: xi,eta
    !-------------------------------------------------------------------------!
    !
    select case(iNeighborSide)
    ! Inside the element itself
    case(0) 
        select case(iSide)
        case(1)
          xi   = chi
          eta  = 0.
        case(2)
          xi   = 1.-chi
          eta  = chi 
        case(3)
          xi   = 0. 
          eta  = 1.-chi
        end select
        !
        return
    ! Inside the neighbor
    case(1)
      chi1 = 1.-chi
      xi   = chi1
      eta  = 0.
    case(2)
      chi1 = 1.-chi
      xi   = 1.-chi1
      eta  = chi1 
    case(3)
      chi1 = 1.-chi
      xi   = 0. 
      eta  = 1.-chi1
    end select
    !
    return
    !
  end subroutine TriTrafoChi2XiEta

  ! ***************************************************************** !
  ! *                                                               * !
  ! * QuadTrafoChi2XiEta maps side-local chi coordinates             * ! 
  ! * to reference xi-eta coordinates for quadrangular elements     * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine QuadTrafoChi2XiEta(xi,eta,chi,iSide,iNeighborSide)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xi, eta                             ! Output: xi eta       !
    real        :: chi                                 ! Input:  chi          !
    integer     :: iSide, iNeighborSide                ! Element orientation  !
    ! Local variable declaration
    real        :: chi1                                ! Intermediate results !
    !-------------------------------------------------------------------------!
    intent(IN)  :: chi, iSide, iNeighborSide
    intent(OUT) :: xi,eta
    !-------------------------------------------------------------------------!
    !
    select case(iNeighborSide)
    ! Inside the element itself
    case(0) 
        select case(iSide)
        case(1)
          xi   = chi
          eta  = 0.
        case(2)
          xi   = 1.
          eta  = chi 
        case(3)
          xi   = 1.-chi
          eta  = 1.
        case(4)
          xi   = 0.
          eta  = 1.-chi
        end select
        !
        return
    ! Inside the neighbor
    case(1)
      chi1 = 1.-chi
      xi   = chi1
      eta  = 0.
    case(2)
      chi1 = 1.-chi
      xi   = 1.
      eta  = chi1
    case(3)
      chi1 = 1.-chi
      xi   = 1.-chi1
      eta  = 1.
    case(4)
      chi1 = 1.-chi
      xi   = 0.
      eta  = 1.-chi1
    end select
    !
    return
    !
  end subroutine QuadTrafoChi2XiEta

  ! ***************************************************************** !
  ! *                                                               * !
  ! * TrafoChiTau2XiEtaZeta maps side-local chi-tau coordinates     * ! 
  ! * to reference xi-eta-zeta coordinates                          * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,iNeighborVertex)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xi, eta, zeta                       ! Output: xi eta zeta  !
    real        :: chi,tau                             ! Input:  chi, tau     !
    integer     :: iSide, iNeighborVertex              ! Element orientation  !
    ! Local variable declaration
    real        :: chi1,tau1                           ! Intermediate results !
    !-------------------------------------------------------------------------!
    intent(IN)  :: chi, tau, iSide, iNeighborVertex
    intent(OUT) :: xi,eta,zeta
    !-------------------------------------------------------------------------!
    !
    select case(iNeighborVertex)
    ! Inside the element itself
    case(0) 
        select case(iSide)
        case(1)
          xi   = tau
          eta  = chi
          zeta = 0. 
        case(2)
          xi   = chi
          eta  = 0. 
          zeta = tau
        case(3)
          xi   = 0. 
          eta  = tau
          zeta = chi
        case(4)
          xi   = 1.-chi-tau
          eta  = chi
          zeta = tau
        end select
        !
        return
    ! Inside the neighbor
    case(1)
      chi1 = tau
      tau1 = chi
    case(2)
      chi1 = 1.-chi-tau
      tau1 = tau
    case(3)
      chi1 = chi
      tau1 = 1.-chi-tau
    end select
    !
    select case(iSide)
    case(1)
      xi   = tau1
      eta  = chi1
      zeta = 0. 
    case(2)
      xi   = chi1
      eta  = 0. 
      zeta = tau1
    case(3)
      xi   = 0. 
      eta  = tau1
      zeta = chi1
    case(4)
      xi   = 1.-chi1-tau1
      eta  = chi1
      zeta = tau1
    end select
    !
    return
    !
  end subroutine TrafoChiTau2XiEtaZeta


  ! ***************************************************************** !
  ! *                                                               * !
  ! * HexaTrafoChiTau2XiEtaZeta maps side-local chi-tau coordinates * ! 
  ! * to reference xi-eta-zeta coordinates                          * !
  ! *                                                               * !
  ! ***************************************************************** !

  subroutine HexaTrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,iNeighborVertex)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: xi, eta, zeta                       ! Output: xi eta zeta  !
    real        :: chi,tau                             ! Input:  chi, tau     !
    integer     :: iSide, iNeighborVertex              ! Element orientation  !
    ! Local variable declaration
    real        :: chi1,tau1                           ! Intermediate results !
    !-------------------------------------------------------------------------!
    intent(IN)  :: chi, tau, iSide, iNeighborVertex
    intent(OUT) :: xi,eta,zeta
    !-------------------------------------------------------------------------!
    !
    select case(iNeighborVertex)
    case(0) ! Inside the element itself        
        select case(iSide)
        case(1)
          xi   = chi
          eta  = 0.
          zeta = tau 
        case(2)
          xi   = 1. 
          eta  = chi 
          zeta = tau
        case(3)
          xi   = 1.-chi 
          eta  = 1.
          zeta = tau
        case(4)
          xi   = 0. 
          eta  = 1.-chi
          zeta = tau
        case(5)
          xi   = 1.-chi
          eta  = tau 
          zeta = 0.
        case(6)
          xi   = chi
          eta  = tau
          zeta = 1. 
        end select
        !
        return
    case(1) ! Inside the neighbor
      chi1 = tau
      tau1 = chi
    case(2)
      chi1 = 1.-chi
      tau1 = tau
    case(3)
      chi1 = chi
      tau1 = 1.-tau
    case(4)
      chi1 = 1.-tau
      tau1 = 1.-chi      
    end select
    !
    select case(iSide)
    case(1)
      xi   = chi1
      eta  = 0.
      zeta = tau1 
    case(2)
      xi   = 1. 
      eta  = chi1 
      zeta = tau1
    case(3)
      xi   = 1.-chi1 
      eta  = 1.
      zeta = tau1
    case(4)
      xi   = 0. 
      eta  = 1.-chi1
      zeta = tau1
    case(5)
      xi   = 1.-chi1
      eta  = tau1 
      zeta = 0.
    case(6)
      xi   = chi1
      eta  = tau1
      zeta = 1. 
    end select
    !
    return
    !
  end subroutine HexaTrafoChiTau2XiEtaZeta

  subroutine SpaceTimeBasis(iDegFr_xi,iDegFr_tau,iDegFrST,nDegFr)
    !-------------------------------------------------------------------------!
    
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    integer         :: iDegFr_xi,iDegFr_tau,iDegFrST,nDegFr
    !-------------------------------------------------------------------------!
    intent(IN)      :: iDegFrST,nDegFr
    intent(OUT) :: iDegFr_xi,iDegFr_tau
    !-------------------------------------------------------------------------!
    ! Input is the number of the space-time degree of freedom
    ! Output is the pair of numbers of the space and the time degrees of freedom
    ! for the tensor-product approach
    iDegFr_tau = ceiling(real(iDegFrST)/real(nDegFr))
    iDegFr_xi  = iDegFrST - (nDegFr)*(iDegFr_tau-1) 
  end subroutine SpaceTimeBasis

  ! ********************************************************************************* !
  ! *                                                                               * !
  ! * HexaTrafoXiEtaZetaGradGLL computes the gradient of the hexahedral mapping     * !
  ! *   for the [-1,+1]^3 reference element of GLL-based shape functinons           * !
  ! *                                                                               * !
  ! ********************************************************************************* !

  subroutine HexaTrafoXiEtaZetaGradGLL(grad,xi,eta,zeta,x,y,z)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: grad(3,3)                           ! Output:  grad        !
    real        :: xi, eta, zeta                       ! Input: xi eta zeta   !
    real        :: x(8), y(8), z(8)                    ! Vertex coordinates   !
    ! Local Variable declaration
    integer     :: i,j                                 ! Loop index           ! 
    real        :: L01xi, L01eta, L01zeta              ! Lagrangian Polynomial!
    real        :: L11xi, L11eta, L11zeta              ! Lagrangian Polynomial!
    real        :: SF(8,0:3)                           ! Shape Function       !
    !-------------------------------------------------------------------------!
    intent(IN)  :: x,y,z,xi,eta,zeta
    intent(OUT) :: grad
    !-------------------------------------------------------------------------!
    !
    ! Compute Lagrangian Polynomials of degree 1 for 8-node hexahedron (straight edges)
    !
    L01xi   = (1 - xi)   * 0.5d0
    L11xi   = (1 + xi)   * 0.5d0
    L01eta  = (1 - eta)  * 0.5d0
    L11eta  = (1 + eta)  * 0.5d0
    L01zeta = (1 - zeta) * 0.5d0
    L11zeta = (1 + zeta) * 0.5d0
    !
    ! Compute the shape functions of an 8-node hexahedron
    !
    SF(1,0) = L01xi * L01eta * L01zeta
    SF(2,0) = L11xi * L01eta * L01zeta
    SF(3,0) = L01xi * L11eta * L01zeta
    SF(4,0) = L11xi * L11eta * L01zeta
    SF(5,0) = L01xi * L01eta * L11zeta
    SF(6,0) = L11xi * L01eta * L11zeta
    SF(7,0) = L01xi * L11eta * L11zeta
    SF(8,0) = L11xi * L11eta * L11zeta
    !
    ! Compute the xi-derivative of the shape functions of an 8-node hexahedron
    !
    SF(1,1) = (-1) * L01eta * L01zeta
    SF(2,1) =        L01eta * L01zeta
    SF(3,1) = (-1) * L11eta * L01zeta
    SF(4,1) =        L11eta * L01zeta
    SF(5,1) = (-1) * L01eta * L11zeta
    SF(6,1) =        L01eta * L11zeta
    SF(7,1) = (-1) * L11eta * L11zeta
    SF(8,1) =        L11eta * L11zeta
    !
    ! Compute the eta-derivative of the shape functions of an 8-node hexahedron
    !
    SF(1,2) = L01xi * (-1) * L01zeta
    SF(2,2) = L11xi * (-1) * L01zeta
    SF(3,2) = L01xi        * L01zeta
    SF(4,2) = L11xi        * L01zeta
    SF(5,2) = L01xi * (-1) * L11zeta
    SF(6,2) = L11xi * (-1) * L11zeta
    SF(7,2) = L01xi        * L11zeta
    SF(8,2) = L11xi        * L11zeta
    !
    ! Compute the zeta-derivative of the shape functions of an 8-node hexahedron
    !
    SF(1,3) = L01xi * L01eta * (-1)
    SF(2,3) = L11xi * L01eta * (-1)
    SF(3,3) = L01xi * L11eta * (-1)
    SF(4,3) = L11xi * L11eta * (-1)
    SF(5,3) = L01xi * L01eta  
    SF(6,3) = L11xi * L01eta  
    SF(7,3) = L01xi * L11eta  
    SF(8,3) = L11xi * L11eta  
    !
    SF(:,:) = SF(:,:) * 0.5d0

    ! Gradient of the coordinate Transformation from xi, eta, zeta to X, Y, Z
    !
    grad(:,:) = 0.
    !
    do j = 1, 3
     do i = 1,8
        grad(1,j) = grad(1,j) + SF(i,j) * x(i)
        grad(2,j) = grad(2,j) + SF(i,j) * y(i)
        grad(3,j) = grad(3,j) + SF(i,j) * z(i)
     enddo
    enddo
    !
  end subroutine HexaTrafoXiEtaZetaGradGLL

  ! Compute the angle between two vectors
  !  from v1 to v2 in radians in the range [0 , 2 pi]
  ! Positive in Counterclockwise 
  subroutine Angle_Vector_Vector_2D(theta,v1,v2)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: theta                               ! Output:  angle       !
    real        :: v1(2), v2(2)                        ! Input :  vectors     !
    ! Local Variable declaration
    real        :: v1_aux(2), v2_aux(2)                ! auxiliary vectors    !
    real        :: v1_mod, v2_mod                      ! module of vectors    !
    real        :: theta1, theta2                      ! angles to the x axe  !
    real        :: cos_theta1, cos_theta2              ! cos(theta)           !
    real, parameter :: pi=3.141592653589793            ! CONSTANT pi
    !-------------------------------------------------------------------------!
    intent(IN)  :: v1,v2
    intent(OUT) :: theta
    !-------------------------------------------------------------------------!
    !
    ! Normalize the vectors
    v1_mod = sqrt(v1(1)**2+v1(2)**2);
    v1_aux = v1 / v1_mod;
    v2_mod = sqrt(v2(1)**2+v2(2)**2);
    v2_aux = v2 / v2_mod;
        
    cos_theta1 = (v1(1));   ! cos(theta) of vector 1 respect to [1,0]
    cos_theta2 = (v2(1));   ! cos(theta) of vector 2 respect to [1,0]
    
    theta1 = acos(cos_theta1);
    theta2 = acos(cos_theta2);
    
    ! If the rotation is higher than (pi rad) 2*pi-theta
    if (v1(2)<0) theta1 = 2*pi-theta1;
    if (v2(2)<0) theta2 = 2*pi-theta2;
    
    theta = theta2-theta1;

  end subroutine Angle_Vector_Vector_2D


  ! Compute the point location p mirrored at a 2D line
  ! given by the base point a and vector b

  subroutine Mirror_P_atLine(p_new,p_old,a,b)
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real        :: p_new(2)                                  ! Output:  vector      !
    real        :: p_old(2), a(2), b(2)                ! Input :  vectors     !
    ! Local Variable declaration
    real        :: t,d(2)                              ! auxiliary vectors    !
    !-------------------------------------------------------------------------!
    intent(IN)  :: p_old, a, b
    intent(OUT) :: p_new
    !-------------------------------------------------------------------------!
    !
    ! mirror point at line               
    
    t    = ((p_old(1) - a(1)) * b(1) + (p_old(2) - a(2)) * b(2)) / (b(1)**2+b(2)**2)
    d(1) = p_old(1) - a(1) - t*b(1)
    d(2) = p_old(2) - a(2) - t*b(2)

    p_new(1) = p_old(1) - 2*d(1)
    p_new(2) = p_old(2) - 2*d(2)

  end subroutine Mirror_P_atLine


end module DGBasis_mod

