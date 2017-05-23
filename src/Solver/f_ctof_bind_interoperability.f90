!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!! @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
!!
!! @section LICENSE
!! Copyright (c) 2015-2017, SeisSol Group
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
!! C++/Fortran-interoperability: Fortran-bindings

#include "Initializer/preProcessorMacros.fpp"
#include "initialization/precision.h"

#if defined(DOUBLE_PRECISION)
#define REAL_TYPE real*8
#elif defined(SINGLE_PRECISION)
#define REAL_TYPE real*4
#else
#error Unknown floating point precision type.
#endif

module f_ctof_bind_interoperability
  implicit none

  interface f_interoperability_faultOutput
    module procedure f_interoperability_faultOutput
  end interface

  interface f_interoperability_evaluateFrictionLaw
    module procedure f_interoperability_evaluateFrictionLaw
  end interface

  interface f_interoperability_writeReceivers
    module procedure f_interoperability_writeReceivers
  end interface

  contains
    subroutine copyDynamicRuptureState(domain, fromMeshId, toMeshId)
      use typesDef
      implicit none
      type(tUnstructDomainDescript), pointer :: domain
      integer     :: fromMeshId, toMeshId
      

      domain%disc%DynRup%output_Mu(:,fromMeshId:toMeshId)             = domain%disc%DynRup%Mu(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Strength(:,fromMeshId:toMeshId)       = domain%disc%DynRup%Strength(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Slip(:,fromMeshId:toMeshId)           = domain%disc%DynRup%Slip(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Slip1(:,fromMeshId:toMeshId)          = domain%disc%DynRup%Slip1(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Slip2(:,fromMeshId:toMeshId)          = domain%disc%DynRup%Slip2(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_rupture_time(:,fromMeshId:toMeshId)   = domain%disc%DynRup%rupture_time(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_PeakSR(:,fromMeshId:toMeshId)         = domain%disc%DynRup%PeakSR(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_dynStress_time(:,fromMeshId:toMeshId) = domain%disc%DynRup%dynStress_time(:,fromMeshId:toMeshId)
        
      domain%disc%DynRup%output_StateVar(fromMeshId:toMeshId,:)       = domain%disc%DynRup%StateVar(fromMeshId:toMeshId,:)
      
    end subroutine
    !
    ! C to fortran bindings
    !
    subroutine f_interoperability_copyDynamicRuptureState( i_domain ) bind (c, name='f_interoperability_copyDynamicRuptureState')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain

      ! convert c to fortran pointers
      call c_f_pointer( i_domain, l_domain)
      
      call copyDynamicRuptureState(l_domain, 1, l_domain%mesh%Fault%nSide)
    end subroutine

    subroutine f_interoperability_faultOutput( i_domain, i_time, i_timeStepWidth ) bind (c, name='f_interoperability_faultOutput')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      use faultoutput_mod
      implicit none

      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain

      type(c_ptr), value                     :: i_time
      real*8, pointer                        :: l_time

      type(c_ptr), value                     :: i_timeStepWidth
      real*8, pointer                        :: l_timeStepWidth

      ! register scorep region dynamic rupture output (receiver)
      SCOREP_USER_REGION_DEFINE( r_dr_output )

      ! convert c to fortran pointers
      call c_f_pointer( i_domain,        l_domain)
      call c_f_pointer( i_time,          l_time  )
      call c_f_pointer( i_timeStepWidth, l_timeStepWidth )

      SCOREP_USER_REGION_BEGIN( r_dr_output, "fault_output_receiver", SCOREP_USER_REGION_TYPE_COMMON )
      ! NOTE: This will only handle faul receivers
      call faultoutput(l_domain%eqn, l_domain%disc, l_domain%mesh, l_domain%io, l_domain%mpi, l_domain%optionalFields%BackgroundValue, l_domain%bnd, l_time, l_timeStepWidth)
      SCOREP_USER_REGION_END( r_dr_output )

      l_domain%disc%iterationstep = l_domain%disc%iterationstep + 1
    end subroutine

    subroutine f_interoperability_evaluateFrictionLaw( i_domain, i_face, i_godunov, i_imposedStatePlus, i_imposedStateMinus, i_numberOfPoints, i_godunovLd, i_time, timePoints, timeWeights, densityPlus, pWaveVelocityPlus, sWaveVelocityPlus, densityMinus, pWaveVelocityMinus, sWaveVelocityMinus ) bind (c, name='f_interoperability_evaluateFrictionLaw')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      use Eval_friction_law_mod
      use JacobiNormal_mod, only: RotationMatrix3D
      implicit none

      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain

      integer(kind=c_int), value             :: i_face
      integer(kind=c_int), value             :: i_numberOfPoints
      integer(kind=c_int), value             :: i_godunovLd

      type(c_ptr), value                     :: i_godunov
      REAL_TYPE, pointer                     :: l_godunov(:,:,:)

      type(c_ptr), value                     :: i_imposedStatePlus
      REAL_TYPE, pointer                     :: l_imposedStatePlus(:,:)

      type(c_ptr), value                     :: i_imposedStateMinus
      REAL_TYPE, pointer                     :: l_imposedStateMinus(:,:)

      type(c_ptr), value                     :: i_time
      real*8, pointer                        :: l_time

      real(c_double), intent(in), dimension(CONVERGENCE_ORDER)  :: timePoints
      real(c_double), intent(in), dimension(CONVERGENCE_ORDER)  :: timeWeights

      real(kind=c_double), value             :: densityPlus, pWaveVelocityPlus, sWaveVelocityPlus, densityMinus, pWaveVelocityMinus, sWaveVelocityMinus

      REAL        :: rho, rho_neig
      REAL        :: w_speed(3),w_speed_neig(3)

      REAL        :: TractionGP_XY(1:i_numberOfPoints,CONVERGENCE_ORDER)
      REAL        :: TractionGP_XZ(1:i_numberOfPoints,CONVERGENCE_ORDER)
      REAL        :: NorStressGP(1:i_numberOfPoints,CONVERGENCE_ORDER)
      REAL        :: XYStressGP(1:i_numberOfPoints,CONVERGENCE_ORDER)
      REAL        :: XZStressGP(1:i_numberOfPoints,CONVERGENCE_ORDER)
      real        :: subTimeStepWidth
      integer :: iSide, iElem, iObject, MPIIndex, MPIIndex_DR, i, j

      ! register scorep region dynamic rupture
      SCOREP_USER_REGION_DEFINE( r_dr )
      SCOREP_USER_REGION_BEGIN( r_dr, "friction_law", SCOREP_USER_REGION_TYPE_COMMON )

      ! convert c to fortran pointers
      call c_f_pointer( i_domain,             l_domain)
      call c_f_pointer( i_godunov,            l_godunov, [i_godunovLd,9,CONVERGENCE_ORDER])
      call c_f_pointer( i_imposedStatePlus,   l_imposedStatePlus, [i_godunovLd,9])
      call c_f_pointer( i_imposedStateMinus,  l_imposedStateMinus, [i_godunovLd,9])
      call c_f_pointer( i_time,               l_time  )
      
      call copyDynamicRuptureState(l_domain, i_face, i_face)

      iElem               = l_domain%MESH%Fault%Face(i_face,1,1)          ! Remark:
      iSide               = l_domain%MESH%Fault%Face(i_face,2,1)          ! iElem denotes "+" side

      rho = densityPlus
      w_speed(:) = (/ pWaveVelocityPlus, sWaveVelocityPlus, sWaveVelocityPlus /)
      rho_neig = densityMinus
      w_speed_neig(:) = (/ pWaveVelocityMinus, sWaveVelocityMinus, sWaveVelocityMinus /)

      do j=1,CONVERGENCE_ORDER
        do i=1,i_numberOfPoints
          NorStressGP(i,j) = l_godunov(i,1,j)
          XYStressGP(i,j) = l_godunov(i,4,j)
          XZStressGP(i,j) = l_godunov(i,6,j)
        enddo
      enddo

      call Eval_friction_law( TractionGP_XY,TractionGP_XZ,        & ! OUT: updated Traction
                              NorStressGP,XYStressGP,XZStressGP,  & ! IN: Godunov status
                              i_face,iSide,iElem,l_time,timePoints,          & ! IN: element ID, time, inv Trafo
                              rho,rho_neig,w_speed,w_speed_neig,  & ! IN: background values
                              l_domain%eqn, l_domain%disc, l_domain%mesh, l_domain%mpi, l_domain%io, l_domain%bnd)

      l_imposedStatePlus = 0.0
      l_imposedStateMinus = 0.0

      do j=1,CONVERGENCE_ORDER
        do i=1,i_numberOfPoints
          l_imposedStateMinus(i,1) = l_imposedStateMinus(i,1) + timeWeights(j) * l_godunov(i,1,j)
          l_imposedStateMinus(i,4) = l_imposedStateMinus(i,4) + timeWeights(j) * TractionGP_XY(i,j)
          l_imposedStateMinus(i,6) = l_imposedStateMinus(i,6) + timeWeights(j) * TractionGP_XZ(i,j)
          l_imposedStateMinus(i,7) = l_imposedStateMinus(i,7) + timeWeights(j) * l_godunov(i,7,j)
          l_imposedStateMinus(i,8) = l_imposedStateMinus(i,8) + timeWeights(j) * (l_godunov(i,8,j) - 1.0D0/(w_speed_neig(2)*rho_neig) * (TractionGP_XY(i,j)-l_godunov(i,4,j)))
          l_imposedStateMinus(i,9) = l_imposedStateMinus(i,9) + timeWeights(j) * (l_godunov(i,9,j) - 1.0D0/(w_speed_neig(2)*rho_neig) * (TractionGP_XZ(i,j)-l_godunov(i,6,j)))

          l_imposedStatePlus(i,1) = l_imposedStatePlus(i,1) + timeWeights(j) * l_godunov(i,1,j)
          l_imposedStatePlus(i,4) = l_imposedStatePlus(i,4) + timeWeights(j) * TractionGP_XY(i,j)
          l_imposedStatePlus(i,6) = l_imposedStatePlus(i,6) + timeWeights(j) * TractionGP_XZ(i,j)
          l_imposedStatePlus(i,7) = l_imposedStatePlus(i,7) + timeWeights(j) * l_godunov(i,7,j)
          l_imposedStatePlus(i,8) = l_imposedStatePlus(i,8) + timeWeights(j) * (l_godunov(i,8,j) + 1.0D0/(w_speed(2)*rho) * (TractionGP_XY(i,j)-l_godunov(i,4,j)))
          l_imposedStatePlus(i,9) = l_imposedStatePlus(i,9) + timeWeights(j) * (l_godunov(i,9,j) + 1.0D0/(w_speed(2)*rho) * (TractionGP_XZ(i,j)-l_godunov(i,6,j)))
        enddo
      enddo

      SCOREP_USER_REGION_END( r_dr )
    end subroutine

    subroutine f_interoperability_calcElementwiseFaultoutput(i_domain, time) bind (c, name="f_interoperability_calcElementwiseFaultoutput")
      use iso_c_binding
      use typesDef
      use faultoutput_mod
      implicit none

      type( c_ptr ), value         :: i_domain
      real( kind=c_double ), value :: time

      type(tUnstructDomainDescript), pointer :: domain
      integer :: OutputPointType

      ! convert c to fortran pointers
      call c_f_pointer(i_domain, domain)

      OutputPointType = domain%DISC%DynRup%OutputPointType
      domain%DISC%DynRup%OutputPointType = 4
      call calc_FaultOutput(domain%DISC%DynRup%DynRup_out_elementwise, domain%DISC, domain%EQN, domain%MESH, &
          domain%optionalFields%BackgroundValue, domain%BND, time)
      domain%DISC%DynRup%OutputPointType = OutputPointType
    end subroutine f_interoperability_calcElementwiseFaultoutput

    subroutine f_interoperability_computePlasticity( i_domain, i_timeStep, &
            i_numberOfAlignedBasisFunctions, i_plasticParameters, i_initialLoading, io_dofs, io_Energy, io_pstrain ) bind( c, name='f_interoperability_computePlasticity')
      use iso_c_binding
      use typesDef
      use plasticity_mod

      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain

      type(c_ptr), value                     :: i_timeStep
      real*8, pointer                        :: l_timeStep

      integer(kind=c_int), value             :: i_numberOfAlignedBasisFunctions

      type(c_ptr), value                     :: i_plasticParameters
      real*8, pointer                        :: l_plasticParameters(:)

      type(c_ptr), value                     :: i_initialLoading
      real*8, pointer                        :: l_initialLoading(:,:)

      type(c_ptr), value                     :: io_dofs
      real*8, pointer                        :: l_dofs(:,:)

      type(c_ptr), value                     :: io_Energy
      real*8, pointer                        :: l_Energy(:)

      type(c_ptr), value                     :: io_pstrain
      real*8, pointer                        :: l_pstrain(:)

      ! convert c to fotran pointers
      call c_f_pointer( i_domain,         l_domain                                         )
      call c_f_pointer( i_timeStep,       l_timeStep                                       )
      call c_f_pointer( i_plasticParameters, l_plasticParameters, [4]                      )
      call c_f_pointer( i_initialLoading, l_initialLoading, [NUMBER_OF_BASIS_FUNCTIONS,6]  )
      call c_f_pointer( io_dofs,          l_dofs,       [i_numberOfAlignedBasisFunctions,9])
      call c_f_pointer( io_Energy,        l_Energy, [3]                             )
      call c_f_pointer( io_pstrain,       l_pstrain,    [7]                                )


      select case(l_domain%eqn%PlastMethod) !two different methods to check plasticity

      case(0) !values at internal GP = high-order points
              call plasticity_3d_high( dgvar    = l_dofs, &
                                  dofStress     = l_initialLoading, & !l_domain%eqn%inistress, & !l_initialLoading is the same as inistress for the high-order case
                                  nDegFr        = NUMBER_OF_BASIS_FUNCTIONS, &
                                  nAlignedDegFr = i_numberOfAlignedBasisFunctions, &
                                  tv            = l_domain%eqn%Tv, &
                                  dt            = l_timeStep, &
                                  mu            = l_domain%eqn%mu, &
                                  lambda        = l_domain%eqn%lambda, &
                                  parameters    = l_plasticParameters, &
                                  Energy        = l_Energy, &
                                  pstrain       = l_pstrain, &
                                  intGaussP     = l_domain%disc%Galerkin%intGaussP_Tet,&
                                  intGaussW     = l_domain%disc%Galerkin%intGaussW_Tet,&
                                  disc          = l_domain%disc, &
                                  nVar          = l_domain%eqn%nVar, &
                                  nIntGP        =l_domain%disc%galerkin%nIntGP)

      case(2) !average of an element
              call plasticity_3d_avg( disc      = l_domain%disc, &
                                  dgvar         = l_dofs, &
                                  dofStress     = l_initialLoading, &
                                  nDegFr        = NUMBER_OF_BASIS_FUNCTIONS, &
                                  nAlignedDegFr = i_numberOfAlignedBasisFunctions, &
                                  tv            = l_domain%eqn%Tv, &
                                  dt            = l_timeStep, &
                                  mu            = l_domain%eqn%mu, &
                                  lambda        = l_domain%eqn%lambda, &
                                  parameters    = l_plasticParameters ,&
                                  Energy        = l_Energy,&
                                  pstrain       = l_pstrain )
      end select


    end subroutine

    subroutine f_interoperability_writeReceivers( i_domain, i_fullUpdateTime, i_timeStepWidth, i_receiverTime, i_numberOfReceivers, i_receiverIds ) bind (c, name='f_interoperability_writeReceivers')
      use iso_c_binding
      use typesDef
      use receiver_mod
      implicit none

      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain

      type(c_ptr), value                     :: i_fullUpdateTime
      real*8, pointer                        :: l_fullUpdateTime

      type(c_ptr), value                     :: i_timeStepWidth
      real*8, pointer                        :: l_timeStepWidth

      type(c_ptr), value                     :: i_receiverTime
      real*8, pointer                        :: l_receiverTime

      type(c_ptr), value                     :: i_numberOfReceivers
      integer, pointer                       :: l_numberOfReceivers

      type(c_ptr), value                     :: i_receiverIds
      integer, pointer                       :: l_receiverIds(:)

      ! convert c to fortran pointers
      call c_f_pointer( i_domain,            l_domain                             )
      call c_f_pointer( i_fullUpdateTime,    l_fullUpdateTime                     )
      call c_f_pointer( i_timeStepWidth,     l_timeStepWidth                      )
      call c_f_pointer( i_receiverTime,      l_receiverTime                       )
      call c_f_pointer( i_numberOfReceivers, l_numberOfReceivers                  )
      call c_f_pointer( i_receiverIds,       l_receiverIds, [l_numberOfReceivers] )

      ! call SeisSol's receiver procedure
      call receiver( i_fullUpdateTime    = l_fullUpdateTime,    \
                     i_timeStepWidth     = l_timeStepWidth,     \
                     i_receiverTime      = l_receiverTime,      \
                     i_numberOfReceivers = l_numberOfReceivers, \
                     i_receiverIds       = l_receiverIds,       \
                     eqn                 = l_domain%eqn,        \
                     mesh                = l_domain%mesh,       \
                     disc                = l_domain%disc,       \
                     mpi                 = l_domain%mpi,        \
                     io                  = l_domain%io )
    end subroutine


    subroutine f_interoperability_computeMInvJInvPhisAtSources( i_domain, i_x, i_y, i_z, i_elem, o_mInvJInvPhisAtSources ) bind( c, name='f_interoperability_computeMInvJInvPhisAtSources')
      use iso_c_binding
      use TypesDef
      use DGBasis_mod

      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain

      real(kind=c_double), value             :: i_x, i_y, i_z
      integer(kind=c_int), value             :: i_elem

      type(c_ptr), value                     :: o_mInvJInvPhisAtSources
      real*8, pointer                        :: l_mInvJInvPhisAtSources(:)

      real                                   :: l_xi, l_eta, l_zeta
      integer                                :: indices(4) ! == MESH%nVertices_Tet
      integer                                :: l_elem, l_dof
      real                                   :: vx(4), vy(4), vz(4)

      call c_f_pointer( i_domain,                 l_domain                                              )
      call c_f_pointer( o_mInvJInvPhisAtSources,  l_mInvJInvPhisAtSources,  [NUMBER_OF_BASIS_FUNCTIONS] )

      ! f_elem = c_elem + 1
      l_elem = i_elem + 1
      indices = l_domain%MESH%ELEM%Vertex(1:l_domain%MESH%nVertices_Tet, l_elem)
      vx = l_domain%MESH%VRTX%xyNode(1, indices)
      vy = l_domain%MESH%VRTX%xyNode(2, indices)
      vz = l_domain%MESH%VRTX%xyNode(3, indices)

      call TrafoXYZ2XiEtaZeta(xi    = l_xi,   &
                              eta   = l_eta,  &
                              zeta  = l_zeta, &
                              xP    = i_x,    &
                              yP    = i_y,    &
                              zP    = i_z,    &
                              x     = vx,     &
                              y     = vy,     &
                              z     = vz,     &
                              vType = l_domain%MESH%GlobalVrtxType )

      do l_dof = 1, NUMBER_OF_BASIS_FUNCTIONS
        call BaseFunc3D(l_mInvJInvPhisAtSources(l_dof),             &
                        l_dof,                                      &
                        l_xi,                                       &
                        l_eta,                                      &
                        l_zeta,                                     &
                        l_domain%DISC%Galerkin%nPoly,               &
                        l_domain%DISC%Galerkin%cPoly3D_Tet,         &
                        l_domain%DISC%Galerkin%NonZeroCPoly_Tet,    &
                        l_domain%DISC%Galerkin%NonZeroCPolyIndex_Tet)

        l_mInvJInvPhisAtSources(l_dof) = l_mInvJInvPhisAtSources(l_dof) / ( 6.0d0 * l_domain%MESH%ELEM%Volume(l_elem) * l_domain%DISC%Galerkin%MassMatrix_Tet(l_dof, l_dof, l_domain%DISC%Galerkin%nPoly) )
      end do
    end subroutine
end module
