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

module f_ctof_bind_interoperability
  implicit none

  interface f_interoperability_faultOutput
    module procedure f_interoperability_faultOutput
  end interface

  interface f_interoperability_evaluateFrictionLaw
    module procedure f_interoperability_evaluateFrictionLaw
  end interface

  contains
    subroutine copyDynamicRuptureState(domain, fromMeshId, toMeshId)
      use typesDef
      implicit none
      type(tUnstructDomainDescript), pointer :: domain
      integer     :: fromMeshId, toMeshId,iFace
      

      domain%disc%DynRup%output_Mu(:,fromMeshId:toMeshId)             = domain%disc%DynRup%Mu(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Strength(:,fromMeshId:toMeshId)       = domain%disc%DynRup%Strength(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Slip(:,fromMeshId:toMeshId)           = domain%disc%DynRup%Slip(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Slip1(:,fromMeshId:toMeshId)          = domain%disc%DynRup%Slip1(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_Slip2(:,fromMeshId:toMeshId)          = domain%disc%DynRup%Slip2(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_rupture_time(:,fromMeshId:toMeshId)   = domain%disc%DynRup%rupture_time(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_PeakSR(:,fromMeshId:toMeshId)         = domain%disc%DynRup%PeakSR(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_dynStress_time(:,fromMeshId:toMeshId) = domain%disc%DynRup%dynStress_time(:,fromMeshId:toMeshId)
      domain%disc%DynRup%output_StateVar(:,fromMeshId:toMeshId)       = domain%disc%DynRup%StateVar(:,fromMeshId:toMeshId)
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
      real(kind=8), pointer                  :: l_time

      type(c_ptr), value                     :: i_timeStepWidth
      real(kind=8), pointer                  :: l_timeStepWidth

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

    subroutine f_interoperability_evaluateFrictionLaw( i_domain, i_face, i_QInterpolatedPlus, i_QInterpolatedMinus, &
      i_imposedStatePlus, i_imposedStateMinus, i_numberOfPoints, i_godunovLd, i_time, timePoints, timeWeights, densityPlus, &
      pWaveVelocityPlus, sWaveVelocityPlus, densityMinus, pWaveVelocityMinus, sWaveVelocityMinus, c_resampleMatrix ) bind (c, name='f_interoperability_evaluateFrictionLaw')
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

      type(c_ptr), value                     :: i_QInterpolatedPlus
      REAL_TYPE, pointer                     :: l_QInterpolatedPlus(:,:,:)

      type(c_ptr), value                     :: i_QInterpolatedMinus
      REAL_TYPE, pointer                     :: l_QInterpolatedMinus(:,:,:)

      type(c_ptr), value                     :: i_imposedStatePlus
      REAL_TYPE, pointer                     :: l_imposedStatePlus(:,:)

      type(c_ptr), value                     :: i_imposedStateMinus
      REAL_TYPE, pointer                     :: l_imposedStateMinus(:,:)

      type(c_ptr), intent(in), value         :: c_resampleMatrix
      REAL_TYPE, pointer                     :: resampleMatrix(:,:)

      type(c_ptr), value                     :: i_time
      real(kind=8), pointer                  :: l_time

      real(c_double), intent(in), dimension(CONVERGENCE_ORDER)  :: timePoints
      real(c_double), intent(in), dimension(CONVERGENCE_ORDER)  :: timeWeights

      real(kind=c_double), value             :: densityPlus, pWaveVelocityPlus, sWaveVelocityPlus, densityMinus, pWaveVelocityMinus, sWaveVelocityMinus

      REAL        :: rho, rho_neig, Zp_inv, Zp_neig_inv, Zs_inv, Zs_neig_inv, eta_p, eta_s
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
      call c_f_pointer( i_QInterpolatedPlus,  l_QInterpolatedPlus, [i_godunovLd,9,CONVERGENCE_ORDER])
      call c_f_pointer( i_QInterpolatedMinus, l_QInterpolatedMinus, [i_godunovLd,9,CONVERGENCE_ORDER])
      call c_f_pointer( i_imposedStatePlus,   l_imposedStatePlus, [i_godunovLd,9])
      call c_f_pointer( i_imposedStateMinus,  l_imposedStateMinus, [i_godunovLd,9])
      call c_f_pointer( c_resampleMatrix,     resampleMatrix, [i_numberOfPoints, i_numberOfPoints])
      call c_f_pointer( i_time,               l_time  )
      
      call copyDynamicRuptureState(l_domain, i_face, i_face)

      iElem               = l_domain%MESH%Fault%Face(i_face,1,1)          ! Remark:
      iSide               = l_domain%MESH%Fault%Face(i_face,2,1)          ! iElem denotes "+" side

      rho = densityPlus
      w_speed(:) = (/ pWaveVelocityPlus, sWaveVelocityPlus, sWaveVelocityPlus /)
      rho_neig = densityMinus
      w_speed_neig(:) = (/ pWaveVelocityMinus, sWaveVelocityMinus, sWaveVelocityMinus /)

      Zp_inv = 1d0 / (rho * pWaveVelocityPlus)
      Zp_neig_inv = 1d0 / (rho_neig * pWaveVelocityMinus)
      Zs_inv = 1d0 / (rho * sWaveVelocityPlus)
      Zs_neig_inv = 1d0 / (rho_neig * sWaveVelocityMinus)

      eta_p = 1d0 / (Zp_inv + Zp_neig_inv)
      eta_s = 1d0 / (Zs_inv + Zs_neig_inv)

      do j=1,CONVERGENCE_ORDER
        do i=1,i_numberOfPoints
        NorStressGP(i,j) = eta_p * (l_QInterpolatedMinus(i,7,j) - l_QInterpolatedPlus(i,7,j) +&
                                    l_QInterpolatedPlus(i,1,j) * Zp_inv + l_QInterpolatedMinus(i,1,j) * Zp_neig_inv)
        XYStressGP(i,j)  = eta_s * (l_QInterpolatedMinus(i,8,j) - l_QInterpolatedPlus(i,8,j) +&
                                    l_QInterpolatedPlus(i,4,j) * Zs_inv + l_QInterpolatedMinus(i,4,j) * Zs_neig_inv)
        XZStressGP(i,j)  = eta_s * (l_QInterpolatedMinus(i,9,j) - l_QInterpolatedPlus(i,9,j) +&
                                    l_QInterpolatedPlus(i,6,j) * Zs_inv + l_QInterpolatedMinus(i,6,j) * Zs_neig_inv)
        enddo
      enddo

      call Eval_friction_law( TractionGP_XY,TractionGP_XZ,        & ! OUT: updated Traction
                              NorStressGP,XYStressGP,XZStressGP,  & ! IN: Godunov status
                              i_face,iSide,iElem,l_time,timePoints,          & ! IN: element ID, time, inv Trafo
                              rho,rho_neig,w_speed,w_speed_neig, resampleMatrix,  & ! IN: background values
                              l_domain%eqn, l_domain%disc, l_domain%mesh, l_domain%mpi, l_domain%io, l_domain%bnd)

      l_imposedStatePlus = 0.0
      l_imposedStateMinus = 0.0

      do j=1,CONVERGENCE_ORDER
        do i=1,i_numberOfPoints
          l_imposedStateMinus(i,1) = l_imposedStateMinus(i,1) + timeWeights(j) * NorStressGP(i,j)
          l_imposedStateMinus(i,4) = l_imposedStateMinus(i,4) + timeWeights(j) * TractionGP_XY(i,j)
          l_imposedStateMinus(i,6) = l_imposedStateMinus(i,6) + timeWeights(j) * TractionGP_XZ(i,j)
          l_imposedStateMinus(i,7) = l_imposedStateMinus(i,7) + timeWeights(j) * (l_QInterpolatedMinus(i,7,j) - Zp_neig_inv * (NorStressGP(i,j)-l_QInterpolatedMinus(i,1,j)))
          l_imposedStateMinus(i,8) = l_imposedStateMinus(i,8) + timeWeights(j) * (l_QInterpolatedMinus(i,8,j) - Zs_neig_inv * (TractionGP_XY(i,j)-l_QInterpolatedMinus(i,4,j)))
          l_imposedStateMinus(i,9) = l_imposedStateMinus(i,9) + timeWeights(j) * (l_QInterpolatedMinus(i,9,j) - Zs_neig_inv * (TractionGP_XZ(i,j)-l_QInterpolatedMinus(i,6,j)))

          l_imposedStatePlus(i,1) = l_imposedStatePlus(i,1) + timeWeights(j) * NorStressGP(i,j)
          l_imposedStatePlus(i,4) = l_imposedStatePlus(i,4) + timeWeights(j) * TractionGP_XY(i,j)
          l_imposedStatePlus(i,6) = l_imposedStatePlus(i,6) + timeWeights(j) * TractionGP_XZ(i,j)
          l_imposedStatePlus(i,7) = l_imposedStatePlus(i,7) + timeWeights(j) * (l_QInterpolatedPlus(i,7,j) + Zp_inv * (NorStressGP(i,j)-l_QInterpolatedPlus(i,1,j)))
          l_imposedStatePlus(i,8) = l_imposedStatePlus(i,8) + timeWeights(j) * (l_QInterpolatedPlus(i,8,j) + Zs_inv * (TractionGP_XY(i,j)-l_QInterpolatedPlus(i,4,j)))
          l_imposedStatePlus(i,9) = l_imposedStatePlus(i,9) + timeWeights(j) * (l_QInterpolatedPlus(i,9,j) + Zs_inv * (TractionGP_XZ(i,j)-l_QInterpolatedPlus(i,6,j)))
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

    subroutine f_interoperability_fitAttenuation( domain, rho, mu, lambda, Qp, Qs, materialFitted) bind( c, name='f_interoperability_fitAttenuation')
      use iso_c_binding
      use TypesDef
      use ini_MODEL_mod

      type(c_ptr), value                     :: domain
      type(tUnstructDomainDescript), pointer :: l_domain

      real(kind=c_double), value             :: rho, mu, lambda, Qp, Qs

      type(c_ptr), value                     :: materialFitted
      real(kind=8), pointer                  :: l_materialFitted(:)

      real                                   :: material(5)

      call c_f_pointer( domain,                   l_domain                                         )
      call c_f_pointer( materialFitted,           l_materialFitted,  [l_domain%EQN%nBackgroundVar] )

      material(:) = (/ rho, mu, lambda, Qp, Qs /)
      call fitAttenuation(material, l_materialFitted, l_domain%EQN)
    end subroutine
end module
