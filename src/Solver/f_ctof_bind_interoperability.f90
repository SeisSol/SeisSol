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


    !!Code added by ADRIAN
    subroutine f_interoperability_getDynRup(i_domain, iFace, i_InitialStressInFaultCS, i_mu, i_slipRate1, i_slipRate2,i_RF) bind (c, name='f_interoperability_getDynRup')
        use iso_c_binding
        use typesDef
        use f_ftoc_bind_interoperability
        implicit none

        integer                                :: nBndGP
        type(c_ptr), value                     :: i_domain
        type(tUnstructDomainDescript), pointer :: l_domain
        integer(kind=c_int), value             :: iFace
        type(c_ptr), value                     :: i_InitialStressInFaultCS
        REAL_TYPE, pointer                     :: l_InitialStressInFaultCS(:,:)
        type(c_ptr), value                     :: i_mu
        REAL_TYPE, pointer                     :: l_mu(:)
        type(c_ptr), value                     :: i_slipRate1
        REAL_TYPE, pointer                     :: l_slipRate1(:)
        type(c_ptr), value                     :: i_slipRate2
        REAL_TYPE, pointer                     :: l_slipRate2(:)
        type(c_ptr), value                     :: i_RF
        logical(kind=C_bool), pointer          :: l_RF(:)

        call c_f_pointer( i_domain,             l_domain)
        nBndGP = l_domain%DISC%Galerkin%nBndGP

        call c_f_pointer( i_InitialStressInFaultCS, l_InitialStressInFaultCS, [nBndGP,6])
        call c_f_pointer( i_mu, l_mu, [nBndGP])
        call c_f_pointer( i_slipRate1, l_slipRate1, [nBndGP])
        call c_f_pointer( i_slipRate2, l_slipRate2, [nBndGP])
        call c_f_pointer( i_RF, l_RF, [nBndGP])

        !DISC%DynRup%StateVar(:,:)      = EQN%IniStateVar
        l_InitialStressInFaultCS(:,:)   = l_domain%EQN%InitialStressInFaultCS(:,:,iFace)
        l_mu(:)                         = l_domain%EQN%IniMu(:,iFace)
        l_slipRate1(:)                  = l_domain%EQN%IniSlipRate1
        l_slipRate2(:)                  = l_domain%EQN%IniSlipRate2
        l_RF(:)                         = l_domain%DISC%DynRup%RF(:,iFace)
    end subroutine

    subroutine f_interoperability_getDynRupStateVar(i_domain, iFace, i_stateVar) bind (c, name='f_interoperability_getDynRupStateVar')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      integer                                :: nBndGP, i
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer(kind=c_int), value             :: iFace
      type(c_ptr), value                     :: i_stateVar
      REAL_TYPE, pointer                     :: l_stateVar(:)

      call c_f_pointer( i_domain,             l_domain)
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_stateVar, l_stateVar, [nBndGP])

      l_stateVar     = l_domain%EQN%IniStateVar(:,iFace)
    end subroutine

    subroutine f_interoperability_getDynRupNucStress(i_domain, iFace ,i_nucleationStressInFaultCS) bind (c, name='f_interoperability_getDynRupNucStress')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      integer                                :: nBndGP
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer(kind=c_int), value             :: iFace
      type(c_ptr), value                     :: i_nucleationStressInFaultCS
      REAL_TYPE, pointer                     :: l_nucleationStressInFaultCS(:,:)

      call c_f_pointer( i_domain,             l_domain)
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_nucleationStressInFaultCS, l_nucleationStressInFaultCS, [nBndGP,6])

      l_nucleationStressInFaultCS(:,:)   = l_domain%EQN%NucleationStressInFaultCS(:,:,iFace)
    end subroutine

    subroutine f_interoperability_getDynRupFL_3(i_domain, iFace ,i_RS_a, i_RS_sl0, i_RS_sr0) bind (c, name='f_interoperability_getDynRupFL_3')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      integer                                :: nBndGP
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer(kind=c_int), value             :: iFace
      type(c_ptr), value                     :: i_RS_a
      real*8, pointer                        :: l_RS_a
      type(c_ptr), value                     :: i_RS_sl0
      real*8, pointer                        :: l_RS_sl0
      type(c_ptr), value                     :: i_RS_sr0
      real*8, pointer                        :: l_RS_sr0

      call c_f_pointer( i_domain,             l_domain)
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_RS_a,      l_RS_a)
      call c_f_pointer( i_RS_sl0,    l_RS_sl0)
      call c_f_pointer( i_RS_sr0,    l_RS_sr0)

      l_RS_a                    = l_domain%DISC%DynRup%RS_a
      l_RS_sl0                  = l_domain%DISC%DynRup%RS_sl0
      l_RS_sr0                  = l_domain%DISC%DynRup%RS_sr0

    end subroutine


    subroutine f_interoperability_getDynRupTP(i_domain, i_TP_grid, i_TP_DFinv) bind (c, name='f_interoperability_getDynRupTP')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      integer                                :: TP_grid_nz
      integer                                :: i_numberOfPoints, iFace
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      type(c_ptr), value                     :: i_TP_grid
      REAL_TYPE, pointer                     :: l_TP_grid(:)
      type(c_ptr), value                     :: i_TP_DFinv
      REAL_TYPE, pointer                     :: l_TP_DFinv(:)

      call c_f_pointer( i_domain,             l_domain)
      TP_grid_nz = l_domain%DISC%DynRup%TP_grid_nz
      i_numberOfPoints = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_TP_grid,       l_TP_grid , [TP_grid_nz])
      call c_f_pointer( i_TP_DFinv,      l_TP_DFinv, [TP_grid_nz])

      l_TP_grid(:)   = l_domain%DISC%DynRup%TP_grid(:)
      l_TP_DFinv(:)  = l_domain%DISC%DynRup%TP_DFinv(:)
    end subroutine

    subroutine f_interoperability_setFrictionOutputGeneral(i_domain, i_face, &
              i_slip, i_slipStrike, i_slipDip, i_ruptureTime, i_dynStressTime,&
              i_PeakSlipRate, i_tractionXY, i_tractionXZ)&
        bind (c, name='f_interoperability_setFrictionOutputGeneral')
        use iso_c_binding
        use typesDef
        use f_ftoc_bind_interoperability
        implicit none

        INTEGER     :: i ,j, k
        type(c_ptr), value                     :: i_domain
        type(tUnstructDomainDescript), pointer :: l_domain
        integer                                :: nSide , nBndGP, iBndGP
        integer(kind=c_int), value             :: i_face

        type(c_ptr), value                     :: i_slip
        REAL_TYPE, pointer                     :: l_slip(:)
        type(c_ptr), value                     :: i_slipStrike
        REAL_TYPE, pointer                     :: l_slipStrike(:)
        type(c_ptr), value                     :: i_slipDip
        REAL_TYPE, pointer                     :: l_slipDip(:)

        type(c_ptr), value                     :: i_ruptureTime
        REAL_TYPE, pointer                     :: l_ruptureTime(:)
        type(c_ptr), value                     :: i_dynStressTime
        REAL_TYPE, pointer                     :: l_dynStressTime(:)
        type(c_ptr), value                     :: i_PeakSlipRate
        REAL_TYPE, pointer                     :: l_PeakSlipRate(:)
        type(c_ptr), value                     :: i_tractionXY
        REAL_TYPE, pointer                     :: l_tractionXY(:)
        type(c_ptr), value                     :: i_tractionXZ
        REAL_TYPE, pointer                     :: l_tractionXZ(:)

        ! convert c to fortran pointers
        call c_f_pointer( i_domain,             l_domain)
        nSide = l_domain%MESH%Fault%nSide
        nBndGP = l_domain%DISC%Galerkin%nBndGP

        call c_f_pointer( i_slip, l_slip, [nBndGP])
        call c_f_pointer( i_slipStrike, l_slipStrike, [nBndGP])
        call c_f_pointer( i_slipDip, l_slipDip, [nBndGP])
        call c_f_pointer( i_ruptureTime, l_ruptureTime, [nBndGP])
        call c_f_pointer( i_dynStressTime, l_dynStressTime, [nBndGP])
        call c_f_pointer( i_PeakslipRate, l_PeakslipRate, [nBndGP])
        call c_f_pointer( i_tractionXY, l_tractionXY, [nBndGP])
        call c_f_pointer( i_tractionXZ, l_tractionXZ, [nBndGP])

        !copy to output
        l_domain%DISC%DynRup%output_Slip(:,i_face)                  = l_slip(:)
        l_domain%DISC%DynRup%output_Slip1(:,i_face)                 = l_slipStrike(:)
        l_domain%DISC%DynRup%output_Slip2(:,i_face)                 = l_slipDip(:)
        l_domain%DISC%DynRup%output_rupture_time(:,i_face)          = l_ruptureTime(:)
        l_domain%DISC%DynRup%output_dynStress_time(:,i_face)        = l_dynStressTime(:)
        l_domain%DISC%DynRup%rupture_time(:,i_face)                 = l_ruptureTime(:)
        l_domain%DISC%DynRup%output_PeakSR(:,i_face)                = l_PeakSlipRate(:)
        l_domain%DISC%DynRup%TracXY(:,i_face)                       = l_tractionXY(:)
        l_domain%DISC%DynRup%TracXZ(:,i_face)                       = l_tractionXZ(:)
    end subroutine

    subroutine f_interoperability_setFrictionOutputSpecific(i_domain, i_face, &
          i_averagedSlip, i_slipRate1, i_slipRate2, i_mu)&
          bind (c, name='f_interoperability_setFrictionOutputSpecific')

      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      INTEGER     :: i ,j, k
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer                                :: nSide , nBndGP
      integer(kind=c_int), value             :: i_face
      type(c_ptr), value                     :: i_averagedSlip
      REAL_TYPE, pointer                     :: l_averagedSlip
      type(c_ptr), value                     :: i_slipRate1
      REAL_TYPE, pointer                     :: l_slipRate1(:)
      type(c_ptr), value                     :: i_slipRate2
      REAL_TYPE, pointer                     :: l_slipRate2(:)
      type(c_ptr), value                     :: i_mu
      REAL_TYPE, pointer                     :: l_mu(:)

      ! convert c to fortran pointers
      call c_f_pointer( i_domain,             l_domain)
      nSide = l_domain%MESH%Fault%nSide
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_averagedSlip, l_averagedSlip)
      call c_f_pointer( i_slipRate2, l_slipRate2, [nBndGP])
      call c_f_pointer( i_slipRate1, l_slipRate1, [nBndGP])
      call c_f_pointer( i_mu, l_mu, [nBndGP])

      !copy to output
      ! average slip only copied back if allocated
      IF (l_domain%DISC%DynRup%magnitude_output_on.EQ.1) THEN
        l_domain%DISC%DynRup%averaged_Slip(i_face)        = l_averagedSlip
      ENDIF
      l_domain%DISC%DynRup%SlipRate1(:,i_face)                    = l_slipRate1(:)
      l_domain%DISC%DynRup%SlipRate2(:,i_face)                    = l_slipRate2(:)
      l_domain%DISC%DynRup%output_Mu(:,i_face)                    = l_mu(:)
    end subroutine

    !!Code added by ADRIAN
    subroutine f_interoperability_setFrictionOutputStateVar(i_domain, i_face, i_stateVar)&
            bind (c, name='f_interoperability_setFrictionOutputStateVar')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      INTEGER     :: i ,j, k
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer                                :: nSide , nBndGP
      integer(kind=c_int), value             :: i_face
      type(c_ptr), value                     :: i_stateVar
      REAL_TYPE, pointer                     :: l_stateVar(:)

      !integer :: nSide
      ! convert c to fortran pointers
      call c_f_pointer( i_domain,             l_domain)
      nSide = l_domain%MESH%Fault%nSide
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_stateVar, l_stateVar, [nBndGP])
      !copy to output

      l_domain%DISC%DynRup%output_StateVar(:,i_face) = l_stateVar(:)  !l_domain%DISC%DynRup%dynStress_time(:,i_face)
    end subroutine



    !!Code added by ADRIAN
    subroutine f_interoperability_setFrictionOutputStrength(i_domain, i_face, i_strength)&
            bind (c, name='f_interoperability_setFrictionOutputStrength')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      INTEGER     :: i ,j, k
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer                                :: nSide , nBndGP
      integer(kind=c_int), value             :: i_face
      type(c_ptr), value                     :: i_strength
      REAL_TYPE, pointer                     :: l_strength(:)

      !integer :: nSide
      ! convert c to fortran pointers
      call c_f_pointer( i_domain,             l_domain)
      nSide = l_domain%MESH%Fault%nSide
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_strength, l_strength, [nBndGP])
      !copy to output

      l_domain%DISC%DynRup%output_Strength(:,i_face) = l_strength(:)  !l_domain%DISC%DynRup%dynStress_time(:,i_face)
    end subroutine

    subroutine f_interoperability_setFrictionOutputThermalPressurization(i_domain, i_face, i_fluidPressure, i_fluidTemperature)&
        bind (c, name='f_interoperability_setFrictionOutputThermalPressurization')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      INTEGER     :: i ,j, k
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer                                :: nSide , nBndGP
      integer(kind=c_int), value             :: i_face
      type(c_ptr), value                     :: i_fluidPressure
      type(c_ptr), value                     :: i_fluidTemperature
      REAL_TYPE, pointer                     :: l_fluidPressure(:)
      REAL_TYPE, pointer                     :: l_fluidTemperature(:)

      !integer :: nSide
      ! convert c to fortran pointers
      call c_f_pointer( i_domain,             l_domain)
      nSide = l_domain%MESH%Fault%nSide
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_fluidPressure, l_fluidPressure, [nBndGP])
      call c_f_pointer( i_fluidTemperature, l_fluidTemperature, [nBndGP])
      !copy to output

      l_domain%DISC%DynRup%TP(:,i_face,2) = l_fluidPressure(:)
      l_domain%DISC%DynRup%TP(:,i_face,1) = l_fluidTemperature(:)
    end subroutine

    !!Code added by ADRIAN
    subroutine f_interoperability_setFrictionOutputInitialStress(i_domain, iFace, i_InitialStressInFaultCS, i_bulkXX, i_bulkYY, i_bulkZZ, i_shearXY, i_shearYZ, i_shearXZ) bind (c, name='f_interoperability_setFrictionOutputInitialStress')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      integer                                :: i, nBndGP
      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain
      integer(kind=c_int), value             :: iFace
      type(c_ptr), value                     :: i_InitialStressInFaultCS
      REAL_TYPE, pointer                     :: l_InitialStressInFaultCS(:,:)
      type(c_ptr), value                     :: i_bulkXX
      REAL_TYPE, pointer                     :: l_bulkXX(:)
      type(c_ptr), value                     :: i_bulkYY
      REAL_TYPE, pointer                     :: l_bulkYY(:)
      type(c_ptr), value                     :: i_bulkZZ
      REAL_TYPE, pointer                     :: l_bulkZZ(:)
      type(c_ptr), value                     :: i_shearXY
      REAL_TYPE, pointer                     :: l_shearXY(:)
      type(c_ptr), value                     :: i_shearYZ
      REAL_TYPE, pointer                     :: l_shearYZ(:)
      type(c_ptr), value                     :: i_shearXZ
      REAL_TYPE, pointer                     :: l_shearXZ(:)

      call c_f_pointer( i_domain,             l_domain)
      nBndGP = l_domain%DISC%Galerkin%nBndGP

      call c_f_pointer( i_InitialStressInFaultCS, l_InitialStressInFaultCS, [6,nBndGP])
      call c_f_pointer( i_bulkXX, l_bulkXX, [nBndGP])
      call c_f_pointer( i_bulkYY, l_bulkYY, [nBndGP])
      call c_f_pointer( i_bulkZZ, l_bulkZZ, [nBndGP])
      call c_f_pointer( i_shearXY, l_shearXY, [nBndGP])
      call c_f_pointer( i_shearYZ, l_shearYZ, [nBndGP])
      call c_f_pointer( i_shearXZ, l_shearXZ, [nBndGP])
      DO i = 1, 6
        l_domain%EQN%InitialStressInFaultCS(:,i,iFace) = l_InitialStressInFaultCS(i,1:nBndGP)
      END DO
      l_domain%eqn%inibulk_xx(:,iFace) = l_bulkXX(:)
      l_domain%eqn%inibulk_yy(:,iFace) = l_bulkYY(:)
      l_domain%eqn%inibulk_zz(:,iFace) = l_bulkZZ(:)
      l_domain%eqn%inishearxy(:,iFace) = l_shearXY(:)
      l_domain%eqn%inishearyz(:,iFace) = l_shearYZ(:)
      l_domain%eqn%inishearxz(:,iFace) = l_shearXZ(:)
    end subroutine

    subroutine f_interoperability_initializeFaultOutput(i_domain) bind (c, name="f_interoperability_initializeFaultOutput")
      use iso_c_binding
      use typesDef
      use ini_faultoutput_mod
      implicit none
      type( c_ptr ), value         :: i_domain
      type(tUnstructDomainDescript), pointer :: domain

      call c_f_pointer(i_domain, domain)

      call ini_fault_subsampled(domain%EQN,domain%MESH,domain%BND,domain%DISC,domain%IO,domain%MPI)

    end subroutine f_interoperability_initializeFaultOutput


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
