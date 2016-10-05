!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2015, SeisSol Group
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

module f_ftoc_bind_interoperability
  implicit none
  !
  ! Fortran to C bindings
  !

  interface c_interoperability_setDomain
    subroutine c_interoperability_setDomain( i_domain ) bind( C, name='c_interoperability_setDomain' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_domain
    end subroutine
  end interface

  interface c_interoperability_setTimeStepWidth
    subroutine c_interoperability_setTimeStepWidth( i_meshId, i_timeStepWidth ) bind( C, name='c_interoperability_setTimeStepWidth' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value :: i_meshId
      real(kind=c_double), value :: i_timeStepWidth
    end subroutine
  end interface

  interface c_interoperability_initializeClusteredLts
    subroutine c_interoperability_initializeClusteredLts( i_clustering ) bind( C, name='c_interoperability_initializeClusteredLts' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value :: i_clustering
    end subroutine
  end interface

  ! Don't forget to add // c_null_char to NRFFileName when using this interface
  interface
    subroutine c_interoperability_setupNRFPointSources( NRFFileName ) bind( C, name='c_interoperability_setupNRFPointSources' )
      use iso_c_binding, only: c_char
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: NRFFileName
    end subroutine
  end interface

  ! Don't forget to add // c_null_char to NRFFileName when using this interface
  interface
    subroutine c_interoperability_setupFSRMPointSources( momentTensor, numberOfSources, centres, strikes, dips, rakes, onsets, areas, timestep, numberOfSamples, timeHistories ) bind( C, name='c_interoperability_setupFSRMPointSources' )
      use iso_c_binding, only: c_double, c_int
      implicit none
      real(kind=c_double), dimension(*), intent(in) :: momentTensor
      integer(kind=c_int), value                    :: numberOfSources
      real(kind=c_double), dimension(*), intent(in) :: centres
      real(kind=c_double), dimension(*), intent(in) :: strikes
      real(kind=c_double), dimension(*), intent(in) :: dips
      real(kind=c_double), dimension(*), intent(in) :: rakes
      real(kind=c_double), dimension(*), intent(in) :: onsets
      real(kind=c_double), dimension(*), intent(in) :: areas
      real(kind=c_double), value                    :: timestep
      integer(kind=c_int), value                    :: numberOfSamples
      real(kind=c_double), dimension(*), intent(in) :: timeHistories
    end subroutine
  end interface

  interface c_interoperability_addReceiver
    subroutine c_interoperability_addReceiver( i_receiverId, i_meshId ) bind( C, name='c_interoperability_addReceiver' )
      use iso_c_binding, only: c_int
      implicit none
      integer(kind=c_int), value :: i_receiverId, i_meshId
    end subroutine
  end interface

  interface c_interoperability_setReceiverSampling
    subroutine c_interoperability_setReceiverSampling( i_receiverSampling ) bind( C, name='c_interoperability_setReceiverSampling' )
      use iso_c_binding, only: c_double
      implicit none
      real(kind=c_double), value :: i_receiverSampling
    end subroutine
  end interface

  interface c_interoperability_setMaterial
    subroutine c_interoperability_setMaterial( i_elem, i_side, i_materialVal, i_numMaterialVals ) bind( C, name='c_interoperability_setMaterial' )
    use iso_c_binding
    implicit none
    integer(kind=c_int), value :: i_elem
    integer(kind=c_int), value :: i_side
    real(kind=c_double), dimension(*), intent(in) :: i_materialVal
    integer(kind=c_int), value :: i_numMaterialVals
    end subroutine
  end interface

  interface c_interoperability_setInitialLoading
    subroutine c_interoperability_setInitialLoading( i_meshId, i_initialLoading ) bind( C, name='c_interoperability_setInitialLoading' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: i_initialLoading
    end subroutine
  end interface


  interface c_interoperability_setPlasticParameters
    subroutine c_interoperability_setPlasticParameters( i_meshId, i_plasticParameters ) bind( C, name='c_interoperability_setPlasticParameters' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: i_plasticParameters
    end subroutine
  end interface

  interface c_interoperability_initializeCellLocalMatrices
    subroutine c_interoperability_initializeCellLocalMatrices() bind( C, name='c_interoperability_initializeCellLocalMatrices' )
    end subroutine
  end interface

  interface c_interoperability_synchronizeCellLocalData
    subroutine c_interoperability_synchronizeCellLocalData() bind( C, name='c_interoperability_synchronizeCellLocalData' )
    end subroutine
  end interface

  interface c_interoperability_synchronizeCopyLayerDofs
    subroutine c_interoperability_synchronizeCopyLayerDofs() bind( C, name='c_interoperability_synchronizeCopyLayerDofs' )
    end subroutine
  end interface

  interface c_interoperability_enableDynamicRupture
    subroutine c_interoperability_enableDynamicRupture() bind( C, name='c_interoperability_enableDynamicRupture' )
    end subroutine
  end interface

  interface
    subroutine c_interoperability_enableWaveFieldOutput( i_waveFieldInterval, i_waveFieldFilename ) bind( C, name='c_interoperability_enableWaveFieldOutput' )
      use iso_c_binding
      implicit none
      real(kind=c_double), value :: i_waveFieldInterval
      character(kind=c_char), dimension(*), intent(in) :: i_waveFieldFilename
    end subroutine

    subroutine c_interoperability_enableCheckPointing( i_checkPointInterval, i_checkPointFilename, i_checkPointBackend ) bind( C, name='c_interoperability_enableCheckPointing' )
      use iso_c_binding
      implicit none
      real(kind=c_double), value :: i_checkPointInterval
      character(kind=c_char), dimension(*), intent(in) :: i_checkPointFilename
      character(kind=c_char), dimension(*), intent(in) :: i_checkPointBackend
    end subroutine

    subroutine c_interoperability_initializeIO( i_mu, i_slipRate1, i_slipRate2, i_slip, i_slip1, i_slip2, i_state, i_strength, &
        i_numSides, i_numBndGP, i_refinement, i_outputMask, i_outputRegionBounds ) &
        bind( C, name='c_interoperability_initializeIO' )
      use iso_c_binding
      implicit none

      real(kind=c_double), dimension(*), intent(in) :: i_mu
      real(kind=c_double), dimension(*), intent(in) :: i_slipRate1
      real(kind=c_double), dimension(*), intent(in) :: i_slipRate2
      real(kind=c_double), dimension(*), intent(in) :: i_slip
      real(kind=c_double), dimension(*), intent(in) :: i_slip1
      real(kind=c_double), dimension(*), intent(in) :: i_slip2
      real(kind=c_double), dimension(*), intent(in) :: i_state
      real(kind=c_double), dimension(*), intent(in) :: i_strength
      real(kind=c_double), dimension(*), intent(out) :: i_outputRegionBounds
      integer(kind=c_int), dimension(*), intent(out) :: i_outputMask
      integer(kind=c_int), value                    :: i_numSides
      integer(kind=c_int), value                    :: i_numBndGP
      integer(kind=c_int), value                    :: i_refinement
    end subroutine

    subroutine c_interoperability_addToDofs( i_meshId, i_update, numberOfQuantities ) bind( C, name='c_interoperability_addToDofs' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value                    :: i_meshId
      real(kind=c_double), dimension(*), intent(in) :: i_update
      integer(kind=c_int), value                    :: numberOfQuantities
    end subroutine

    subroutine c_interoperability_getTimeDerivatives( i_meshId, o_timeDerivatives ) bind( C, name='c_interoperability_getTimeDerivatives' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value :: i_meshId
      real(kind=c_double), dimension(*), intent(out) :: o_timeDerivatives
    end subroutine

    subroutine c_interoperability_getFaceDerInt( i_meshId, i_faceId, i_timeStepWidth, o_timeDerivativesCell, o_timeDerivativesNeighbor, o_timeIntegratedCell, o_timeIntegratedNeighbor ) bind( C, name='c_interoperability_getFaceDerInt' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value :: i_meshId
      integer(kind=c_int), value :: i_faceId
      real(kind=c_double), value :: i_timeStepWidth
      real(kind=c_double), dimension(*), intent(out) :: o_timeDerivativesCell
      real(kind=c_double), dimension(*), intent(out) :: o_timeDerivativesNeighbor
      real(kind=c_double), dimension(*), intent(out) :: o_timeIntegratedCell
      real(kind=c_double), dimension(*), intent(out) :: o_timeIntegratedNeighbor
    end subroutine

    subroutine c_interoperability_getDofs( i_meshId, o_dofs ) bind( C, name='c_interoperability_getDofs' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value                     :: i_meshId
      real(kind=c_double), dimension(*), intent(out) :: o_dofs
    end subroutine

    subroutine c_interoperability_getDofsFromDerivatives( i_meshId, o_dofs ) bind( C, name='c_interoperability_getDofsFromDerivatives' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value :: i_meshId
      real(kind=c_double), dimension(*), intent(out) :: o_dofs
    end subroutine

    subroutine c_interoperability_getNeighborDofsFromDerivatives( i_meshId, i_faceId, o_dofs ) bind( C, name='c_interoperability_getNeighborDofsFromDerivatives' )
      use iso_c_binding
      implicit none
      integer(kind=c_int), value :: i_meshId
      integer(kind=c_int), value :: i_faceId
      real(kind=c_double), dimension(*), intent(out) :: o_dofs
    end subroutine
  end interface

  interface c_interoperability_simulate
    subroutine c_interoperability_simulate( i_finalTime ) bind( C, name='c_interoperability_simulate' )
      use iso_c_binding, only: c_double
      implicit none
      real(kind=c_double), value :: i_finalTime
    end subroutine
  end interface

  interface c_interoperability_finalizeIO
    subroutine c_interoperability_finalizeIO() bind( C, name='c_interoperability_finalizeIO' )
    end subroutine
  end interface
end module
