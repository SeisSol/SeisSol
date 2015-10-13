!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) 2006-2014, SeisSol Group
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

MODULE ini_OptionalFields_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE ini_OptionalFields
     MODULE PROCEDURE ini_OptionalFields
  END INTERFACE
  INTERFACE close_OptionalFields
     MODULE PROCEDURE close_OptionalFields
  END INTERFACE
  
  !----------------------------------------------------------------------------
  PUBLIC  :: ini_OptionalFields
  PUBLIC  :: close_OptionalFields
  !----------------------------------------------------------------------------

CONTAINS
  
  SUBROUTINE ini_OptionalFields(OptionalFields,SOURCE,EQN,MESH,DISC,IO)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tUnstructOptionalFields)     :: OptionalFields
    TYPE (tSource)                     :: SOURCE
    TYPE (tEquations)                  :: EQN
    TYPE (tUnstructMesh)               :: MESH
    TYPE (tDiscretization)             :: DISC
    TYPE (tInputOutput)                :: IO
    !local Variables
    INTEGER                            :: allocstat(20)
    INTEGER                            :: RK_size
    ! -------------------------------------------------------------------------
    INTENT(IN)                         :: SOURCE,EQN,MESH,DISC,IO
    INTENT(INOUT)                      :: OptionalFields
    ! -------------------------------------------------------------------------
    !                                                                          
    allocstat(:) = 0                                                           
    !                                                                          
    ! Je nach Verfahren, Gleichung oder Zusatzfeatures   
    ! werden hier optionale Felder angelegt.                                                         
    !                                                                           
    ! --- Zeitschritt in jeder Zelle ----------------------------------------- 
    !                                                                          
    ALLOCATE(OptionalFields%dt(MESH%nElem),    & 
             STAT=allocstat(1)                 ) 
    OptionalFields%dt(:) = 0.                                                  
    !                                                                           
    ! --- Hintergrundswerte f�r ALLE linearen Gleichungen -------------------- 
    !                                                                          
    IF(EQN%linearized) THEN                                                   
       ALLOCATE(OptionalFields%BackgroundValue(MESH%nElem,EQN%nBackgroundVar),& 
                STAT=allocStat(4)                                          ) 
       OptionalFields%BackgroundValue(:,:) = 0.                                           
    ELSE
       NULLIFY(OptionalFields%BackgroundValue) 
    ENDIF                                                                             
    !                                                                        !                                                                          
    ! ========================================================================  
    !                 One common error handler for all allocations             
    ! ========================================================================  
    !                                                                          
    IF (SUM(ABS(allocStat(:))).NE.0) THEN                                       
       logError(*) 'Allocation error in ini_opt_fields. '
       logError(*) 'Status: ', allocstat(:)
       STOP                                                                    
    END IF                                                                              
    !
  END SUBROUTINE ini_OptionalFields

   SUBROUTINE close_OptionalFields(OptionalFields,SOURCE,EQN,MESH,DISC,IO)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tUnstructOptionalFields)     :: OptionalFields
    TYPE (tSource)                     :: SOURCE
    TYPE (tEquations)                  :: EQN
    TYPE (tUnstructMesh)               :: MESH
    TYPE (tDiscretization)             :: DISC
    TYPE (tInputOutput)                :: IO
    ! -------------------------------------------------------------------------
    INTENT(IN)                         :: SOURCE,EQN,MESH,DISC,IO
    INTENT(INOUT)                      :: OptionalFields
    ! -------------------------------------------------------------------------
    !                                                                         
    !                                                                           
    ! --- Zeitschritt in jeder Zelle ----------------------------------------- 
    !                                                                          
    DEALLOCATE(    OptionalFields%dt                                         )
    !                                                                           
    ! --- Hintergrundswerte f�r ALLE linearen Gleichungen -------------------- 
    !                                                                          
    IF(EQN%linearized) THEN                                                    
       DEALLOCATE( OptionalFields%BackgroundValue                            )
    ENDIF                                                                             
    !
  END SUBROUTINE close_OptionalFields
 
END MODULE ini_OptionalFields_mod
