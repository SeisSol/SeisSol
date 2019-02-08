Attenuation
===========

General comments
----------------

-  Q is measuring the decay of wave amplitude during its propagation in
   the medium and is physically defined as the ratio of wave energy to
   the energy dissipated per cycle of oscillation
-  The actual dissipation mechanisms, as scattering,heat generation,
   frictional losses in the vibrating crystal lattice, etc., "causing Q"
   can be manifold and do not have to follow similar physics
-  Q is not very frequency dependent over a large range of (lower)
   frequencies but observations are ambiguous at higher frequencies
-  Common choices are Qs ~ 0.5 Qp and Qs as 40-50*Vs (in km/s) (e.g.
   Olsen et al., 2009)(given Qs = 100-200 for shallow sediments)
-  In damaged fault zones Qs would be expected to be as low as on the
   order of 10 - and potentially smear out slip rate on the fault
   without affecting rupture speeds

Numerical implementing in SeisSol
---------------------------------

-  Implementing attenuation in SeisSol follows the ideas of a general
   Maxwell or Zener body
-  Each damping mechanism can be parametrized by its own relaxation
   frequency
-  It aims at resolving a frequency independent Q with an adequate
   number of anelastic functions and relaxation frequencies to cover the
   frequency range under interest. Usually, 3 Maxwell bodies is enough
   for 5% error.

Local time stepping and CFL
---------------------------

Using attenuation and local time stepping (LTS) seems to strongly impact
the stability. The CFL might need to be strongly decreased to recover
stability. This may remove most of the speed-up of the LTS. Examples:

-  Two order 3 elastic simulations (different meshes) could run with a
   CFL of 0.8 using LTS. The CFL had to be decreased to 0.15 (resp. 0.2)
   for a stable run using LTS and attenuation.
-  Using one of the two previous meshes, an order 5 elastic simulations
   could run with a CFL of 0.6 using LTS. The CFL had to be decreased to
   0.3 for a stable run using LTS and attenuation.

Compiling
---------

The attenuation implementation requires libxsmm. To install it:

-  ``git clone https://github.com/hfp/libxsmm``
-  ``cd libxsmm``
-  ``make generator``
   (on supermuc, the python module should be loaded)
-  Add libxsmm/bin to the PATH environment variable
-  Add in SeisSol build configuration script:
   ``equations = 'viscoelastic2'``
   ``numberOfMechanisms = 3``

Dispersion
----------

| The attenuation implementation implies dispersion in P and S wave
  velocities: that why we need to define a central frequency at which
  Vp/Vs are exact. In addition, the effective Q values are not exactly
  equal to the desired Q, but are oscillating around those values. The
  variation of Vp, Vs, Qp and Qs with frequency can be visualized using
  this script:
| `https://github.com/SeisSol/SeisSol/blob/master/preprocessing/science/ViscoelasticModComp.m <https://github.com/SeisSol/SeisSol/blob/master/preprocessing/science/ViscoelasticModComp.m>`__

Additionnal parameters
----------------------

| The relaxation frequencies have to be logarithmically equispaced, i.e.
  *log(wi+1)-log(wi) = constant*. In the parameter file one has to give
  a frequency ratio of maximum to minimum frequency and a central
  frequency. For example, in the case of 3 mechanisms the following
  relations define the relaxation frequencies:
| *w2 = FreqCentral*
| *log(w3)-log(w2) = log(w2) - log(w1)*
| *w3 / w1 = FreqRatio*
| Outside of the frequency band *w1 - w3*, Q goes to infinity, yielding
  elastic behaviour.

easi branch
~~~~~~~~~~~

| Add Inside the parameter file of SeisSol, in the '&equations' section
  (frequencies values to be adapted to the source frequency content):
| ``FreqCentral=2.5``
| ``FreqRatio=100``

The spatial variation of Qs and Qp are defined with easi in the
MaterialFileName. Here is an example of easi file, in which Qs and Qp
are directly related to the shear wave speed Vs:

.. code:: cpp

   !ASAGI
   file: ../material/vmodel_500.nc
   parameters: [rho, mu, lambda]
   var: data
   components: !FunctionMap
     map:
       rho:    return rho;
       mu:     return mu;
       lambda: return lambda;
       Qs:     return 0.1 * sqrt(mu/rho);
       Qp:     return 0.2 * sqrt(mu/rho);

master branch
~~~~~~~~~~~~~

-  | Add Inside the parameter file of SeisSol, in the '&equations'
     section:
   | ``Anelasticity=1``
   | ``nMechanisms = 3``
   | ``FreqCentral=2.5``
   | ``FreqRatio=100``
   | ``MaterialFileName = test.def``

-  | The material file (here test.def) is optional. Here we can assign
     the material values rho, mu, lambda, as well as Qp and Qs for each
     layer of the mesh (specified by the mesher).
   | ``1          !number of layers``
   | ``3          ! Mechanisms``
   | ``2.5        ! Central frequency``
   | ``100.       ! Frequency Ratio``
   | ``1 2670.0 3.203812032e10 3.204375936e10 69.3 155.9 !rho, mu, lambda, Qp, Qs``

.. _hard-coded-attenuation-model-(master-branch):

Hard-coded attenuation model (master branch)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another approach is to add an attenuation model in the source code
directly (in case we don't have specified layers in the mesh or just
want to apply a layer-independent attenuation model) example (in
ini_model.f90)

.. code:: fortran

           if (EQN%Anelasticity.EQ.1) THEN
              DO iElem=1, MESH%nElem
                 MaterialTmp(:) = EQN%MODEL(1,:)
                 MaterialVal(iElem,1:3) = MaterialTmp(1:3)
                 EQN%LocAnelastic(iElem) = 1                                        ! Mark element with anelastic material
                 CALL ini_ATTENUATION(Theta,w_freq,Material_INF,MaterialTmp,EQN)    ! Initialize anelastic coefficients for this zone     
                 MaterialVal(iElem,2:EQN%AneMatIni-1) = Material_INF(:)             ! Set unrelaxed material properties for this zone.                                                                      !
                 ! Fill MaterialVal vector for each element with anelastic coefficients w_freq and theta 
                 DO iMech = 1, EQN%nMechanisms
                    MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1))             = w_freq(iMech)
                    MaterialVal(iElem,EQN%AneMatIni+4*(iMech-1)+1:EQN%AneMatIni+4*(iMech-1)+3) = Theta(iMech,:)
                 ENDDO
              ENDDO
           ELSE
              MaterialVal(:,1) = EQN%rho0
              MaterialVal(:,2) = EQN%mu
              MaterialVal(:,3) = EQN%lambda
           ENDIF

.. _hard-coded-in-readpar.f90-(master-branch):

Hard-coded in readpar.f90 (master branch)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

example

.. code:: fortran

       CASE(XXX)
         if (EQN%Anelasticity.EQ.1) THEN
         logInfo0(*) 'Material properties are read from file : ', TRIM(EQN%MaterialFileName)
         CALL OpenFile(                                        &
               UnitNr       = IO%UNIT%other01                , &
               Name         = EQN%MaterialFileName           , &
               create       = .FALSE.                          )
         logInfo(*) 'Reading material property file ...  '
         READ(IO%UNIT%other01,'(i10,a)') EQN%nLayers, cdummy             ! Number of different material zones
         READ(IO%UNIT%other01,'(i10,a)') EQN%nMechanisms, cdummy         ! Number of different attenuation mechanisms
         logInfo(*) 'Model has ',EQN%nMechanisms,' attenuation mechanisms.'
         READ(IO%UNIT%other01,*) EQN%FreqCentral                             ! Central frequency of the absorption band (in Hertz)
         logInfo(*) 'with central frequency ',EQN%FreqCentral
         READ(IO%UNIT%other01,*) EQN%FreqRatio                               ! The ratio between the maximum and minimum frequencies of our bandwidth
         logInfo(*) 'and frequency ratio ',EQN%FreqRatio

         EQN%nBackgroundVar  = 3 + EQN%nMechanisms * 4
         EQN%nAneMaterialVar = 5        ! rho, mu, lambda, Qp, Qs
         EQN%nVarTotal = EQN%nVar + EQN%nAneFuncperMech*EQN%nMechanisms                                                    !
         EQN%AneMatIni = 4                                                  ! indicates where in MaterialVal begin the anelastic parameters 

         ALLOCATE(EQN%MODEL(1:EQN%nLayers,EQN%nAneMaterialVar))
         DO i = 1,EQN%nLayers
              READ(IO%UNIT%other01,*) intDummy, EQN%MODEL(i,:)
         ENDDO
         CLOSE(IO%UNIT%other01)
         else
           logInfo0(*) 'Jacobians are globally constant with rho0, mu, lambda:'
           logInfo0(*) ' rho0 = ', EQN%rho0     ! (1)
           logInfo0(*) ' mu = ', EQN%mu       ! (2)
           logInfo0(*) ' lambda = ', EQN%lambda   ! (3)
         endif

