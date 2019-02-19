Attenuation
===========

General comments
----------------

-  Q is measuring the decay of wave amplitude during its propagation in
   the medium and is physically defined as the ratio of wave energy to
   the energy dissipated per cycle of oscillation.
-  The actual dissipation mechanisms, as scattering,heat generation,
   frictional losses in the vibrating crystal lattice, etc., "causing Q"
   can be manifold and do not have to follow similar physics.
-  Q is not very frequency dependent over a large range of (lower)
   frequencies but observations are ambiguous at higher frequencies.
-  Common choices are :math:`Q_s \sim 0.5 Q_p` and :math:`Q_s \sim 40-50V_s` (in km/s) (e.g.
   Olsen et al., 2009). :math:`Q_s = 100-200` for shallow sediments.
-  In damaged fault zones :math:`Q_s` would be expected to be as low as on the
   order of 10 and potentially smear out slip rate on the fault
   without affecting rupture speeds.

Numerical implementing in SeisSol
---------------------------------

-  Implementing attenuation in SeisSol follows the ideas of a general
   Maxwell or Zener body.
-  Each damping mechanism can be parametrized by its own relaxation
   frequency.
-  It aims at resolving a frequency independent Q with an adequate
   number of anelastic functions and relaxation frequencies to cover the
   frequency range under interest. Usually, 3 Maxwell bodies is enough
   for 5% error.

Local time stepping and CFL
---------------------------

Using attenuation and local time stepping (LTS) seems to strongly impact
the stability. The CFL might need to be strongly decreased to recover
stability. This may remove most of the speed-up of the LTS. 

Examples:

.. |checkmark| unicode:: U+2713

.. list-table::
   :widths: 20 20 20 20
   :header-rows: 1

   * - order
     - mesh id
     - max CFL elastic
     - max CFL attenuation
   * - 3
     - 1
     - 0.8
     - 0.15
   * - 3
     - 2
     - 0.8
     - 0.2
   * - 5
     - 1
     - 0.6
     - 0.3

Nevertheless, it seems that the unstable elements are most of the time very large elements, rarely updated by the LTS.
Therefore, stability can probably be achieved without a too strong impact on the LTS speed-up by using the following patch, 
which degrade the LTS clustering in a more progressive way than a change of CFL.

.. code:: cpp

   diff --git a/src/Initializer/time_stepping/MultiRate.hpp b/src/Initializer/time_stepping/MultiRate.hpp
   index dce883a..6a2bfb9 100644
   --- a/src/Initializer/time_stepping/MultiRate.hpp
   +++ b/src/Initializer/time_stepping/MultiRate.hpp
   @@ -76,8 +76,11 @@ class seissol::initializers::time_stepping::MultiRate {
                                      unsigned int &o_clusterId ) {
          // first multi-rate interval
          double l_lower = i_minimumTimeStepWidth;
   +#if NUMBER_OF_QUANTITIES > 9
   +      double l_upper = i_multiRate*l_lower*(1.35*i_multiRate/2.);
   +#else
          double l_upper = i_multiRate*l_lower;
   -
   +#endif
          for( unsigned int l_id = 0; ; l_id++ ) {
            // the first cluster with an upper bound above the time step width is our
            if( l_upper > i_timeStepWidth ) {
   @@ -88,7 +91,11 @@ class seissol::initializers::time_stepping::MultiRate {

            // update interval and continue searching
            l_lower = l_upper;
   +#if NUMBER_OF_QUANTITIES > 9
   +        l_upper = i_multiRate * l_lower*(1.35*i_multiRate/2.);
   +#else
            l_upper = i_multiRate * l_lower;
   +#endif
          }
        }


Compiling
---------

The attenuation implementation requires libxsmm. To install it:

.. code:: c

   git clone https://github.com/hfp/libxsmm
   cd libxsmm
   make generator

(on supermuc, the python module should be loaded)

-  Add libxsmm/bin to the PATH environment variable
-  Add in SeisSol build configuration script:

.. code:: python

   equations = 'viscoelastic2'
   numberOfMechanisms = 3

Dispersion
----------

The attenuation implementation implies dispersion in P and S wave
velocities: that why we need to define a central frequency at which
:math:`V_p/V_s` are exact. In addition, the effective Q values are not exactly
equal to the desired Q, but are oscillating around those values. The
variation of :math:`V_p`, :math:`V_s`, :math:`Q_p` and :math:`Q_s` with frequency can be visualized using
`ViscoelasticModComp.m <https://github.com/SeisSol/SeisSol/blob/master/preprocessing/science/ViscoelasticModComp.m>`__.

Additionnal parameters
----------------------

| The relaxation frequencies have to be logarithmically equispaced, i.e.

| :math:`log(w_{i+1})-log(w_i) =` constant. 

In the parameter file one has to give a frequency ratio of maximum to minimum frequency and a central frequency. 
For example, in the case of 3 mechanisms the following relations define the relaxation frequencies:

| :math:`w_2 = FreqCentral`  

| :math:`log(w_3)-log(w_2) = log(w_2) - log(w_1)`  

| :math:`w_3 / w_1 = FreqRatio`  

Outside of the frequency band :math:`w_1 - w_3`, Q goes to infinity, yielding
elastic behaviour.

easi branch
~~~~~~~~~~~

Add Inside the parameter file of SeisSol, in the '&equations' section 
(frequencies values to be adapted to the source frequency content):

.. code:: fortran

   FreqCentral=2.5
   FreqRatio=100

The spatial variation of :math:`Q_s` and :math:`Q_p` are defined with easi in the
MaterialFileName. Here is an example of easi file, in which :math:`Q_s` and :math:`Q_p`
are directly related to the shear wave speed :math:`V_s`:

.. code:: yaml

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


Hardcoded_ini branch (outdated)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  | Add Inside the parameter file of SeisSol, in the '&equations'
     section:


.. code:: fortran

   Anelasticity=1
   nMechanisms = 3
   FreqCentral=2.5
   FreqRatio=100
   MaterialFileName = test.def

-  | The material file (here test.def) is optional. Here we can assign
     the material values rho, mu, lambda, as well as Qp and Qs for each
     layer of the mesh (specified by the mesher).

.. code:: fortran

   1          !number of layers
   3          ! Mechanisms
   2.5        ! Central frequency
   100.       ! Frequency Ratio
   1 2670.0 3.203812032e10 3.204375936e10 69.3 155.9 !rho, mu, lambda, Qp, Qs


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

