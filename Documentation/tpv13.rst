SCEC TPV13
==========

TPV13 is similar with TPV12 except for that material properties are **non-associative Drucker-Prager plastic**. TPV13 is using non-associative Drucker-Prager plasticity with yielding in shear. The material is characterized by six constitutive parameters:

 *Bulk friction = 0.85*
  *Fluid pressure = 1000 kg/m3*
  Based on the SCEC benchmark, the stress components
  
  :math:`\sigma_{ij}` are defined as :
  
  *Mean stress*:
  
  :math:`\sigma_m = (\sigma_{11}+\sigma_{22}+\sigma_{33})/3`
  
  Stress deviator: :math:`s_{ij} = \sigma_{ij} - \sigma_m \delta_{ij}`
  Second invariant of the stress deviator:
  
  :math:`J_2(\sigma) = 1/2 *\sum_{ij} s_{ij} s_{ji}`
  
  Drucker-Prager yield stress:
  
  :math:`Y(\sigma) =\max(0,c\cos \phi - (\sigma_m +P_f)\sin \phi)`
  
  Drucker-Prager yield function:
  
  :math:`F(\sigma)=\sqrt{J_s(\sigma)-Y(\sigma)}`

The Drucker-Prager material is required to satisfy the yield equation:

  :math:`F(\sigma)\leq 0`
  
When :math:`F(\sigma) < 0`, the material behaves like a linear isotropic elastic material, 
with Lame paramters :math:`\lambda` and  :math:`\mu`.

Wen :math:`F(\sigma) = 0`, if the material is subjected to a strain that 
tends to cause an increase in :math:`F(\sigma)`, then the material
yields. For TPV13, we assume that the material yields in shear. Yielding
in shear means that when the material yields, the stress tensor
:math:`\sigma_{ij}` changes by an amount proportional to the stress
deviator :math:`s_{ij}`, so as to preserve the condition
:math:`F(\sigma)` with no change in mean stress :math:`\sigma_m` .

Nucleation
~~~~~~~~~~

TPV13 uses the same nucleation method as TPV12

Plasticity parameters
~~~~~~~~~~~~~~~~~~~~~

To turn on plasticity in SeisSol, add the following lines in
*parameter.par*:


&SourceType
Plasticity = 1 ! default = 0
Tv = 0.03 ! Plastic relaxation
/


In the **material.yaml**, add plasticity parameters:

.. literalinclude:: tpv12/tpv12_13_material.yaml
   :language: fortran


Results
~~~~~~~

Figure [fig:tpv13compare] shows the comparison between TPV12 (elastic)
and TPV13 (plastic). The peak of slip rate in TPV12 is higher than
TPV13. This difference attributes to the response of the off-fault
plasticity. Refer to Wollherr et al. [2018] for detailed
discussions.

.. figure:: LatexFigures/SRs_12_13.png
   :alt: Diagram of along-strike slip rate
   :width: 12.00000cm

   Diagram of along-strike slip rate (left) and along-dip slip rate
   (right) in TPV12 (blue) and TPV13 (orange). 
   
