Slip-rate imposed on a DR boundary
===================================

This "pseudo" friction law allows imposing slip-rate on a dynamic rupture boundary.
The FL id for this friction law is 33.
The advantage of this approach compared to a multi point-source representation is that the fault slip is not condensed to points. 
Therefore the discontinuity of the displacement across the fault can be accurately accounted for.
In addition, no spurious waves, related to the multi point-sources representation, are generated.

The current implementation allows imposing kinematic models parameterized by regularized Yoffe function (see Tinti et al., 2005, doi:10.1029/2005JB003644).
The Yoffe functions are parametrized by ``rupture_onset``, ``acc_time`` and ``effective_rise_time``, where ``rupture_onset`` is the onset time of the rupture, 
``acc_time`` is the time to the peak slip-rate, and ``effective_rise_time`` is the duration of slip.
The slip distribution is defined using easi by the ``strike_slip`` and ``dip_slip`` variables.  

Warning: the direction of positive ``strike_slip`` and ``dip_slip`` is based on the convention of Seissol (e.g. positive strike_slip for right-lateral faulting).   


A fully working example based on Northridge is given at https://github.com/SeisSol/Examples/tree/master/Northridge_FL33.
