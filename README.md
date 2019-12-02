phase-precession
=================

Library for quantification of hippocampal phase precession translated to python from matlab scripts written by Richard Kempter

described in Kempter R, Leibold C, Buzsáki G, Diba K, Schmidt R (2012) Quantifying circular-linear associations: hippocampal phase precession. J Neurosci Methods 207:113–124. http://dx.doi.org/10.1016/j.jneumeth.2012.03.007


The main usage is

```python
from phase_precession import cl_corr
circ_lin_corr, ci, slope, phi0, RR = cl_corr(lin_x, circ_y_deg, min_slope, max_slope, 0.05, 10000)
```

Also, see Philipp Beren's CircStat toolbox, which this is based on.  https://github.com/circstat/pycircstat
