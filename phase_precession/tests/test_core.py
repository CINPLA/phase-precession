import pytest
import numpy as np

def test_cl_corr():
    from phase_precession import cl_corr
    x = np.linspace(0, 1, 100)
    phase = np.linspace(-np.pi, np.pi, 100)
    circ_lin_corr, ci_out, slope, phi0, RR = cl_corr(x, phase, 0, 1)
    assert circ_lin_corr - 1 < 1e-10
    assert slope - 1 < 1e-10
