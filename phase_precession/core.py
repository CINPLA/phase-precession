import numpy as np
from scipy.optimize import fminbound
from pycircstat import rayleigh, corrcc, corrcl


def goodness(x, phase, slope):
    '''
    Determines the goodness of a fit of data in the cyclic domain
    with respect to some linear function of a specific slope
    '''
    g_cos, g_sin = summed_vector(x, phase, slope)
    return np.sqrt(g_cos**2 + g_sin**2) / len(x)


def summed_vector(x, phase, slope):
    '''
    Returns the two components of a vector,
    which is the summed difference between a linear function and the data points.
    '''
    ph = phase - 2 * np.pi * x * slope  # phase difference
    g_cos = np.sum(np.cos(ph))       # sum up phase differences
    g_sin = np.sum(np.sin(ph))
    return g_cos, g_sin


def cl_regression(x, phase, min_slope, max_slope):
    '''Determines the best linear fit to data on the surface of a cylinder

    Translated to python from Richard Kempter (October-5, 2007)
    (December-1, 2019) by Mikkel Lepperød

    Parameters
    ----------
    x : array
        real-valued vector of `linear' instances
        (e.g. places, attenuations, frequencies, etc.)
    phase : array
        vector of phases at instances x in rad (values need NOT be
       restricted to the interval [0, 2pi)
    [min_slope, max_slope]:  float
        interval of slopes in which the best_slope is determined.
        In contrast to linear regression, we MUST restrict the range of
        slopes to some `useful' range determined by prior knowledge.
         ATTENTION ! because this is a possible source of errors,
        in particular for low numbers of data points.

    Returns
    -------
    R : float
        mean resultant lenght of the residual distribution
        is a measure of the goodness of the fit
        (also called vector strength).
        Small R indicates a bad fit (worst case R=0)
        Large R indicates a good fit (best case R=1)
    slope :  float
        slope (at the maximum of R) within the interval [min_slope, max_slope]
    phi0: float
        initial phase (or phase offset) of a cyclic regression line;
        values of phi0 are always restricted to the interval [0, 2pi]
    '''
    if len(x) != len(phase):
          raise ValueError('The lengths of x and phase must match.')

    if len(x) < 2:
          raise ValueError('The length of x is too small: len(x) < 2.')

    if not isinstance(min_slope, (float, int)):
       raise ValueError('The ''min_slope'' parameter must be a scalar');

    if not isinstance(max_slope, (float, int)):
       raise ValueError('The ''max_slope'' parameter must be a scalar');

    assert min_slope < max_slope, 'min_slope < max_slope'


    # We determine the value of the best `slope' using the MATLAB function fminbnd.
    # Please note that we have a minus sign in front of the
    # function `goodness' because we need the maximum (not the minimum)
    def func(opt_slope):
        return - goodness(x, phase, opt_slope)

    slope = fminbound(func, min_slope, max_slope)

    # Given the best `slope' we can explicitly calculate the goodness of the
    # fit ...
    R = goodness(x, phase, slope)

    # ... and the phase offset phi0:
    g_cos, g_sin = summed_vector(x, phase, slope)

    phi0 = np.arctan2(g_sin, g_cos)
    if phi0 < 0:   # restrict phases to the interval [0, 2pi]
      phi0 = phi0 + 2 * np.pi

    return phi0, slope, R


def cl_corr(x, phase, min_slope, max_slope, ci=.05, bootstrap_iter=1000):
    '''
    Function to (1) fit a line to circular-linear data and (2) determine
    the circular-linear correlation coefficient

    Parameters
    ----------
    x : array
        real-valued vector of `linear' instances
        (e.g. places, attenuations, frequencies, etc.)
    phase : array
        vector of phases at instances x in rad (values need NOT be
       restricted to the interval [0, 2pi)
    [min_slope, max_slope]:  float
        interval of slopes in which the best_slope is determined.
        In contrast to linear regression, we MUST restrict the range of
        slopes to some `useful' range determined by prior knowledge.
         ATTENTION ! because this is a possible source of errors,
        in particular for low numbers of data points.
    ci : float
        level of confidence desired, e.g. .05 for 95 % confidence
    bootstrap_iter : int
        number of bootstrap iterations (number of samples if None)

    See also
    --------
    pycircstat.corrcc

    Returns
    -------
        circ_lin_corr : float
            circular-linear correlation coefficient
        ci : array
            confidence interval
        slope : float
            slope of the fitted line in rad
        phi0_deg : float
            phase offset of the fitted line in deg
        RR : float
            goodness of fit
    '''
    phi0, slope, RR = cl_regression(x, phase, min_slope, max_slope) # fit line to data
    circ_x = np.mod(2 * np.pi * abs(slope) * x, 2 * np.pi) # convert linear variable to circular one

    circ_lin_corr, ci_out = corrcc(circ_x, phase, ci=ci, bootstrap_iter=bootstrap_iter)
    # circ_lin_corr, ci_out = corrcl(phase, x, ci=ci, bootstrap_iter=bootstrap_iter)

    return circ_lin_corr, ci_out, slope, phi0, RR
