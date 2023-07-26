
# Auxiliary .py files used in the MetSta interface

import numpy as np
from scipy.stats import norm, chi2
import holoviews as hv

# Functions slightly changed from the elips.py file to be ran with plotly

def plot_confidence_ellipse(points, q=None, nstd=2):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        q : float, optional
            Confidence level, should be in (0, 1)
        nstd : int, optional
            Confidence level in unit of standard deviations. 
        E.g. 1 stands for 68.3% and 2 stands for 95.4%.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    if q != 0:
        q = np.asarray(q)
    elif nstd is not None:
        q = 2 * norm.cdf(nstd) - 1
    else:
        raise ValueError('One of `q` and `nsig` should be specified.')

    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, q)

def plot_cov_ellipse(cov, pos, q):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        q : float, optional
            Confidence level, should be in (0, 1)
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    r2 = chi2.ppf(q, 2)

    vals, vecs = eigsorted(cov)
    # Width and height are "full" widths, not radius
    # width, height = 2 * nstd * np.sqrt(vals)
    width, height = 2 * np.sqrt(vals * r2)
    theta = np.arctan2(*vecs[:,0][::-1])

    ellip = hv.Ellipse(pos[0], pos[1], (width,height), orientation=theta)

    return ellip