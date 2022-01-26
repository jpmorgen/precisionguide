"""Utilities for precisionguide system"""

import matplotlib.pyplot as plt
import numpy as np

def hist_of_im(im, binsize=1, show=False):
    """Returns a tuple of the histogram of image and index into *centers* of
bins."""
    # Code from west_aux.py, maskgen.
    # Histogram bin size should be related to readnoise
    hrange = (im.data.min(), im.data.max())
    nbins = int((hrange[1] - hrange[0]) / binsize)
    hist, edges = np.histogram(im, bins=nbins,
                               range=hrange, density=False)
    # Convert edges of histogram bins to centers
    centers = (edges[0:-1] + edges[1:])/2
    if show:
        plt.plot(centers, hist)
        plt.show()
        plt.close()
    return (hist, centers)

def iter_linfit(x, y, max_resid=None):
    """Performs least squares linear fit iteratively to discard bad points

    If you actually know the statistical weights on the points,
    just use polyfit directly.

    """
    # Let polyfit report errors in x and y
    coefs = np.polyfit(x, y, 1)
    # We are done if we have just two points
    if len(x) == 2:
        return coefs
        
    # Our first fit may be significantly pulled off by bad
    # point(s), particularly if the number of points is small.
    # Construct a repeat until loop the Python way with
    # while... break to iterate to squeeze bad points out with
    # low weights
    last_redchi2 = None
    iterations = 1
    while True:
        # Calculate weights roughly based on chi**2, but not going
        # to infinity
        yfit = x * coefs[0] + coefs[1]
        resid = (y - yfit)
        if resid.all == 0:
            break
        # Add 1 to avoid divide by zero error
        resid2 = resid**2 + 1
        # Use the residual as the variance + do the algebra
        redchi2 = np.sum(1/(resid2))
        coefs = np.polyfit(x, y, 1, w=1/resid2)
        # Converge to a reasonable epsilon
        if last_redchi2 and last_redchi2 - redchi2 < np.finfo(float).eps*10:
            break
        last_redchi2 = redchi2
        iterations += 1

    # The next level of cleanliness is to exclude any points above
    # max_resid from the fit (if specified)
    if max_resid is not None:
        goodc = np.where(np.abs(resid) < max_resid)
        # Where returns a tuple of arrays!
        if len(goodc[0]) >= 2:
            coefs = iter_linfit(x[goodc], y[goodc])
    return coefs
    
