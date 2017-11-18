"""
General plotting routines for the FMS scripts.

This should be a stand-alone module which handles all calls to
matplotlib. The Figure object can be used to handle multiple
plots on a single or multiple canvases.
"""
import numpy as np
import matplotlib.pyplot as plt


def Figure(object):
    """An object which handles matplotlib objects needed for figures."""
    def __init__(self, subplots=(1,1), **kwargs):
        self.subplots = subplots
        self.fig, axarr = plt.subplots(*subplots, **kwargs)
        self.axarr = np.atleast_2d(axarr)

    def scatter(self, x, y, isub=(0,0), **kwargs):
        """Plots a scatter plot of x and y positions in the subplot
        position isub."""
        ax = self.axarr[isub[0], isub[1]]
        ax.scatter(x, y, **kwargs)

    def heatmap(self, x, y, z, isub=(0,0), **kwargs):
        """Plots a heatmap of z vs. x and y in the subplot position isub."""
        ax = self.axarr[isub[0], isub[1]]
        ax.pcolormesh(x, y, z, **kwargs)

    def lineplot(self, x, y, isub=(0,0), **kwargs):
        """Plots a line from data y vs. x in the subplot position isub."""
        ax = self.axarr[isub[0], isub[1]]
        ax.plot(x, y, **kwargs)


def scatter(x, y, **kwargs):
    """Plots a scatter plot of x and y positions in a single frame."""
    plt.scatter(x, y, **kwargs)


def heatmap(x, y, z, **kwargs):
    """Plots a heatmap of z vs. x and y in a single frame."""
    plt.pcolormesh(x, y, z, **kwargs)


def lineplot(x, y, **kwargs):
    """Plots a line from data y vs. x in a single frame."""
    plt.plot(x, y, **kwargs)
