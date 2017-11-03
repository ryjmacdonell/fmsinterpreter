"""
General plotting routines for the FMS scripts.

This should be a stand-alone script which handles all calls to
matplotlib. The Figure object can be used to handle multiple
plots on a single or multiple canvases.
"""
import numpy as np
import matplotlib.pyplot as plt


def Figure(object):
    """An object which handles matplotlib objects needed for figures."""
    def __init__(self, subplots=None, sharex=False, sharey=False):
        self.subplots = subplots
        if subplots is None:
            self.fig, self.ax = plt.subplots(sharex=sharex, sharey=sharey)
        else:
            self.fig, self.axarr = plt.subplots(sharex=sharex, sharey=sharey)

    def scatter():
        pass

    def heatmap():
        pass

    def lineplot():
        pass


def scatter():
    pass


def heatmap():
    pass


def lineplot():
    pass
