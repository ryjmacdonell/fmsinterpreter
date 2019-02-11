"""
Module for interpreting internal coordinate densities.

At the moment, this is meant to use Monte-Carlo density files generated
by mcdensity.f90.
"""
import numpy as np
from fmsinterpreter import fileio


def read_density(stub, tinc, tmax, cnum, cmax, tmin=0., cmin=0.):
    """Reads the mcdensity output files."""
    times = np.arange(0, tmax+0.1*tinc, tinc)
    cbin = np.linspace(cmin, cmax, cnum+1)
    coord = (cbin[1:] + cbin[:-1]) / 2
    dens = np.zeros((len(times), len(coord)))

    for i, t in enumerate(times):
        fname = stub + '{:.1f}'.format(t)
        try:
            rawc, rawd = fileio.read_dat(fname, skiprow=1).T
        except IndexError:
            # this is messy, but it catches files with no density
            continue
        dens[i] = np.histogram(rawc, bins=cbin, weights=rawd)[0]

    return times, coord, dens
