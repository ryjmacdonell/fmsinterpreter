"""
FMS analysis routine for visualizing the adiabatic population as a
function of time and fitting the population to a mono-, bi- or
triexponential curve.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from fmsinterpreter import fileio


def read_amps(fnames, times, states):
    """Reads the amplitudes of each trajectory within a time window
    and returns an array of total amplitudes."""
    amps = np.zeros((len(states), len(times)))

    for fname in fnames:
        with open(fname) as f:
            lines = [line.split() for line in f.readlines() if
                     'Time' not in line and line != '\n']

        alldata = np.array(lines, dtype=float)
        tlocal = alldata[:,0]
        st = states.index(int(alldata[0,-1]) - 1)
        if st in states:
            for i, t in enumerate(times):
                if t > tlocal[-1]:
                    amps[st, i:] += alldata[len(tlocal)-1, -2]
                    break
                elif t >= tlocal[0]:
                    amps[st, i] += alldata[np.argmin(np.abs(tlocal - t)), -2]

    return amps


def plot_amps(ax, times, amps, states, tconv=1.):
    """Plots the amplitudes as a function of time."""
    for a, s in zip(amps, states):
        ax.plot(times*tconv, a, label=r'S{:d}'.format(s))
    ax.set_xlim(min(times*tconv), max(times*tconv))
    ax.set_ylim(0, 1)
    ax.set_xlabel('Time / fs')
    ax.set_ylabel('Adiabatic population')
    ax.legend(loc=7)

    return ax


def fit_function(func, times, decay, p0, tconv=1.):
    """Fits amplitudes to a given exponential decay function.

    For now, this is just fitting the return to the ground state. With
    some sort of string parsing it could be changed to accept any
    state or combination of states (e.g. 'S1 + S2').
    """
    #decay = 1 - amps[0]
    #decay = (decay - min(decay)) / (max(decay) - min(decay))
    t = times * tconv
    popt, pcov = curve_fit(globals()[func], t, decay, p0=p0)
    perr = np.sqrt(np.diag(pcov))

    return popt, perr


def write_fit(func, popt, perr, outfname):
    """Writes fit information to an output file.

    This should be generalized to accept more than one set of fit values.
    """
    if func == 'exp_func':
        fitvals = ['t0', 'tau1']
    elif func == 'biexp_func':
        fitvals = ['t0', 'amp1', 'tau1', 'amp2', 'tau2']
    elif func == 'triexp_func':
        fitvals = ['t0', 'amp1', 'tau1', 'amp2', 'tau2', 'amp3', 'tau3']

    with open(outfname, 'w') as f:
        f.write('Curve  ')
        f.write(''.join(['{:>10s}'.format(fv) for fv in fitvals]) + '\n')
        f.write('1-S0   ')
        f.write(''.join(['{:10.4f}'.format(p) for p in popt]) + '\n')
        f.write('Error  ')
        f.write(''.join(['{:10.4f}'.format(p) for p in perr]) + '\n')


def plot_fit(ax, func, times, decay, popt, tconv=1.):
    """Plots the fit results for comparison to the raw data.

    This should be generalized to accept more than one set of fit values.
    """
    t = times * tconv
    ax.plot(t, decay, 'k-')
    ax.plot(t, globals()[func](t, *popt), 'r-')
    ax.set_ylim(0, 1)
    ax.set_xlabel('Time / fs')
    ax.set_ylabel('1-S0 Amplitude')

    return ax


def exp_func(x, x0, b):
    """Returns an exponential function for curve fitting purposes."""
    return np.exp(-(x - x0) / b)


def biexp_func(x, x0, a1, b1, a2, b2):
    """Returns a biexponential function for curve fitting purposes."""
    return a1 * np.exp(-(x - x0) / b1) + a2 * np.exp(-(x - x0) / b2)


def triexp_func(x, x0, a1, b1, a2, b2, a3, b3):
    """Returns a triexponential function for curve fitting purposes."""
    return (a1 * np.exp(-(x - x0) / b1) + a2 * np.exp(-(x - x0) / b2) +
            a3 * np.exp(-(x - x0) / b3))
