"""
Plotting routine for visualizing critical points on a potential energy
surface relative to a given minimum energy.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import fileio as fileio


def read_nrg(fname):
    """Reads the energy file and returns the data as a list of labels
    and an array of energies.

    File format:
    lbl1 e01 e11 e21 ...
    lbl2 e02 e12 e22 ...
    lbl3 e03 e13 e23 ...
    """
    with open(fname, 'r') as f:
        lines = f.readlines()
        lbl = [line.split()[0] for line in lines]
        e = np.array([line.split()[1:] for line in lines], dtype=float)

    return lbl, e


def conv_nrg(lbl, e, order=None, conv=1., base=None, states=None):
    """Reorders labels and energies according to pltorder and converts
    energies to desired units."""
    if order is None:
        order = range(len(lbl))
    if base is None:
        base = np.amin(e)
    if states is None:
        states = range(len(e[0]))

    new_lbl = [lbl[i] for i in order]
    new_e = conv * (e[order][:,states] - base)

    return new_lbl, new_e


def plot_nrg(ax, lbl, e, wid=1, sep=1, rot=90, maxe=None):
    """Plots the potential energy surface at critical geometries as a
    set of bars with connecting lines."""
    eplot = np.repeat(e, 2, axis=0)
    x1 = np.arange(len(e)) * (wid + sep)
    x2 = np.insert(x1 + wid, range(len(x1)), x1)

    ax.plot(x2, eplot, zorder=0)
    for i in range(len(e[0])):
        ax.hlines(e[:,i], x1, x1 + wid, linewidth=2)

    ax.set_xticks(x1 + 0.5*wid)
    ax.set_xticklabels(lbl, rotation=rot)

    ax.set_xlim(x2[0] - wid, x2[-1] + wid)
    if maxe is None:
        ax.set_ylim(-0.1)
    else:
        ax.set_ylim(-0.1, maxe)
    ax.set_ylabel('Energy / eV')

    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')

    return ax
