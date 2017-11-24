"""
Plotting routine for visualizing critical points on a potential energy
surface relative to a given minimum energy.
"""
import os
import numpy as np
from fmsinterpreter import fileio


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
