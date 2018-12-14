"""
Module for determining branching between geometrically similar conical
intersections.

Based largely on the Kabsch geometry matching routine in GeomTools.
"""
import os
import numpy as np
import geomtools.molecule as gm
from fmsinterpreter import fileio
from fmsinterpreter import populations


def branching(reffnames, trajfnames, states, plist=None):
    """Gets the conical intersection branching ratios by matching to a
    set of reference geometries."""
    refs = gm.import_bundle(reffnames)
    testfnames = np.array([fn.replace('TrajDump', 'Spawn') for
                           fn in trajfnames])

    traj_info = populations.read_tjinfo(trajfnames)
    mask = populations.state_mask(traj_info, states)
    testfnames = testfnames[mask]
    tests = gm.import_bundle(testfnames, hascom=True)

    amps, times = populations.read_amps(trajfnames)
    pops = populations.integ_spawn_pop(amps, times, traj_info)
    pops = pops[mask]
    pops[pops < 0] = 0

    for mol, p in zip(tests.molecules, pops):
        t = mol.get_comment().split()[0]
        mol.set_comment('time = {:s}  amp = {:.4e}'.format(t, p))

    matched_geoms = tests.match_to_ref(refs, weighted=True,
                                       plist=plist, invert=invert)
    return matched_geoms


def trim_states(bundle, state0, state1):
    """Returns a MoleculeBundle only containing only the geometries with the
    correct state labels."""
    comments = [mol.comment.split() for mol in bundle.molecules]
    istate = np.array([int(com[9])-1 for com in comments])
    fstate = np.array([int(com[3])-1 for com in comments])
    mask = np.logical_and(istate != state0, fstate != state1)

    bundle.rm_molecules(np.where(mask))
    return bundle, np.logical_not(mask)
