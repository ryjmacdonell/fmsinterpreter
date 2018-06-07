"""
Module for determining branching between geometrically similar conical
intersections.

Based largely on the Kabsch geometry matching routine in GeomTools.
"""
import os
import numpy as np
from geomtools import molecule
from geomtools import kabsch
from fmsinterpreter import fileio
from fmsinterpreter import populations


def branching(reffnames, trajfnames, states, plist=[], invert=True):
    """Gets the conical intersection branching ratios by matching to a
    set of reference geometries."""
    refs = molecule.import_bundle(reffnames)
    testfnames = np.array([fn.replace('TrajDump', 'Spawn') for
                           fn in trajfnames])
    #for mol, fname in zip(tests.molecules, test_fnames):
    #    # don't think this is actually working now...
    #    mol.set_comment(fname + mol.get_comment())

    traj_info = populations.read_tjinfo(trajfnames)
    mask = populations.state_mask(traj_info, states)
    testfnames = testfnames[mask]
    tests = molecule.import_bundle(testfnames, hascom=True)

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


def trim_states(bundle, states):
    """Returns a MoleculeBundle only containing only the geometries with the
    correct state labels."""
    comments = [mol.comment.split() for mol in bundle.molecules]
    istate = np.array([int(com[9])-1 for com in comments])
    fstate = np.array([int(com[3])-1 for com in comments])
    mask = np.logical_and(istate != states[0], fstate != states[1])

    bundle.rm_molecules(np.where(mask))
    return np.logical_not(mask)
