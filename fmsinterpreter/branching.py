"""
FMS analysis routine to read spawn geometries, find the closest matching
MECI and find the total population transferred through each MECI region.
"""
import os
import numpy as np
from geomtools import molecule
from geomtools import kabsch
import fileio as fileio


def pop_estim(fnamelist):
    """Returns a list of populations calculated from the squares of the
    norms of the amplitudes."""
    pass


def pop_accur(fnamelist):
    """Returns a list of populations calculated from the product A* S A
    for amplitudes A and overlaps S."""
    pass


def branching(reffnames, testfnames, states, pop_func=pop_estim):
    """Gets the conical intersection branching ratios by matching to a
    set of reference geometries."""
    refs = molecule.import_bundle(ref_fnames)
    tests = molecule.import_bundle(test_fnames, hascom=True)
    for mol, fname in zip(tests.molecules, test_fnames):
        # don't think this is actually working now...
        mol.set_comment(fname + mol.get_comment())

    tests = trim_states(tests, states)

    pop_func(testfnames)

    matched_geoms = tests.match_to_ref(refs, weighted=True,
                                       plist=inp['permute'],
                                       invert=inp['invert'])
    return matched_geoms # and populations


def trim_states(bundle, states):
    """Returns a MoleculeBundle only containing only the geometries with the
    correct state labels."""
    comments = [mol.comment.split() for mol in bundle.molecules]
    istate = np.array([int(com[9])-1 for com in comments])
    fstate = np.array([int(com[3])-1 for com in comments])
    mask = np.logical_and(istate != states[0], fstate != states[1])

    bundle.rm_molecules(np.where(mask))
    return bundle
