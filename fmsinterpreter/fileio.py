"""
Routines for reading in FMS files and configuration files.
"""
import os
import numpy as np
from glob import glob


def convert_str(string):
    """If possible, converts a string to an int, float, list of ints
    or list of floats."""
    if string == 'None':
        return None
    elif ';' in string and ',' in string:
        # 2D list delimited by ';' followed by ','
        slist = [ln.split(',') for ln in string.split(';')]
        try:
            return [[int(el) for el in ln] for ln in slist]
        except ValueError:
            pass
        try:
            return [[float(el) for el in ln] for ln in slist]
        except ValueError:
            return slist
    elif ';' in string or ',' in string:
        # 1D list delimited by either ',' or ';'
        slist = string.replace(';',',').split(',')
        try:
            return [int(el) for el in slist]
        except ValueError:
            pass
        try:
            return [float(el) for el in slist]
        except ValueError:
            return slist
    else:
        # not a list, try int or float
        try:
            return int(string)
        except ValueError:
            pass
        try:
            return float(string)
        except ValueError:
            return string


def read_cfg(fname, reqvars=[]):
    """Reads the configuration file."""
    fvars = dict()
    with open(fname, 'r') as f:
        lines = f.readlines()

    # remove spaces and comments and get inputs
    lines = [ln.partition('#')[0] for ln in lines]
    for ln in lines:
        if ln != '\n':
            vv = [string.strip() for string in ln.split('=', 1)]
            if len(vv) < 2:
                raise ValueError('Variables must be set by \'=\'')
            else:
                fvars[vv[0]] = convert_str(vv[1])

    for ivar in reqvars:
        if ivar not in fvars:
            raise KeyError('Missing input variable: {}'.format(ivar))

    return fvars


def cfg_update(defdict, fname, reqvars=[]):
    """Updates a dictionary of default variables with values from
    a configuration file."""
    if os.path.exists(fname):
        defdict.update(read_cfg(fname))
    elif reqvars != []:
        raise KeyError('Missing required input variables')


def get_fnames(matchex):
    """Returns a list of filenames based on a matching expression
    that may include wildcards."""
    fnames = []
    if isinstance(matchex, list):
        for reffn in matchex:
            fnames += glob(reffn)
    else:
        fnames += glob(matchex)

    return fnames


def read_dat(fname, labelrow=1, labelcol=0):
    """Reads an array of data from an input file."""
    for i in range(labelrow):
        labels = np.readline().split()
    return labels, np.loadtxt(fname, skiprows=labelrow)


def write_dat(fname, data, labels=None, charwid=10, decwid=4):
    """Writes an array of data to an output file."""
    with open(fname, 'w') as f:
        if labels is not None:
            f.write(''.join(['{:>{w}s}'.format(lbl, w=charwid) for
                             lbl in labels]) + '\n')
        for line in data:
            f.write(''.join(['{:{w}.{d}f}'.format(num, w=charwid, d=decwid) for
                             num in line]) + '\n')
