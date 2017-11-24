"""
The default dictionary values for different routines.

When a console routine it calls, the default dictionary is called from this
file. Alternatively, the command getdefault [routine_name] will create
an input file with the default set of inputs.
"""
import os

inpname = dict(
    branchget = 'branch.inp',
    contourplot = 'contour.inp',
    denplot = 'den.inp',
    histplot = 'hist.inp',
    pesplot = 'pes.inp',
    popplot = 'pop.inp'
               )

branchget = dict(
    states = [1, 0],
    ref_geoms = '*.xyz',
    invert = True,
    permute = None,
    test_geoms = '../seed.*/Spawn.[0-9]'
                 )

contourplot = dict()

denplot = dict()

histplot = dict()

trpesplot = dict()

pesplot = dict(
    states = None,
    infname = 'nrg.out',
    plotorder = None,
    gwid = 1,
    gsep = 1,
    lblrot = 90,
    maxener = None,
    minener = None,
    econv = 27.21138505,
    data_name = 'pes.dat',
    plot_name = 'pes.pdf'
               )

popplot = dict(
    states = [0, 1, 2],
    fms_time_increment = 100,
    time_conv = 0.02418884326505,
    tmin = 0.0,
    tmax = 1000.0,
    directory_stem = '../seed.*',
    trajectory_files = '../seed.*/TrajDump.*',
    amplitude_data_name = 'pop.dat',
    amplitude_plot_name = 'pop.pdf',
    fit_function = None,
    p0 = [10, 50],
    fit_data_name = 'pop.fit',
    fit_plot_name = None
               )

def write_default(routine):
    """Writes a default input file based on the routine name."""
    routine_dict = globals()[routine]
    fname = inpname[routine]

    # get existing input if possible
    if os.path.exists(fname):
        with open(fname, 'r') as f:
            orig_inp = f.readlines()
    else:
        orig_inp = []

    with open(fname, 'w') as f:
        # keep old input but comment it out
        for line in orig_inp:
            if '# default' not in line:
                f.write('#' + line)

        for key in routine_dict:
            val = routine_dict[key]
            if isinstance(val, list):
                if isinstance(val[0], list):
                    # 2D list, delimit with ',' and ';'
                    val = '; '.join([', '.join([str(i) for i in j]) for j in val])
                else:
                    # 1D list, delimit with ','
                    val = ', '.join([str(i) for i in val])
            else:
                # str() will handle ints, floats, bools and NoneType
                val = str(val)
            f.write('{:s} = {:s} # default\n'.format(key, val))
