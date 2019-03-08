"""
Module for reading and interpreting adiabatic populations, including
finding transferred population and fitting population curves.
"""
import os
import numpy as np
from scipy.optimize import curve_fit
from fmsinterpreter import fileio


def read_amps(fnames, aconv=1.):
    """Reads the amplitudes of each trajectory at given times and returns
    list of arrays of amplitudes, times and trajectory info."""
    amps = [np.empty(0) for fn in fnames]
    times = [np.empty(0) for fn in fnames]
    for i, fname in enumerate(fnames):
        times[i], amps[i] = fileio.read_dat(fname, skiprow=1, usecols=(0,-2)).T
        amps[i] *= aconv

    return amps, times


def read_ndat(fnames, nstates, aconv=1.):
    """Reads the amplitudes directly from N.dat and returns a list of
    amplitudes and times."""
    amps = [np.empty((nstates, 0)) for fn in fnames]
    times = [np.empty(0) for fn in fnames]
    for i, fname in enumerate(fnames):
        rawdat = fileio.read_dat(fname, skiprow=1).T
        times[i] = rawdat[0]
        amps[i] = np.empty((nstates, len(times[i])))
        for j in range(nstates):
            amps[i][j] = rawdat[j+1] * aconv

    return amps, times


def total_amps(fnames, times, nstates, aconv=1.):
    """Reads the trajectory amplitudes and sums them by their
    adiabatic state label for set time bins."""
    aseed, tseed = read_ndat(fnames, nstates, aconv=aconv)
    stamps = np.zeros((len(fnames), nstates, len(times)))
    for i in range(len(fnames)):
        tlocal = tseed[i]
        for j, t in enumerate(times):
            if t > tlocal[-1]:
                stamps[i,:,j:] += aseed[i][:,len(tlocal)-1:]
                break
            elif t > tlocal[0] - 1e-6:
                stamps[i,:,j] += aseed[i][:,np.argmin(np.abs(tlocal - t))]

    return stamps


def error_amps(times, stamps, nstates, nboot=1000, bthrsh=1e-3):
    """Calculates the amplitude errors using the bootstrap method.

    A random set of seeds are sampled until the bootstrap average
    converges to the true average or the maximum number of iterations
    is reached.
    """
    nseed = len(stamps)
    bootsamp = np.random.randint(nseed, size=(nboot, nseed))
    tavg = np.sum(stamps, axis=0)
    bavg = np.sum(stamps[bootsamp[0]], axis=0)
    bdel = np.zeros_like(stamps[0])
    for i in range(1, nboot):
        bstamp = np.sum(stamps[bootsamp[i]], axis=0)
        ei = bstamp - bavg
        bavg += ei / (i + 1)
        bdel += ei * (bstamp - bavg)
        if np.all(np.abs(bavg - tavg) < bthrsh):
            break
    print(i+1)
    return bavg, np.sqrt(bdel / i)


def read_tjinfo(fnames, noparent=False):
    """Reads the trajectory information for all files and returns
    an array."""
    tjinfo = np.empty((len(fnames), 5), dtype=int)
    for i, fname in enumerate(fnames):
        tjinfo[i] = get_spawn_info(fname, noparent=noparent)

    return tjinfo


def get_spawn_info(fname, first_lbl=1, noparent=False):
    """Reads the trajectory labels and states for parent and child of
    a given trajectory as well as the seed."""
    dsplit = fname.split('/')
    seed = int(dsplit[-2].split('.')[1])
    cid = int(dsplit[-1].split('.')[1])
    if cid > first_lbl and not noparent:
        try:
            spawnf = fname.replace('TrajDump', 'Spawn').replace('trajectory',
                                                                'spawn')
            with open(spawnf, 'r') as f:
                f.readline()
                comment = f.readline().split()
            cstate = int(comment[3]) - 1
            pid = int(comment[6])
            pstate = int(comment[9]) - 1
        except FileNotFoundError:
            # assume it is also a parent
            with open(fname, 'r') as f:
                f.readline()
                alldat = f.readline().split()
            cstate = int(alldat[-1].split('.')[0]) - 1
            pid = -1
            pstate = first_lbl
    else:
        with open(fname, 'r') as f:
            f.readline()
            alldat = f.readline().split()
        cstate = int(alldat[-1].split('.')[0]) - 1
        pid = -1
        pstate = first_lbl

    return np.array([seed, cid, cstate, pid, pstate], dtype=int)


def state_mask(traj_info, statei=None, statef=None):
    """Returns a boolean mask for states that have trajectory info
    matching given state indices and aren't initial trajectories."""
    ntraj = len(traj_info)
    mask0 = traj_info[:,3] >= 0
    if statei is not None:
        mask1 = traj_info[:,4] == statei
    else:
        mask1 = np.ones(ntraj, dtype=bool)
    if statef is not None:
        mask2 = traj_info[:,2] == statef
    else:
        mask2 = np.ones(ntraj, dtype=bool)
    mask12 = np.logical_and(mask1, mask2)
    return np.logical_and(mask0, mask12)


def integ_spawn_pop(amps, times, traj_info, statei=None, statef=None, maxt=400):
    """Returns the population transferred at each spawning event using
    the integration method.

    The change in population of the child is the negative of the change
    in population of the parent (assuming no effects from other
    trajectories). Thus, taking the negative product of the two
    derivatives gives a correlated function,
        (dp/dt)^2 = -(dp_c/dt)*(dp_p/dt),
    where the subscripts c and p correspond to child and parent. If
    negative values of (dp/dt)^2 are neglected, the total transferred
    population can be found by:
        p_trans = int[ dt * sgn(dp_c/dt) * sqrt((dp/dt)^2) ],
    where int represents integration over all time and sgn represents
    the sign function.
    """
    ntraj = len(amps)
    pops = np.zeros(ntraj)
    mask = state_mask(traj_info, statei=statei, statef=statef)
    irange = np.arange(ntraj)[mask]

    for i in irange:
        # get parent trajectory index and match times
        pmask = np.logical_and(traj_info[:,0] == traj_info[i,0],
                               traj_info[:,1] == traj_info[i,3])
        if not np.any(pmask):
            raise ValueError('Parent with index {:d}'.format(traj_info[i,3]) +
                             ' not found for trajectory ' +
                             '{:d} in seed {:d}.'.format(traj_info[i,1],
                                                         traj_info[i,0]))
        j = np.argmax(pmask)
        tshft = np.argmin(np.abs(times[j] - times[i][0]))
        lm = min(len(times[i]), len(times[j])-tshft, maxt)

        # get the child and parent populations
        pc = amps[i][:lm]
        pp = amps[j][tshft:lm+tshft]

        # take the derivatives
        dpc = pc[1:] - pc[:-1]
        dpp = pp[1:] - pp[:-1]

        # multiply the derivatives, take the sqrt and unwrap the sign
        dp2 = -dpc*dpp
        dp2[dp2 < 0] = 0
        dp = np.sign(dpc) * np.sqrt(dp2)

        # integrate from the inital child population
        pops[i] = pc[0] + np.sum(dp)

    return pops, mask


def thrsh_spawn_pop(amps, times, traj_info, inthrsh=5e-4, fithrsh=1e-4, nbuf=4):
    """Returns the population tranferred at each spawning event using
    the threshold method.

    A population event generally corresponds to a nearly stepwise
    change in the population, and thus a spike in its derivative. When
    dp/dt goes above a specified value, that is marked at the inital time
    and a final time is found where dp/dt falls below another given
    threshold.

    To correct for possible turning points, a buffer window
    can be chosen such that all points in the buffer must fall below the
    threshold before the final time is chosen.
    """
    ntraj = len(amps)
    pops = np.zeros(ntraj)

    for i in range(ntraj):
        if traj_info[i,3] >= 0:
            p = amps[i]
            t = times[i]

            # get the derivative, population range and number of time steps
            dp = (p[1:] - p[:-1]) / (t[1:] - t[:-1])
            rng = np.ptp(p)
            nt = len(dp)

            # find the time where dp exceeds inthrsh
            for iin in range(nt):
                if dp[iin] > inthrsh*rng:
                    break

            # find the time where dp falls below fithrsh for nbuf steps
            for iout in range(iin, nt-nbuf):
                comp = dp[iout:iout+nbuf]
                if max(comp) < fithrsh*rng:
                    iout += 1
                    break

            pops[i] = p[iout] - p[iin]

    return pops


def fit_function(func, times, decay, p0, tconv=1., err=None):
    """Fits amplitudes to a given exponential decay function.

    For now, this is just fitting the return to the ground state. With
    some sort of string parsing it could be changed to accept any
    state or combination of states (e.g. 'S1 + S2').
    """
    if err is not None:
        mask = err > 1e-5
        err = err[mask]
        times = times[mask]
        decay = decay[mask]
        abssig = True
    else:
        abssig = False
    t = times * tconv
    popt, pcov = curve_fit(globals()[func], t, decay, p0=p0, sigma=err,
                           absolute_sigma=abssig)
    perr = np.sqrt(np.diag(pcov))

    return popt, perr


def write_fit(func, popt, perr, outfname):
    """Writes fit information to an output file.

    This should be generalized to accept more than one set of fit values.
    """
    if func == 'exp':
        fitvals = ['t0', 'tau1']
    if func == 'expc':
        fitvals = ['t0', 'tau1', 'c']
    elif func == 'biexp':
        fitvals = ['t0', 'amp1', 'tau1', 'amp2', 'tau2']
    elif func == 'triexp':
        fitvals = ['t0', 'amp1', 'tau1', 'amp2', 'tau2', 'amp3', 'tau3']

    with open(outfname, 'w') as f:
        f.write('Curve  ')
        f.write(''.join(['{:>10s}'.format(fv) for fv in fitvals]) + '\n')
        f.write('1-S0   ')
        f.write(''.join(['{:10.4f}'.format(p) for p in popt]) + '\n')
        f.write('Error  ')
        f.write(''.join(['{:10.4f}'.format(p) for p in perr]) + '\n')


def exp(x, x0, b):
    """Returns an exponential function for curve fitting purposes."""
    return np.exp(-(x - x0) / b)


def expc(x, x0, b, c):
    """Returns an exponential function plus a constant for curve fitting
    purposes."""
    return np.exp(-(x - x0) / b) + c


def biexp(x, x0, a1, b1, a2, b2):
    """Returns a biexponential function for curve fitting purposes."""
    return a1 * np.exp(-(x - x0) / b1) + a2 * np.exp(-(x - x0) / b2)


def triexp(x, x0, a1, b1, a2, b2, a3, b3):
    """Returns a triexponential function for curve fitting purposes."""
    return (a1 * np.exp(-(x - x0) / b1) + a2 * np.exp(-(x - x0) / b2) +
            a3 * np.exp(-(x - x0) / b3))
