"""
FMS analysis routine for visualizing the adiabatic population as a
function of time and fitting the population to a mono-, bi- or
triexponential curve.
"""
import os
import numpy as np
from scipy.optimize import curve_fit
from fmsinterpreter import fileio


def read_amps(fnames, times, states, aconv=1.):
    """Reads the amplitudes of each trajectory at given times and returns
    arrays of amplitudes, child state and parent state."""
    amps = np.zeros((len(fnames), len(times)))
    tjinfo = np.empty((len(fnames), 5), dtype=int)

    for i, fname in enumerate(fnames):
        tlocal, alocal = fileio.read_dat(fname, skiprow=1, usecols=(0,-2)).T
        tjinfo[i] = get_spawn_info(fname)
        #with open(fname) as f:
        #    lines = [line.split() for line in f.readlines() if
        #             'Time' not in line and line != '\n']

        #alldata = np.array(lines, dtype=float)
        #tlocal = alldata[:,0]
        #st = states.index(int(alldata[0,-1]) - 1)
        #if st in states:
        for j, t in enumerate(times):
            if t > tlocal[-1]:
                amps[i, j:] += alocal[len(tlocal)-1]
                break
            elif t >= tlocal[0]:
                amps[i, j] += alocal[np.argmin(np.abs(tlocal - t))]

    return aconv*amps, tjinfo


def get_spawn_info(fname, first_lbl=1):
    """Reads the trajectory labels and states for parent and child of
    a given trajectory."""
    dsplit = fname.split('/')
    seed = int(dsplit[-2].split('.')[1])
    cid = int(dsplit[-1].split('.')[1])
    if cid > first_lbl:
        spawnf = fname.replace('TrajDump', 'Spawn').replace('trajectory',
                                                            'spawn')
        with open(spawnf, 'r') as f:
            f.readline()
            comment = f.readline().split()
        cstate = int(comment[3]) - 1
        pid = int(comment[6])
        pstate = int(comment[9])
    else:
        with open(fname, 'r') as f:
            f.readline()
            alldat = f.readline().split()
        cstate = int(alldat[-1].split('.')[0]) - 1
        pid = -1
        pstate = first_lbl

    return np.array([seed, cid, cstate, pid, pstate], dtype=int)


def total_amps(amps, states):
    """Sums a set of amplitudes from different trajectories by their
    adiabatic state label."""
    stamps = np.zeros((max(states)+1, amps.shape[1]))
    for i, s in enumerate(states):
        stamps[s] += amps[i]
    return stamps


def integ_spawn_pop(amps, traj_info):
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
    ti = np.argmax(amps > 1e-10, axis=1) - 1
    ntraj = len(amps)
    pops = np.zeros(ntraj)

    for i in range(ntraj):
        if traj_info[i,3] >= 0:
            # get the child and parent populations
            pc = amps[i,ti:ti+400]
            mask = np.logical_and(traj_info[:,0] == traj_info[i,0],
                                  traj_info[:,1] == traj_info[i,3])
            pp = amps[np.argmax(mask),ti:ti+400]

            # take the derivatives
            dpc = pc[1:] - pc[:-1]
            dpp = pp[1:] - pp[:-1]

            # multiply the derivatives, take the sqrt and unwrap the phase
            dp2 = -dpc*dpp
            dp2[dp2 < 0] = 0
            dp = np.sign(dpc) * np.sqrt(dp2)

            # integrate from the inital child population
            pops[i] = np.sum(dp)

    return pops


def thrsh_spawn_pop(amps, traj_info, inthrsh=5e-3, fithrsh=1e-3, nbuf=4):
    """Returns the population tranferred at each spawning event using
    the threshold method.

    A population event generally corresponds to a nearly stepwise
    change in the population, and thus a spike in its derivative. When
    dp/dt goes above a specified value, that is marked at the inital time
    and a final time is found where dp/dt falls below another given
    threshold. To correct for possible turning points, a buffer window
    can be chosen such that all points in the buffer must fall below the
    threshold before the final time is chosen.
    """
    ti = np.argmax(amps > 1e-10, axis=1) - 1
    ntraj = len(amps)
    pops = np.zeros(ntraj)

    for i in range(ntraj):
        if traj_info[i,3] >= 0:
            p = amps[i,ti:ti+400]

            # get the derivative, population range and number of time steps
            dp = p[1:] - p[:-1]
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
    if func == 'exp':
        fitvals = ['t0', 'tau1']
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


def biexp(x, x0, a1, b1, a2, b2):
    """Returns a biexponential function for curve fitting purposes."""
    return a1 * np.exp(-(x - x0) / b1) + a2 * np.exp(-(x - x0) / b2)


def triexp(x, x0, a1, b1, a2, b2, a3, b3):
    """Returns a triexponential function for curve fitting purposes."""
    return (a1 * np.exp(-(x - x0) / b1) + a2 * np.exp(-(x - x0) / b2) +
            a3 * np.exp(-(x - x0) / b3))
